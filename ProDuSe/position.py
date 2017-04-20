from scipy.stats import fisher_exact

# If not installed, or running in python2, this works fine
try:
    import nucleotide
except ImportError:
    # If installed and running in Python3
    from ProDuSe import nucleotide

import bisect
import sys
import itertools
import numpy


class ID:

    def __init__(self, ref, start):
        self.ref = ref
        self.start = start


class Tally:

    def __init__(self, list_of_bases):
        self.a = sum('A' == b for b in list_of_bases)
        self.c = sum('C' == b for b in list_of_bases)
        self.g = sum('G' == b for b in list_of_bases)
        self.t = sum('T' == b for b in list_of_bases)

    def __str__(self):
        return ':'.join([str(self.a), str(self.c), str(self.g), str(self.t)]);

    def get_count(self, base):
        if base == 'A':
            return self.a
        elif base == 'C':
            return self.c
        elif base == 'G':
            return self.g
        elif base == 'T':
            return self.t
        else:
            sys.exit("Base not available in Tally")

    def sum(self):
        return (self.a + self.c + self.g + self.t);


class Qname:

    def __init__(self, pysam_read):
        self.qname = pysam_read.qname
        self.mismatches = 0
        if pysam_read.has_tag("NM"):
            self.mismatches = pysam_read.get_tag("NM")
        split = self.qname.split(":")
        self.adapter_sequence = split[0]
        self.ref_start = split[1]
        self.pos_start = int(split[2])
        self.ref_end = split[3]
        self.pos_end = int(split[4])
        self.support = int(split[5])
        self.strand = split[6]
        self.duplex_id = int(split[7])
        self.is_paired = pysam_read.is_paired

        if pysam_read.is_reverse:
            self.mapstrand = "-"
        else:
            self.mapstrand = "+"
        #override strand using map information for FLASH pairs (for this to work, the negative FLASH reads need to be rev-complemented before aligning
        if not self.is_paired:
            self.strand = self.mapstrand
    def __eq__(self, other):
        return self.qname == other.qname

class Order:

    def __init__(self, chrom, start):
        self.chrom = chrom
        self.start = start

    def __lt__(self, other):
        if self.chrom < other.chrom:
            return True
        elif self.chrom == other.chrom and self.start < other.start:
            return True
        else:
            return False

    def lessthen(self, other, buffer):
        if self.chrom < other.chrom:
            return True
        elif self.chrom == other.chrom and self.start < (other.start - buffer):
            return True
        else:
            return False

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start

    def __hash__(self):
        return hash((self.chrom, self.start))

    def __str__(self):
        return ':'.join([str(self.chrom), str(self.start)])


class Pos:

    def __init__(self, base, qual, qname, order, minQual=3):
        self.order = order
        #convert extremely low quality bases to N
        if ord(qual)-33 < minQual:
            #print "low qual"
            #print ord(qual)-33
            self.base = "N"
        else:
            self.base = base
        self.qual = qual
        self.qname = qname

    def is_positive(self):
        #return self.qname.strand == "+"
        if self.is_paired:
            return self.qname.strand == "+"
        else:
            return self.qname.mapstrand == "+" #added by Ryan. Now that Flash is being used, the strand is consistent with the mapped strand for collapsed reads
    def __eq__(self, other):
        return self.qname.duplex_id == qname.other.duplex_id

    def __lt__(self, other):
        return self.qname.duplex_id < other.qname.duplex_id

    def __le__(self, other):
        return self.qname.duplex_id <= other.qname.duplex_id

    def __gt__(self, other):
        return self.qname.duplex_id > other.qname.duplex_id

    def __ge__(self, other):
        return self.qname.duplex_id >= other.qname.duplex_id


class Pos2(Pos):

    def __init__(self, pos):
        self.count = pos.count * -1
        self.base = pos.base
        self.adapter = pos.adapter

    def __lt__(self, other):
        return self.count < other.count

    def __eq__(self, other):
        return self.count == other.count

    def __gt__(self, other):
        return self.count > other.count

class PosCollection:


    def __init__(self, list_of_pos, order, chrom, start, ref):

        self.base_array = {"DPN":{}, "DPn":{}, "DpN":{}, "Dpn":{}, "SP":{}, "Sp":{}, "SN":{} ,"Sn":{},"StP":{},"Stp":{},"StN":{},"Stn":{},"StrBiasP":{}}
        for categ in self.base_array.keys():
            self.base_array[categ] = {"A":0, "C":0, "T":0,"G":0,"N":0}
        self.list_of_pos = list_of_pos
        self.order = order
        self.chrom = chrom
        self.start = start
        self.ref = nucleotide.makeCapital(ref)
        self.alt = []
        self.alt_reason = []
        self.variant_status = None
        self.bases = {}
        self.duplex_bases = []
        self.good_singleton_bases = []
        self.bad_singleton_bases = []
        self.good_duplex_bases = []
        self.bad_duplex_bases = []
        self.good_conflicting_bases = []
        self.bad_conflicting_bases = []
        self.good_conflicting_duplex_bases = []
        self.bad_conflicting_duplex_bases = []
        self.skipped_reads = 0 #number of consensus reads discarded at summary stage (calc_base_stats)
        self.pos_counts = {"A":0, "C":0, "T":0, "G":0, "N":0}
        self.neg_counts = {"A":0, "C":0, "T":0, "G":0, "N":0}
        self.pos_strand_counts = {"A":0, "C":0, "T":0, "G":0, "N":0}
        self.neg_strand_counts = {"A":0, "C":0, "T":0, "G":0, "N":0}


    def calc_base_stats(self, min_reads_per_uid, max_mismatch_per_aligned_read=5):
        #to be compatible with FLASH outputs, read strand information should only be considered from un-merged reads (i.e. mapped as pairs)
        for pos in self.list_of_pos:
            if pos.qname.duplex_id not in self.bases:
                self.bases[pos.qname.duplex_id] = []
            self.bases[pos.qname.duplex_id].append(pos)
        skipped = 0
        for duplex_id in self.bases:

            column = ""
            colnew = ""
            if len(self.bases[duplex_id]) == 1:

                passes_filter = self.bases[duplex_id][0].qname.support >= min_reads_per_uid
                too_many_mismatches = self.bases[duplex_id][0].qname.mismatches > max_mismatch_per_aligned_read
                is_positive = self.bases[duplex_id][0].qname.strand == "+"

                is_positive_strand = self.bases[duplex_id][0].qname.mapstrand == "+"
                is_paired = self.bases[duplex_id][0].qname.is_paired
                if too_many_mismatches:
                    skipped += 1
                    continue

                if passes_filter and is_positive:
                    column = "SP"

                elif passes_filter and not is_positive:
                    column = "SN"

                elif not passes_filter and is_positive:
                    column = "Sp"

                else:
                    column = "Sn"

                self.base_array[column][self.bases[duplex_id][0].base] += 1
                if passes_filter and is_positive_strand:
                    colnew = "StP"
                elif passes_filter and not is_positive_strand:

                    colnew = "StN"

                elif not passes_filter and is_positive_strand:
                    colnew = "Stp"

                else:
                    colnew = "Stn"

                self.base_array[colnew][self.bases[duplex_id][0].base] += 1
                if is_positive_strand:
                    self.pos_strand_counts[self.bases[duplex_id][0].base] += 1
                else:
                    self.neg_strand_counts[self.bases[duplex_id][0].base] += 1
                #the above few lines were added by Ryan to compensate for an issue that led to variants being called with almost exclusively collapsed reads mapped only to one strand
                if is_positive:
                    self.pos_counts[self.bases[duplex_id][0].base] += 1
                else:
                    self.neg_counts[self.bases[duplex_id][0].base] += 1

            elif len(self.bases[duplex_id]) == 2:
                too_many_mismatches1 = self.bases[duplex_id][0].qname.mismatches > max_mismatch_per_aligned_read
                too_many_mismatches2 = self.bases[duplex_id][1].qname.mismatches > max_mismatch_per_aligned_read
                is_positive_strand = self.bases[duplex_id][0].qname.mapstrand == "+"
                if too_many_mismatches1:
                    #print "skipping DS %s because %i mismatchess" % (self.bases[duplex_id][0].qname.qname,self.bases[duplex_id][0].qname.mismatches)
                    skipped += 1
                    continue
                if too_many_mismatches2:
                    #print "skipping DS %s because %i mismatchess" % (self.bases[duplex_id][1].qname.qname,self.bases[duplex_id][1].qname.mismatches)
                    skipped += 1
                    continue
                if self.bases[duplex_id][0].base == self.bases[duplex_id][1].base:

                    column = ""

                    pos_passes_filter = self.bases[duplex_id][0].qname.support >= min_reads_per_uid
                    neg_passes_filter = self.bases[duplex_id][1].qname.support >= min_reads_per_uid

                    if self.bases[duplex_id][0].qname.strand == "-":
                        pos_passes_filter = self.bases[duplex_id][1].qname.support >= min_reads_per_uid
                        neg_passes_filter = self.bases[duplex_id][0].qname.support >= min_reads_per_uid
                    if is_positive_strand:
                        colnew = "StP"
                        self.pos_strand_counts[self.bases[duplex_id][0].base] += 1
                    else:
                        self.neg_strand_counts[self.bases[duplex_id][0].base] += 1
                        colnew = "StN"  #count all duplex by the mapstrand (consider them all high-conf)

                    self.base_array[colnew][self.bases[duplex_id][0].base] += 1
                    if pos_passes_filter and neg_passes_filter:
                        #if self.bases[duplex_id][0].base == "C":
                            #print "DPN"
                            #print self.bases[duplex_id][0].qname.qname
                            #print self.bases[duplex_id][1].qname.qname
                        column = "DPN"

                    elif pos_passes_filter and not neg_passes_filter:
                        #if self.bases[duplex_id][0].base == "C":
                         #   print "DPn"
                          #  print self.bases[duplex_id][0].qname.qname
                           # print self.bases[duplex_id][1].qname.qname
                        column = "DPn"

                    elif not pos_passes_filter and neg_passes_filter:
                        #if self.bases[duplex_id][0].base == "C":
                         #   print "DpN"
                          #  print self.bases[duplex_id][0].qname.qname
                           # print self.bases[duplex_id][1].qname.qname
                        column = "DpN"

                    else:
                        #if self.bases[duplex_id][0].base == "C":
                         #   print "Dpn"
                          #  print self.bases[duplex_id][0].qname.qname
                           # print self.bases[duplex_id][1].qname.qname
                        column = "Dpn"

                    self.base_array[column][self.bases[duplex_id][0].base] += 1
                    self.pos_counts[self.bases[duplex_id][0].base] += 1
                    self.neg_counts[self.bases[duplex_id][0].base] += 1

                else:

                    is_positive = self.bases[duplex_id][0].qname.strand == "+"

                    if is_positive:
                        self.pos_counts[self.bases[duplex_id][0].base] += 1
                        self.neg_counts[self.bases[duplex_id][1].base] += 1
                    else:
                        self.neg_counts[self.bases[duplex_id][0].base] += 1
                        self.pos_counts[self.bases[duplex_id][1].base] += 1

                    if self.bases[duplex_id][0].qname.support >= min_reads_per_uid and self.bases[duplex_id][1].qname.support >= min_reads_per_uid:
                        self.good_conflicting_duplex_bases.append(self.bases[duplex_id][0].base + self.bases[duplex_id][1].base)

                    else:
                        self.bad_conflicting_duplex_bases.append(self.bases[duplex_id][0].base + self.bases[duplex_id][1].base)
        self.skipped_reads = skipped

    def position_stats_header(self, outFile):
        """
        Writes a header string coresponding to molecule abundances to the provided file

        Args:
            outFile: An open file object, which the header line is written to
        """
        outFile.write("##Alternate allele(s) molecule counts for each locus\n")
        outFile.write("#TotalMol\tDPN\tDPn\tDpN\tDpn\tSP\tSp\tSN\tSn\n")

    def position_stats(self, outFile):
        """
        Prints out alt molecule count information at this locus

        Args:
            outFile: An open file object, which the molecule counts are written to
        """
        pos_info = ""
        # Calculates overall molecule abundance at this position (representing depth)
        categs = ["DPN", "DPn", "DpN", "Dpn", "SP", "Sp", "SN", "Sn"]
        total_molecules = 0
        allBases = ["A", "C", "G", "T"]
        refBases = list(x for x in allBases if x not in self.alt)
        for categ in categs:
            counts = sum(self.base_array[categ][base] for base in self.alt)
            molecule_counts = counts + sum(self.base_array[categ][base] for base in refBases)
            # If there are no molecules of this type at this locus, that is informative
            # Set the VAF to 0
            if molecule_counts == 0:
                molecVAF = 0
            else:
                molecVAF = float(counts) / float(molecule_counts)
            pos_info += str(molecVAF) + "\t"
            total_molecules += molecule_counts

        pos_info = str(total_molecules) + "\t" + pos_info[:-1] + "\n"
        outFile.write(pos_info)

    def is_variant(self, min_alt_vaf, min_molecule_count, enforce_dual_strand, mutant_molecules):

        base_sum = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}
        alt_bases = []
        for categ in self.base_array.keys():
            for base in self.base_array[categ].keys():
                if not base == self.ref:
                    alt_bases.append(base)
                base_sum[base] += self.base_array[categ][base]

        alt_bases = numpy.unique(alt_bases)
        # print self.base_array

        for alt_base in alt_bases:
            if alt_base == "N":
                continue
            if self.base_array["DPN"][alt_base] >= 1:
                self.alt.append(alt_base)
                self.alt_reason.append("DPN")

            if (self.base_array["DPn"][alt_base] + self.base_array["DpN"][alt_base]) >= 2:
                self.alt.append(alt_base)
                self.alt_reason.append("DPn/DpN")

            if self.base_array["Dpn"][alt_base] >= 2:
                self.alt.append(alt_base)
                self.alt_reason.append("Dpn")
            pval = None
            if sum(base_sum.values()) >= min_molecule_count and (float(base_sum[alt_base]) / float(sum(base_sum.values()))) >= min_alt_vaf:
                pass_min = False
                pass_dual = True
                if enforce_dual_strand:
                    pass_dual = False

                    if self.pos_strand_counts[alt_base] > 0 and self.neg_strand_counts[alt_base] > 0:
                        pass_dual = True
                        #determine the p value for the bias of reads mapped to + vs - strand
                        pval = fisher_exact([[self.pos_strand_counts[alt_base], self.neg_strand_counts[alt_base]], [self.pos_strand_counts[self.ref], self.neg_strand_counts[self.ref]]])[1]  # doesn't work with FLASH
                        #print "pvalue %s %s %s %s" % (self.pos_strand_counts[alt_base],self.neg_strand_counts[alt_base],self.pos_strand_counts[self.ref],self.neg_strand_counts[self.ref])
                        #print pval.two_tail
                        #pval = pvalue(self.pos_counts[alt_base],self.neg_counts[alt_base],self.pos_counts[self.ref],self.neg_counts[self.ref])
                        self.base_array["StrBiasP"][alt_base] = pval

                    else:
                        pass_dual = False
                else:
                    pval = fisher_exact([[self.pos_strand_counts[alt_base],self.neg_strand_counts[alt_base]],[self.pos_strand_counts[self.ref],self.neg_strand_counts[self.ref]]])[1] #doesn't work with FLASH
                    #print "pvalue %s %s %s %s" % (self.pos_strand_counts[alt_base],self.neg_strand_counts[alt_base],self.pos_strand_counts[self.ref],self.neg_strand_counts[self.ref])
                    #print pval.two_tail
                    self.base_array["StrBiasP"][alt_base] = pval
                if self.pos_counts[alt_base] + self.neg_counts[alt_base] >= mutant_molecules:
                    pass_min = True

                if pass_min and pass_dual:
                    self.alt.append(alt_base)
                    self.alt_reason.append("VAF")

        if len(self.alt) == 0:
            return False

        else:
            return True

    def write_header(self, file_handler, contigs, reference_file):
        file_handler.write('##fileformat=VCFv4.2\n')
        file_handler.write('##reference=%s\n' % (reference_file))
        for contig, length in contigs.items():
            file_handler.write('##contig=<ID=%s,length=%s,>\n' % (contig, length))
        file_handler.write('##INFO=<ID=DPN,Number=R,Type=Integer,Description="Duplex Support with Strong Positive and Strong Negative Consensus">\n')
        file_handler.write('##INFO=<ID=DPn,Number=R,Type=Integer,Description="Duplex Support with Strong Positive and Weak Negative Consensus">\n')
        file_handler.write('##INFO=<ID=DpN,Number=R,Type=Integer,Description="Duplex Support with Weak Positive and Strong Negative Consensus">\n')
        file_handler.write('##INFO=<ID=Dpn,Number=R,Type=Integer,Description="Duplex Support with Weak Positive and Weak Negative Consensus">\n')
        file_handler.write('##INFO=<ID=SP,Number=R,Type=Integer,Description="Singleton Support with Strong Positive Consensus">\n')
        file_handler.write('##INFO=<ID=Sp,Number=R,Type=Integer,Description="Singleton Support with Weak Positive Consensus">\n')
        file_handler.write('##INFO=<ID=SN,Number=R,Type=Integer,Description="Singleton Support with Strong Negative Consensus">\n')
        file_handler.write('##INFO=<ID=Sn,Number=R,Type=Integer,Description="Singleton Support with Weak Negative Consensus">\n')
        file_handler.write('##INFO=<ID=MC,Number=R,Type=Integer,Description="Molecule Counts">\n')
        file_handler.write('##INFO=<ID=PC,Number=R,Type=Integer,Description="Molecule Count in Positive Strand Including Bases of Disagreeance\n')
        file_handler.write('##INFO=<ID=NC,Number=R,Type=Integer,Description="Molecule Count in Negative Strand Including Bases of Disagreeance\n')
        file_handler.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    def write_variant(self, file_handler):

        categs = ["DPN","DPn","DpN","Dpn","SP","Sp","SN","Sn","StP","Stp","StN","Stn","StrBiasP"]
        bases = ["A", "C", "G", "T"]
        info_column = ';'.join(['='.join([categ, ','.join([str(self.base_array[categ][base]) for base in bases])]) for categ in categs])

        molecule_counts = {"A":0, "C":0,"G":0,"T":0}
        for base in bases:
            for categ in categs:
                if categ in ("StP","Stp","StN","Stn","StrBiasP"):
                    continue
                molecule_counts[base] += self.base_array[categ][base]

        info_column = ';'.join([ info_column, '='.join(["MC", ','.join([ str(molecule_counts[base]) for base in bases])])])
        info_column = ';'.join([ info_column, '='.join(["PC", ','.join([ str(self.pos_counts[base]) for base in bases])])])
        info_column = ';'.join([ info_column, '='.join(["NC", ','.join([ str(self.neg_counts[base]) for base in bases])])])
        sr_detail = "SR=" + str(self.skipped_reads)
        info_column = ';'.join([ info_column, sr_detail])

        reason = "."
        if len(self.alt_reason) >= 0:
            reason = ','.join(self.alt_reason)

            file_handler.write('\t'.join([
                str(self.chrom),
                str(self.start + 1),
                ".",
                self.ref,
                ','.join(numpy.unique(self.alt)),
                ".",
                ".",
                info_column
                ]))
            file_handler.write("\n")

    def coords(self):

        return ''.join([ str(self.chrom), ':', str(self.start+1) ]);

    def coords2(self):

        return ''.join([ str(self.chrom), ':', str(self.start+1), '-', str(self.start+1) ]);


class PosCollectionCreate:

    def __init__(self, pysam_alignment_file, pysam_fasta_file, filter_overlapping_reads=True, target_bed=None, min_reads_per_uid=2, min_base_qual=3):
        self.pysam_alignment_file = pysam_alignment_file
        self.min_reads_per_uid = min_reads_per_uid
        self.pysam_fasta_file = pysam_fasta_file
        self.pos_collections = {}
        self.filter_overlapping_reads = filter_overlapping_reads
        self.qname_collections = {}
        self.base_buffer = 5000
        self.order = []
        self.read_processed = False
        self.qnames = {}
        self.first = True
        self.target_bed = target_bed
        self.pysam_alignment_generator = None
        self.min_base_qual = min_base_qual

        # READS WHICH OVERLAP TWO REGIONS MAY BE COUNTED TWICE (VERY VERY VERY VERY RARE, will fix later).
        if self.target_bed is not None:
            for region in self.target_bed.regions:
                if self.pysam_alignment_generator is None:
                    self.pysam_alignment_generator = self.pysam_alignment_file.fetch(region=region)
                else:
                    self.pysam_alignment_generator = itertools.chain(self.pysam_alignment_generator, self.pysam_alignment_file.fetch(region=region))

    def __iter__(self):
        return self.next()

    def __next__(self):
        return self.next()

    def return_pos(self):
        chrom = self.pysam_alignment_file.getrname(self.order[0].chrom)
        start = int(self.order[0].start)
        item = PosCollection(
            self.pos_collections[self.order[0]],
            self.order[0],
            chrom,
            start,
            self.pysam_fasta_file.fetch(chrom, start, start+1)
            )

        del self.pos_collections[self.order[0]]
        del self.qname_collections[self.order[0]]
        del self.order[0]
        return item

    def next(self):
        reads_visited = 0
        while True:

            # Get the next Sequence from the pysam object
            read = None

            try:
                if self.target_bed is None:
                    read = next(self.pysam_alignment_file)

                else:
                    read = next(self.pysam_alignment_generator)

            except StopIteration:

                break

            read_order = Order(read.reference_id, read.reference_start)

            qname = Qname(read)

            if qname.support < self.min_reads_per_uid:
                continue

            # For each base in the alignment, add to the collection structure
            for pairs in read.get_aligned_pairs(matches_only=True):

                if pairs[0] == None or pairs[1] == None:
                    continue

                #IGNORE SOFT CLIPPED BASES
                #print cigarstuff[pairs[0]]
                #print pairs
		#if cigarstuff[pairs[0]]:
                #    continue

                reads_visited+=1
                order = Order(read.reference_id, pairs[1])

                base = read.seq[pairs[0]]
                qual = read.qual[pairs[0]]

                current_pos = Pos(base, qual, qname, order, self.min_base_qual)

                if not order in self.pos_collections:
                    self.pos_collections[order] = []
                    self.qname_collections[order] = {}
                    bisect.insort(self.order, order)

                if not self.filter_overlapping_reads:
                    bisect.insort(self.pos_collections[order], current_pos)

                elif not read.qname in self.qname_collections[order]:
                    self.qname_collections[order][read.qname] = True
                    bisect.insort(self.pos_collections[order], current_pos)

                elif read.qname in self.qname_collections[order]:

                    left = bisect.bisect_left(self.pos_collections[order], current_pos)
                    right = bisect.bisect_right(self.pos_collections[order], current_pos);
                    index = None
                    for i in range(left, right):
                        if self.pos_collections[order][i].qname == current_pos.qname:
                            index = i

                    if self.first:
                        self.first = False

                    elif not self.first:
                        self.first = True
                        self.pos_collections[order][index] = current_pos

            while not len(self.order) == 0 and self.order[0].lessthen(read_order, self.base_buffer):

                yield self.return_pos();
            #print "visited %s reads" % reads_visited
        while not len(self.order) == 0:
            yield self.return_pos();

        raise StopIteration
