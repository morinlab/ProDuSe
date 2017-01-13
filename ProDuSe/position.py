from fisher import pvalue
import pysam
import bisect
import nucleotide
import collections
import sys
import itertools
from collections import Counter
import bed

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

    def __init__(self, base, qual, qname, order):
        self.order = order
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
        #self.base_array = pandas.DataFrame(
        #    numpy.zeros((5, 8)),
        #    index = ["A", "C", "G", "T", "N"],
        #    columns = ["DPN", "DPn", "DpN", "Dpn", "SP", "Sp", "SN", "Sn"]
        #    )
        self.list_of_pos = list_of_pos
        self.order = order
        self.chrom = chrom
        self.start = start
        self.ref = nucleotide.makeCapital(ref)
        self.alt = []
        self.alt_reason = []
        # self.molecule_bases = []
        # self.pos_coords = {}
        # self.neg_coords = {}
        # self.pos_bases = []
        # self.neg_bases = []
        # self.true_bases = []
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
        self.pos_counts = {"A":0, "C":0, "T":0, "G":0, "N":0}
        self.neg_counts = {"A":0, "C":0, "T":0, "G":0, "N":0}
        self.pos_strand_counts = {"A":0, "C":0, "T":0, "G":0, "N":0}
        self.neg_strand_counts = {"A":0, "C":0, "T":0, "G":0, "N":0}
        
    def calc_base_stats(self, min_reads_per_uid):
        #to be compatible with FLASH outputs, read strand information should only be considered from un-merged reads (i.e. mapped as pairs)
        for pos in self.list_of_pos:
            if not pos.qname.duplex_id in self.bases:
                self.bases[pos.qname.duplex_id] = []
            self.bases[pos.qname.duplex_id].append(pos)

        for duplex_id in self.bases:

            column = ""
            colnew = ""
            if len(self.bases[duplex_id]) == 1:
                
                passes_filter = self.bases[duplex_id][0].qname.support >= min_reads_per_uid
                is_positive = self.bases[duplex_id][0].qname.strand == "+"
                
                is_positive_strand = self.bases[duplex_id][0].qname.mapstrand == "+"
                is_paired = self.bases[duplex_id][0].qname.is_paired
                #print "pairing for %s is %s" % (self.bases[duplex_id][0].qname,is_paired)
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

                if self.bases[duplex_id][0].base == self.bases[duplex_id][1].base:

                    column = ""

                    pos_passes_filter = self.bases[duplex_id][0].qname.support >= min_reads_per_uid
                    neg_passes_filter = self.bases[duplex_id][1].qname.support >= min_reads_per_uid

                    if self.bases[duplex_id][0].qname.strand == "-":
                        pos_passes_filter = self.bases[duplex_id][1].qname.support >= min_reads_per_uid
                        neg_passes_filter = self.bases[duplex_id][0].qname.support >= min_reads_per_uid
                    
                    if pos_passes_filter and neg_passes_filter:
                        column = "DPN"

                    elif pos_passes_filter and not neg_passes_filter:
                        column = "DPn"

                    elif not pos_passes_filter and neg_passes_filter:
                        column = "DpN"
                    
                    else:
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

    def is_variant(self, min_alt_vaf, min_molecule_count, enforce_dual_strand, mutant_molecules):
	    
        
        base_sum = {"A":0,"T":0,"C":0,"G":0,"N":0}
        alt_bases = []
	for categ in self.base_array.keys():
            for base in self.base_array[categ].keys():
                if not base == self.ref:
                    alt_bases.append(base)
                base_sum[base] += self.base_array[categ][base]

        alt_bases = numpy.unique(alt_bases) 
        #print self.base_array
        
        for alt_base in alt_bases:
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
                    if self.pos_strand_counts[alt_base] >0 and self.neg_strand_counts[alt_base]>0:
                        pass_dual = True
                        #determine the p value for the bias of reads mapped to + vs - strand
                        #pval = pvalue(self.pos_strand_counts[alt_base],self.neg_strand_counts[alt_base],self.pos_strand_counts[self.ref],self.neg_strand_counts[self.ref])
                        pval = pvalue(self.pos_counts[alt_base],self.neg_counts[alt_base],self.pos_counts[self.ref],self.neg_counts[self.ref])
                        self.base_array["StrBiasP"][alt_base] = pval.two_tail
                        
                    else:
                        pass_dual = False
                if self.pos_counts[alt_base] + self.neg_counts[alt_base] >= mutant_molecules:
                    pass_min = True
                else:
                    print "Skipping this because not enough mutant molecules %s and %s" % (self.pos_counts[alt_base],self.neg_counts[alt_base])
                if pass_min and pass_dual:
                    self.alt.append(alt_base)
                    self.alt_reason.append("VAF")

        if len(self.alt) == 0:
            return False

        else:
            return True
                



        # def is_variant(self, alt_base_count_threshold, strand_bias_threshold, variant_allele_fraction_threshold):

        # self.base_array["DPN"]

        # self.pos_bases = [ pos.base for pos in self.list_of_pos if pos.is_positive() ];
        # self.neg_bases = [ pos.base for pos in self.list_of_pos if not pos.is_positive() ];
        # print self.list_of_pos[0].order 
        # print self.pos_bases
        # print self.neg_bases
        # return (False, [])
        # self.neg_tally = Tally(self.neg_bases);
        # self.pos_tally = Tally(self.pos_bases);
        # self.all_bases = [ pos.base for pos in self.list_of_pos ];
        # self.norm_pos_tally = Tally(list(itertools.chain.from_iterable( [ [pos.base] * pos.count if pos.is_positive else [] for pos in self.list_of_pos ] )))
        # self.norm_neg_tally = Tally(list(itertools.chain.from_iterable( [ [pos.base] * pos.count if not pos.is_positive else [] for pos in self.list_of_pos ] )))
        # self.stuff = self.find_true_bases(adapter_indexes, max_mismatch);
        # self.true_bases = self.stuff[0]
        # self.wrong_bases = self.stuff[1]
        # self.is_variant = True;
        # self.reason = 'VAR';
        # self.alt = '.';

        # if not len(self.true_bases) == 0: 
        #     if not all(self.ref == x for x in self.true_bases):
        #         tmp = Counter(self.true_bases);
        #         most = tmp.most_common(2);

        #         if most[0][0] == self.ref:
        #             self.alt = most[1][0];
        #         else:
        #             self.alt = most[0][0];

        #         self.alt_count = sum([ True if x == self.alt else False for x in self.true_bases]);
        #         self.ref_count = sum([ True if x == self.ref else False for x in self.true_bases]);

        #         if not self.neg_tally.get_count(self.alt) >= alt_base_count_threshold:
        #             self.is_variant = False;
        #             self.reason = 'ABC';
    
        #         elif not self.pos_tally.get_count(self.alt) >= alt_base_count_threshold:
        #             self.is_variant = False;
        #             self.reason = 'ABC';

        #         elif not (float( min( self.neg_tally.sum(), self.pos_tally.sum() ) ) / float( self.neg_tally.sum() + self.pos_tally.sum() ) ) >= strand_bias_threshold:
        #             self.is_variant = False;
        #             self.reason = 'SB';

        #         elif not (float(self.alt_count) / float(self.ref_count + self.alt_count)) >= variant_allele_fraction_threshold:
        #             self.is_variant = False;
        #             self.reason = 'VAF';

        #         return (self.is_variant, self.wrong_bases);

        #     else:
        #         self.is_variant = False;
        #         self.reason = 'ABC';

        # else:
        #     self.is_variant = False;
        #     self.reason = 'COV';

        # return (self.is_variant, self.wrong_bases)



    def write_header(self, file_handler):
        file_handler.write('##fileformat=VCFv4.2\n')
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
	info_column = ';'.join(['='.join([categ, ','.join([ str(self.base_array[categ][base]) for base in bases])]) for categ in categs])
        
        molecule_counts = {"A":0, "C":0,"G":0,"T":0}	
        for base in bases:
            for categ in categs:
                if categ in ("StP","Stp","StN","Stn","StrBiasP"):
                    continue
                molecule_counts[base] += self.base_array[categ][base]
        
        info_column = ';'.join([ info_column, '='.join(["MC", ','.join([ str(molecule_counts[base]) for base in bases])])])
        info_column = ';'.join([ info_column, '='.join(["PC", ','.join([ str(self.pos_counts[base]) for base in bases])])])
        info_column = ';'.join([ info_column, '='.join(["NC", ','.join([ str(self.neg_counts[base]) for base in bases])])])

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

    def __init__(self, pysam_alignment_file, pysam_fasta_file, filter_overlapping_reads = True, target_bed = None, min_reads_per_uid = 2):
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

        ### READS WHICH OVERLAP TWO REGIONS MAY BE COUNTED TWICE (VERY VERY VERY VERY RARE, will fix later).
        if not self.target_bed == None:
            for region in self.target_bed.regions:
                if self.pysam_alignment_generator == None:
                    self.pysam_alignment_generator = self.pysam_alignment_file.fetch(region=region);
                else:
                    self.pysam_alignment_generator = itertools.chain(self.pysam_alignment_generator, self.pysam_alignment_file.fetch(region=region));
        
    def __iter__(self):
        return self.next()

    def __next__(self):
        return self.next()

    def return_pos(self):
        chrom = self.pysam_alignment_file.getrname( self.order[0].chrom );
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
        return item;

    def next(self):
        reads_visited = 0
        counts = 0
        while True:

            # Get the next Sequence from the pysam object
            read = None;

            try:
                if self.target_bed == None:
                    read = self.pysam_alignment_file.next()
                    
                else:
                    read = self.pysam_alignment_generator.next()

            except StopIteration:

                break;

            read_order = Order(read.reference_id, read.reference_start)

            qname = Qname(read)

            if qname.support < self.min_reads_per_uid:
                  continue

            #For each base in the alignment, add to the collection structure
            #cigarstuff = list(itertools.chain.from_iterable(numpy.repeat(val[0],val[1]) for val in read.cigartuples))
            #print "support for %s is %s, >= %s" % (qname,qname.support,self.min_reads_per_uid)
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
                

		#if pairs[1] == 148508727:
		#    if base == "C":
                #        print read.qname

                current_pos = Pos(base, qual, qname, order)

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
                            index = i;

                    if self.first:
                        self.first = False;

                    elif not self.first:
                        self.first = True;
                        self.pos_collections[order][index] = current_pos;

            while not len(self.order) == 0 and self.order[0].lessthen(read_order, self.base_buffer):

                yield self.return_pos();
            #print "visited %s reads" % reads_visited
        while not len(self.order) == 0:
            yield self.return_pos();

        raise StopIteration
