import pysam
import bisect
import nucleotide
import collections
import sys
import itertools
from collections import Counter
import bed


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
        return self.qname.strand == "+"

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
        self.list_of_pos = list_of_pos
        self.order = order
        self.chrom = chrom
        self.start = start
        self.ref = nucleotide.makeCapital(ref)
        self.alt = None
        self.molecule_bases = []
        self.pos_coords = {}
        self.neg_coords = {}
        self.pos_bases = []
        self.neg_bases = []
        self.true_bases = []
        self.bases = {}
        self.duplex_bases = []
        self.good_singleton_bases = [] 
        self.bad_singleton_bases = []
        self.good_duplex_bases = []
        self.bad_duplex_bases = []
        self.good_conflicting_bases = []
        self.bad_conflicting_bases = []

    def calc_base_stats(self, min_reads_per_uid):

        for pos in self.list_of_pos:
            if not pos.qname.duplex_id in self.bases:
                self.bases[pos.qname.duplex_id] = []
            self.bases[pos.qname.duplex_id].append(pos)

        for duplex_id in self.bases:

            if len(self.bases[duplex_id]) == 1:

                if self.bases[duplex_id][0].qname.support >= min_reads_per_uid:
                    self.good_singleton_bases.append(self.bases[duplex_id][0].base)

                else:
                    self.bad_singleton_bases.append(self.bases[duplex_id][0].base)

            elif len(self.bases[duplex_id]) == 2:
                if self.bases[duplex_id][0].base == self.bases[duplex_id][1].base:
                    
                    if self.bases[duplex_id][0].qname.support >= min_reads_per_uid and \
                       self.bases[duplex_id][1].qname.support >= min_reads_per_uid:
                        self.good_duplex_bases.append(self.bases[duplex_id][0].base)
                    
                    else:
                        self.bad_duplex_bases.append(self.bases[duplex_id][0].base)

                else:

                    if self.bases[duplex_id][0].qname.support >= min_reads_per_uid and \
                       self.bases[duplex_id][1].qname.support >= min_reads_per_uid:
                        self.good_conflicting_duplex_bases.append(self.bases[duplex_id][0].base + self.bases[duplex_id][1].base)

                    else:
                        self.bad_conflicting_duplex_bases.append(self.bases[duplex_id][0].base + self.bases[duplex_id][1].base)


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


    def __str__(self):

        return '\t'.join([
            str(self.chrom),
            str(self.start + 1),
            str(Tally(self.good_duplex_bases)),
            str(Tally(self.good_singleton_bases)),
            ':'.join(self.good_conflicting_bases),
            ]);

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
        self.qnames = {};
        self.first = True;
        self.target_bed = target_bed;
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

            # For each base in the alignment, add to the collection structure
            for pairs in read.get_aligned_pairs():

                if pairs[0] == None or pairs[1] == None:
                    continue

                order = Order(read.reference_id, pairs[1])

                base = read.seq[pairs[0]]
                qual = read.qual[pairs[0]]

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

        while not len(self.order) == 0:
            yield self.return_pos();

        raise StopIteration
