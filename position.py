import pysam
import bisect
import nucleotide
import collections
import fasta
from collections import Counter

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

    def __init__(self, base, qual, qname):
        self.base = base
        self.qual = qual
        split = qname.split(':')
        self.adapter = split[0]
        self.start = int(split[2])
        self.end = int(split[4])
        self.sense = split[6]

    def __eq__(self, other):
        if self.start == other.start and self.end == other.end and self.base == other.base and self.sense != other.sense:
            return True
        else:
            return False

    def __lt__(self, other):
        if self.start < other.start:
            return True
        elif self.start == other.start:
            if self.end < other.end:
                return True
        else:
            return False

    def __le__(self, other):
        return self < other or self == other

    def __gt__(self, other):
        return not self < other and not self == other

    def __ge__(self, other):
        return not self < other


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


    def which_equal(self):
        which = []
        if len(self.list_of_pos) == 1:
            return which

        else:
            for i in range(len(self.list_of_pos[:-1])):
                if (self.list_of_pos[i] == self.list_of_pos[i+1]):
                    which.append(i);
            return which    


    def collect_true_bases(self, which, adapter_sequence, max_mismatch):

        for i in which:
            indexes = nucleotide.which_ambiguous(adapter_sequence)
            if nucleotide.distance(self.list_of_pos[i].adapter, self.list_of_pos[i+1].adapter, indexes) <= max_mismatch :
                self.true_bases.append(self.list_of_pos[i].base);


    def is_variant(self, adapter_sequence, max_mismatch):

        self.neg_bases = [ pos.base if pos.sense == "-" else None for pos in self.list_of_pos ];
        self.pos_bases = [ pos.base if pos.sense == "+" else None for pos in self.list_of_pos ];
        self.all_bases = self.neg_bases + self.pos_bases;
        which = self.which_equal();

        if not len(which) == 0: 
            self.collect_true_bases(which, adapter_sequence, max_mismatch);
            if len(self.true_bases) > 0:
                if not all(x == self.true_bases[0] for x in self.true_bases):
                    tmp = Counter(self.true_bases);
                    most = tmp.most_common(2);
                    if most[0][0] == self.ref:
                        self.alt = most[1][0];
                    else:
                        self.alt = most[0][0];
                    return True;
        else:
            return False


    def __str__(self):

        return '\t'.join([
            str(self.chrom),
            str(self.start),
            str(self.start),
            '/'.join([
                self.ref,
                self.alt
                ]),
            "1",
            str(Tally(self.true_bases)),
            str(Tally(self.pos_bases)),
            str(Tally(self.neg_bases))
            ]);

class PosCollectionCreate:

    def __init__(self, pysam_alignment_file, pysam_fasta_file):
        self.pysam_alignment_file = pysam_alignment_file
        self.pysam_fasta_file = pysam_fasta_file
        self.pos_collections = {}
        self.base_buffer = 5000
        self.order = []
        self.read_processed = False

    def __iter__(self):
        return self.next()

    def __next__(self):
        return self.next()

    def next(self):

        while True:
           
            # Get the next Sequence from the pysam object
            read = self.pysam_alignment_file.next()
            read_order = Order(read.reference_id, read.reference_start)

    
            # For each base in the alignment, add to the collection structure
            for pairs in read.get_aligned_pairs():

                if pairs[0] == None or pairs[1] == None:
                    continue

                order = Order(read.reference_id, pairs[1])
                
                base = read.seq[pairs[0]]
                qual = read.qual[pairs[0]]
                
                current_pos = Pos(base, qual, read.qname)

                if not order in self.pos_collections:
                    self.pos_collections[order] = []
                    bisect.insort(self.order, order)

                bisect.insort(self.pos_collections[order], current_pos)


            # Blah Blah
            while not len(self.order) == 0 and self.order[0].lessthen(read_order, self.base_buffer):

                chrom = self.pysam_alignment_file.getrname( self.order[0].chrom );
                start = int(self.order[0].start)
                yield PosCollection(
                    self.pos_collections[self.order[0]], 
                    self.order[0], 
                    chrom,
                    start,
                    self.pysam_fasta_file.fetch(chrom, start, start+1)
                    )

                del self.pos_collections[self.order[0]]
                del self.order[0]