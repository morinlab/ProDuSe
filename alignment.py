import operator
import nucleotide
import collections
import fastq
import sys
import itertools

NUC_TO_INDEX = {
    'A':0,
    'C':1,
    'G':2,
    'T':3,
    'N':4
    }

INDEX_TO_NUC = ['A','C','G','T','N']

class ID:

    def __init__(self, ref_one, start, ref_two, end):
        self.ref_one = ref_one
        self.start = start
        self.ref_two = ref_two
        self.end = end

    def __str__(self):
        return ':'.join([
            str(self.ref_one),
            str(self.start),
            str(self.ref_two),
            str(self.end)
            ])

    def __eq__(self, other):
        if other == None:
            return False
        else:
            return self.ref_one == other.ref_one and self.start == other.start and self.ref_two == other.ref_two and self.end == other.end

    def __lt__(self, other):
        return self.ref_one < other.ref_one or (self.end < other.start and self.ref_one <= other.ref_one)

    def lessthan(self, bp_buffer, other):
        return self.ref_one < other.ref_one or (self.end + bp_buffer < other.start and self.ref_one <= other.ref_one)


class Read:

    def __init__(self, pysam_read):
        self.qname = pysam_read.qname
        self.seq = pysam_read.seq
        self.qual = pysam_read.qual
        self.qual_int = pysam_read.query_qualities
        self.ref_id = pysam_read.reference_id
        self.ref_len = pysam_read.reference_length
        self.misformed = False
        self.start = 0
        self.end = 0
        self.start_clip = 0
        self.end_clip = 0

        if pysam_read.cigartuples == None:
            self.misformed = True

        else :
            if len(pysam_read.cigartuples) > 1 and pysam_read.cigartuples[0][0] == 4:
                self.start_clip = pysam_read.cigartuples[0][1]

            if len(pysam_read.cigartuples) > 1 and pysam_read.cigartuples[-1][0] == 4:
                self.end_clip = pysam_read.cigartuples[-1][1]

            # Get the True start and end position of the sequence
            self.start = pysam_read.reference_start - self.start_clip 
            self.end = pysam_read.reference_end + self.end_clip
            
        # Determine if the read is a split alignment
        self.is_split = any(['SA' in x for x in pysam_read.tags])
                
        # Is this read the first in pair?
        self.is_read1 = pysam_read.is_read1


class Alignment:

    def __init__(self, pysam_read_one, pysam_read_two):
        self.r1 = Read(pysam_read_one)
        self.r2 = Read(pysam_read_two)
        self.adapter = self.r1.qname.split(":")[0]
        self.id = ID(self.r1.ref_id, self.r1.start, self.r2.ref_id, self.r2.end)

    def __str__(self):
        return "{0}\t{1}".format(
            self.r1.qname,
            self.id
            )
    # Return the start of read one
    def start(self):
        return self.r1.start

    # Return the end of read two
    def end(self):
        return self.r2.end

    def qname(self):
        return self.r1.qname

    # Check if either reads are split alignments
    def is_split(self):
        return self.r1.is_split or self.r2.is_split

    # Check if reads have different references
    def is_multi_ref(self):
        return not self.r1.ref_id == self.r2.ref_id

    # Check if the first read seen is read one, i.e. positive
    def is_positive(self):
        return self.r1.is_read1

    # Check if the Alignment is misformed
    def is_misformed(self):
        return self.r1.misformed or self.r2.misformed


class AlignmentCollection:

    def __init__(self, list_of_alignments, collection_type='', collection_id=''):
        self.data = list_of_alignments
        self.collection_type = collection_type
        self.collection_id = collection_id
        self.length = len(list_of_alignments)

    def __len__(self):
        return self.length


    def adapter_average_consensus(self, forward_fastq, reverse_fastq, max_mismatch=3, adapter_sequence="WSWSWGACT"):

        # We are trying to identify all of the adapter classes in the collection of Reads
        adapter = ''.join([adapter_sequence, adapter_sequence])
        adapter_to_adapter_class = {}
        adapter_class = []
        adapter_indexes = nucleotide.which_ambiguous(adapter)

        for key in collections.Counter([alignment.adapter for alignment in self.data]):

            # If the adapter_class list is empty, add the most frequent adapter as the first adapter class
            if len(adapter_class) == 0:
                adapter_class.append(key)
                adapter_to_adapter_class[key] = key

            # Otherwise we need to determine if this new adapter is its own class or already in a class
            elif not key in adapter_to_adapter_class:
                
                # Calculate the distance of the adapter to each known adapter classes
                distance_to_adapter_class = [ nucleotide.distance(x, key, adapter_indexes) for x in adapter_class ]

                # Find the minimum distance and its associated index in the list of distances to the adapter classes
                min_index, min_distance = min(enumerate(distance_to_adapter_class), key=operator.itemgetter(1))

                # If this value is larger then the max mismatch value, it must be the case that all values are larger
                # then the max mismatch value, and thus this key is a new adapter class
                if min_distance > max_mismatch:
                    adapter_class.append(key)
                    adapter_to_adapter_class[key] = key

                # Otherwise, update adapter to adapter class dict
                else:
                    adapter_to_adapter_class[key] = adapter_class[min_index]


        # Here we are grouping all alignments to their respected adapter class
        adapter_class_to_reads = {}
        for alignment in self.data:
            if adapter_to_adapter_class[alignment.adapter] in adapter_class_to_reads:
                adapter_class_to_reads[adapter_to_adapter_class[alignment.adapter]].append(alignment)
            else:
                adapter_class_to_reads[adapter_to_adapter_class[alignment.adapter]] = [alignment]

        # Here we are creating the consensus of each adapter class grouping
        for key in adapter_class:

            adapter = key

            value = adapter_class_to_reads[key]
            
            if self.collection_type == "-":
                adapter = key[(len(key)/2):(len(key))] + key[0:(len(key)/2)]

            id = ':'.join([
                ''.join(['@', adapter]),
                str(value[0].id),
                str(len(value)),
                str(self.collection_type)
                ])

            forward_consensus = consensus([x.r1 for x in value], 'first')
            reverse_consensus = consensus([x.r2 for x in value], 'second')

            # Write to Fastq
            if self.collection_type == "+":
                forward_fastq.next(fastq.Fastq(str(id), forward_consensus[0], '+', forward_consensus[1]))
                reverse_fastq.next(fastq.Fastq(str(id), nucleotide.reverseComplement(reverse_consensus[0]), '+', reverse_consensus[1][::-1]))

            if self.collection_type == "-":
                forward_fastq.next(fastq.Fastq(str(id), nucleotide.reverseComplement(reverse_consensus[0]), '+', reverse_consensus[1][::-1]))
                reverse_fastq.next(fastq.Fastq(str(id), forward_consensus[0], '+', forward_consensus[1]))             


def index_max(values):
    return max(xrange(len(values)),key=values.__getitem__)

def consensus(list_of_reads, strand):

    length = min([len(x.seq) for x in list_of_reads])
    seq = [None] * length
    qual = [None] * length

    for i in range(length):
        counts = [0] * 5
        quals = [0] * 5
        for read in list_of_reads:
            if strand == 'first':
                counts[NUC_TO_INDEX[read.seq[i]]] += 1
                if quals[NUC_TO_INDEX[read.seq[i]]] < read.qual_int[i]:
                    quals[NUC_TO_INDEX[read.seq[i]]] = read.qual_int[i]
            elif strand == 'second':
                counts[NUC_TO_INDEX[read.seq[::-1][i]]] += 1
                if quals[NUC_TO_INDEX[read.seq[::-1][i]]] < read.qual_int[i]:
                    quals[NUC_TO_INDEX[read.seq[::-1][i]]] = read.qual_int[i]
        seq[i] = INDEX_TO_NUC[index_max(counts)]
        qual[i] = min(quals[index_max(counts)] + max(counts) + 32, 72)

    if strand == 'first':
        return [''.join(seq), ''.join([str(unichr(x)) for x in qual])]
    else:
        return [''.join(seq)[::-1], ''.join([str(unichr(x)) for x in qual])[::-1]]

class AlignmentCollectionCreate:

    def __init__(self, pysam_alignment_file, sense_matters=True, verbose=None):
        self.pysam_alignment_file = pysam_alignment_file
        self.sense_matters = sense_matters
        self.qname_to_read = {}
        self.misformed = {}
        self.split = {} 
        self.multi_ref = {}
        self.tmp_collections = {}
        self.next_id = []
        self.id = None
        self.base_buffer = 5000

    def __iter__(self):
        return self.next()

    def __next__(self):
        return self.next()

    def next(self):

        while True:
            
            try:
                read = self.pysam_alignment_file.next()
            
            except StopIteration:

                while len(self.next_id) >= 0:

                    list_of_alignments = []
                    list_of_positive_alignments = []
                    list_of_negative_alignments = []

                    collection_id = str(self.id)

                    # for each alignment in the collection, determine the sense
                    for qname in self.tmp_collections[collection_id]:

                        tmp_align = self.tmp_collections[collection_id][qname]

                        # Add that alignment to their respected sense list
                        if not self.sense_matters:
                            list_of_alignments.append(tmp_align)

                        elif tmp_align.is_positive():
                            list_of_positive_alignments.append(tmp_align)

                        else:
                            list_of_negative_alignments.append(tmp_align)

                    del self.tmp_collections[collection_id]

                    if not list_of_alignments == []:
                        yield AlignmentCollection(list_of_alignments, collection_type="?", collection_id=collection_id)

                    if not list_of_positive_alignments == []:
                        yield AlignmentCollection(list_of_positive_alignments, collection_type='+', collection_id=collection_id)

                    if not list_of_negative_alignments == []:
                        yield AlignmentCollection(list_of_negative_alignments, collection_type='-', collection_id=collection_id)


                    if len(self.next_id) == 0:

                        raise StopIteration

                    else:

                        self.id = self.next_id.pop(0)

            # Have we seen this read before?
            if not read.qname in self.qname_to_read:
                
                # If we haven't we should add this to the "qname to read" hash
                self.qname_to_read[read.qname] = read
                continue

            # Form and Alignment on the read and mate
            align = Alignment(self.qname_to_read[read.qname], read)

            # Delete the read from the qname to read hash
            del self.qname_to_read[read.qname]

            # Create the alignment ID
            id = align.id

            # Ignore misformed Alignments
            if align.is_misformed():
                self.misformed[str(id)] = align
                continue

            # Ignore Alignments that contain a Split Alignment tag
            if align.is_split():
                self.split[str(id)] = align
                continue

            # Ignore Alignments that map to different chromosomes
            if align.is_multi_ref():
                self.multi_ref[str(id)] = align
                continue

            # If this ID is new, push it to process
            if not str(id) in self.tmp_collections:
                self.next_id.append(id)

            # Process the first alignment if never done before
            if self.id == None:
                self.id = self.next_id.pop(0)
            
            # Enter this loop where we will begin processing alignments
            while self.id.lessthan(self.base_buffer, id):

                list_of_alignments = []
                list_of_positive_alignments = []
                list_of_negative_alignments = []

                collection_id = str(self.id)

                # for each alignment in the collection, determine the sense
                print 'hello'
                for qname in self.tmp_collections[collection_id]:

                    tmp_align = self.tmp_collections[collection_id][qname]

                    # Add that alignment to their respected sense list
                    if not self.sense_matters:
                        list_of_alignments.append(tmp_align)

                    elif tmp_align.is_positive():
                        list_of_positive_alignments.append(tmp_align)

                    else:
                        list_of_negative_alignments.append(tmp_align)

                del self.tmp_collections[collection_id]

                if not list_of_alignments == []:
                    yield AlignmentCollection(list_of_alignments, collection_type="?", collection_id=collection_id)

                if not list_of_positive_alignments == []:
                    yield AlignmentCollection(list_of_positive_alignments, collection_type='+', collection_id=collection_id)

                if not list_of_negative_alignments == []:
                    yield AlignmentCollection(list_of_negative_alignments, collection_type='-', collection_id=collection_id)

                
                self.id = self.next_id.pop(0)

            if not str(id) in self.tmp_collections:
                self.tmp_collections[str(id)] = { align.qname : align }
                continue

            else:
                self.tmp_collections[str(id)][align.qname] = align
                continue