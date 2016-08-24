import operator
import nucleotide
import collections
import fastq
import sys
import itertools

NUC_TO_INDEX = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3,
    'N': 4
    }

INDEX_TO_NUC = ['A', 'C', 'G', 'T', 'N']


class ID:

    def __init__(self, ref_one, start, ref_two, end, ref_name_one, ref_name_two):
        self.ref_one = ref_one
        self.start = start
        self.ref_two = ref_two
        self.end = end
        self.ref_name_one = ref_name_one
        self.ref_name_two = ref_name_two

    def __str__(self):
        return ':'.join([
            str(self.ref_name_one),
            str(self.start),
            str(self.ref_name_two),
            str(self.end)
            ])

    def __hash__(self):
        return hash((self.ref_name_one, self.start, self.ref_name_two, self.end))

    def __eq__(self, other):
        if other == None:
            return False
        else:
            return self.ref_one == other.ref_one and self.start == other.start and self.ref_two == other.ref_two and self.end == other.end

    def __lt__(self, other):
        return self.ref_one < other.ref_one or (
            self.end < other.start and self.ref_one <= other.ref_one)

    def lessthan(self, bp_buffer, other):
        return self.ref_one < other.ref_one or (
            self.end + bp_buffer < other.start and self.ref_one <= other.ref_one)


class Read:

    def __init__(self, pysam_read):
        self.qname = pysam_read.qname
        self.seq = pysam_read.seq
        self.qual = pysam_read.qual
        self.qual_int = pysam_read.query_qualities
        self.ref_id = pysam_read.reference_id
        self.ref_name = pysam_read.reference_name
        self.ref_len = pysam_read.reference_length
        self.misformed = False
        self.start = 0
        self.end = 0
        self.start_clip = 0
        self.end_clip = 0

        if pysam_read.cigartuples == None:
            self.misformed = True

        else:
            if len(pysam_read.cigartuples) > 1 and pysam_read.cigartuples[
                   0][0] == 4:
                self.start_clip = pysam_read.cigartuples[0][1]

            if len(
                pysam_read.cigartuples) > 1 and pysam_read.cigartuples[-1][0] == 4:
                self.end_clip = pysam_read.cigartuples[-1][1]

            # Get the True start and end position of the sequence
            self.start = pysam_read.reference_start - self.start_clip
            self.end = pysam_read.reference_end + self.end_clip

        # Determine if the read is a split alignment
        self.is_split = any(['SA' in x for x in pysam_read.tags])

        # Is this read the first in pair?
        self.is_read1 = pysam_read.is_read1
        self.is_read2 = pysam_read.is_read2
        self.is_reverse = pysam_read.is_reverse


class Alignment:

    def __init__(self, pysam_read_one, pysam_read_two):
        self.r1 = Read(pysam_read_one)
        self.r2 = Read(pysam_read_two)
        self.adapter = self.r1.qname.split(":")[0]
        self.qname = ':'.join(self.r1.qname.split(":")[1:])
        self.id = ID(self.r1.ref_id, min(self.r1.start, self.r2.start),
                     self.r2.ref_id, max(self.r1.end, self.r2.end),
                     self.r1.ref_name, self.r2.ref_name)

    def __str__(self):
        return "{0}\t{1}".format(
            self.r1.qname,
            self.id
            )
    # Return the start of read one

    def start(self):
        return min(self.r1.start, self.r2.start)

    # Return the end of read two
    def end(self):
        return max(self.r1.end, self.r2.end)

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
        if self.r1.is_read1:
            if not self.r1.is_reverse:
                return True

        if self.r1.is_read2:
            if self.r1.is_reverse:
                return True

        else:
            return False

    # Check if the Alignment is misformed
    def is_misformed(self):
        return self.r1.misformed or self.r2.misformed

duplex_ids = {"+": {}, "-": {}}
duplex_stats = collections.OrderedDict()
duplex_uid = 0
prev_id = None

def index_max(values):
    return max(xrange(len(values)), key=values.__getitem__)


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
        return [''.join(seq)[::-1], ''.join([str(unichr(x))
                        for x in qual])[::-1]]

def adapter_flip(adapter):
    return adapter[(len(adapter) / 2):] + adapter[:(len(adapter) / 2)]


class AlignmentCollection:


    def __init__(
        self, list_of_alignments, collection_type='', collection_id=''):
        self.data = list_of_alignments
        self.collection_type = collection_type
        self.collection_id = collection_id
        self.length = len(list_of_alignments)


    def __len__(self):
        return self.length


    def _pair_strand_alignments(
        self, strand_mismatch = 3, strand_indexes = None):

        # A dictionary which directs each adapter to its
        # adapter class
        adapter_to_adapter_class = {}

        # A dictionary which directs each adapter class to 
        # the list of alignments belonging to that class
        adapter_class_to_alignments = collections.OrderedDict()

        # Create a generator of all the adapters present at this
        # position
        adapters = [alignment.adapter for alignment in self.data]
        
        # Tally these adapters using collections
        adapters_counter = collections.Counter(adapters)

        # And process them using the most common first
        for key,value in adapters_counter.most_common():
            
            # If adapter_to_adapter_class dictionary is empty, no adapter
            # families have been defined. Since most_common is returning
            # keys in order of abundance, the first key will be the most
            # frequent at this position and thus we will consider it the
            # first adapter familiy
            if len(adapter_to_adapter_class) == 0:
                adapter_to_adapter_class[key] = key
                adapter_class_to_alignments[key] = []

            # If the current key is not in the adapter_to_adapter class
            # dictionary, we need to add it. 
            elif not key in adapter_to_adapter_class:

                # First, we calculate the distance between the current key and
                # every known adapter class. 
                distance_to_adapter_class = [ 
                    nucleotide.distance(adapter_class, key, strand_indexes)
                    for adapter_class in adapter_class_to_alignments ]

                # Second, we will find the earliest index in 
                # distance_to_adapter_class which minimizes the distance
                # between the current key (adapter).
                # The earliest is selected because this represents
                # the class which existed in the highest frequency 
                min_index, min_distance = min(
                    enumerate(distance_to_adapter_class), 
                    key=operator.itemgetter(1) )

                # If the distance between the closely related adapter class 
                # and the current key is within the strand mismatch threshold,
                # we will group the current key with the closely related 
                # adapter class
                if min_distance <= strand_mismatch:
                    item = adapter_class_to_alignments.keys()[min_index]
                    adapter_to_adapter_class[key] = adapter_to_adapter_class[item]

                # Otherwise, this key is too distant from any current adapter
                # class and should be considered a new adapter class
                else:
                    adapter_to_adapter_class[key] = key
                    adapter_class_to_alignments[key] = []

        # Now that we have populated the adapter_to_adapter_class, we can
        # group the alignments into the adapter_class_to_alignments dictionary
        for alignment in self.data:
            key = adapter_to_adapter_class[alignment.adapter]
            adapter_class_to_alignments[key].append(alignment)

        return adapter_class_to_alignments


    def _pair_duplex_alignments(
        self, adapter_class_to_alignments, duplex_mismatch, duplex_indexes):

        global duplex_ids, duplex_uid, prev_id, stats

        for key in adapter_class_to_alignments:

            # If the collection type is negative, we need to swap the alpha
            # and beta adapter sequences as they are currently reversed.
            # This is necessary for proper strand comparisons but was not
            # necessary for strand collapsing.
            adapter = key
            if self.collection_type == "-":
                adapter = adapter_flip(key)

            # Fetch the identifier of these alignments
            # These will be the same for every iteration of this for loop
            # per function call
            value = adapter_class_to_alignments[key]


            id = value[0].id

            # If the previous id has changed, that is we have changed positions,
            # We should free up this space in the duplex_ids dictionary
            # because this will never be accessed again
            if not prev_id == id:

                if prev_id in duplex_ids["-"]:
                    del duplex_ids["-"][prev_id]
                
                if prev_id in duplex_ids["+"]:
                    del duplex_ids["+"][prev_id]
                
                prev_id = id

            # Create the ordered Dictionary for the id if it does not exist
            if not id in duplex_ids[self.collection_type]:
                duplex_ids[self.collection_type][id] = \
                    collections.OrderedDict()

            # If the collection type is positive, we can make the assumption
            # that we haven't seen any adapters of this type in the respected
            # minus collection because AlignmentCollectionCreater returns
            # the positive alignments at a position first, so we can just
            # assign the next uid for each adapter class
            if self.collection_type == "+":
                duplex_ids[self.collection_type][id][adapter] = duplex_uid
                duplex_uid += 1
                continue

            # If the collection type is negative, and we have not seen
            # this id in the respected positive collection, we can assume that
            # no positive alignments existed at this position because
            # AlignmentCollectionCreator would have already returned these.
            # Therefore, we can just assign the next uid for each adapter
            # class.
            # If there was a matching id in the positive collection but the
            # length of this list is empty, we can assume that a more fitting
            # duplex pairing has already occured, and thus deleted from the 
            # positive collection.
            elif not id in duplex_ids["+"] or len(duplex_ids["+"][id]) == 0:
                duplex_ids[self.collection_type][id][adapter] = duplex_uid
                duplex_uid += 1
                continue

            # Otherwise, there are items in the positive collection and
            # negative collection at this position. We need to check if any
            # of these pairs could arise from the same molecule using the
            # duplex thresholds.
            else:

                # Calculate the distance between the current negative adapter
                # and all the positive adapters for the current position.
                distance_to_adapter = [
                    nucleotide.distance(x, adapter, duplex_indexes) 
                    for x in duplex_ids["+"][id] ]

                # Fetch the earliest index which contains the minimal distance
                # as this position had the highest support for the adapter
                # class
                min_index, min_distance = min(
                    enumerate(distance_to_adapter), 
                    key=operator.itemgetter(1) )

                # If the minimum distance is within our mismatch threshold,
                # We consider this positive strand and negative strand to
                # arise from the same molecule, and thus we must associate
                # the positive uid with the negative uid.
                if min_distance <= duplex_mismatch:

                    # Fetch the adapter at min_index
                    positive_key = duplex_ids["+"][id].keys()[min_index]

                    # Group these strands by uid
                    duplex_ids["-"][id][adapter] = \
                        duplex_ids["+"][id][positive_key]

                    # If a stats container was passed, update the annealing
                    # artefacts statistics
                    # Delete the positive uid, firstly it has already
                    # been processed, and secondly, we do not want to 
                    # mistakenly pair this with another negative
                    # strand

                    del duplex_ids["+"][id][positive_key]

                # If our threshold was not met, we need to create a new
                # uid for this negative strand
                else:
                    duplex_ids[self.collection_type][id][adapter] = duplex_uid
                    duplex_uid += 1


    # def _fetch_write_data(self, adapter_class_to_alignments, consensus=False):

    #     # Here we are creating the consensus of each adapter class grouping
    #     for key in adapter_class_to_alignments:

    #         alignments = adapter_class_to_alignments[key]
    #         pos_id = str(alignments[0].id)
    #         adapter = key

    #         if self.collection_type == "-":
    #             adapter = adapter_flip(key)

    #         id = ':'.join([
    #             ''.join(['@', adapter]),
    #             pos_id,
    #             str(len(alignments)),
    #             str(self.collection_type),
    #             str(duplex_ids[self.collection_type][pos_id][adapter])
    #             ])

    #         forward_consensus = ""
    #         reverse_consensus = ""
    #         if consensus:
    #             forward_consensus = consensus([x.r1 for x in alignments], 'first')
    #             reverse_consensus = consensus([x.r2 for x in alignments], 'second')

    #         yield( (alignments, id, adapter, forward_consensus, reverse_consensus) )


    # def _update_stats(self, adapter_class_to_alignments, stats_file):

    #     global duplex_ids, duplex_stats

    #     for key in adapter_class_to_alignments:

    #         alignments = adapter_class_to_alignments[key]
    #         pos_id = str(alignments[0].id)
    #         adapter = key

    #         if self.collection_type == "-":
    #             adapter = adapter_flip(key)

    #         duplex_id = str(duplex_ids[self.collection_type][pos_id][adapter])

    #         print duplex_stats

    #         if not duplex_id in duplex_stats:
    #             duplex_stats[duplex_id] = {}
    #             duplex_stats[duplex_id]["pos_id"] = pos_id
    #             duplex_stats[duplex_id]["pos_count"] = 0
    #             duplex_stats[duplex_id]["pos_adapter"] = "NA"
    #             duplex_stats[duplex_id]["neg_count"] = 0
    #             duplex_stats[duplex_id]["neg_adapter"] ="NA"

    #         if self.collection_type == "+":
    #             duplex_stats[duplex_id]["pos_count"] = len(alignments)
    #             duplex_stats[duplex_id]["pos_adapter"] = adapter

    #         else:
    #             duplex_stats[duplex_id]["neg_count"] = len(alignments)
    #             duplex_stats[duplex_id]["neg_adapter"] = adapter


    #     if self.collection_type == "-":

    #         for duplex_id in duplex_stats:

    #             output = "\t".join([ 
    #                 duplex_stats[duplex_id]["pos_id"], \
    #                 duplex_id, \
    #                 duplex_stats[duplex_id]["pos_adapter"], \
    #                 duplex_stats[duplex_id]["neg_adapter"], \
    #                 str(duplex_stats[duplex_id]["pos_count"]), \
    #                 str(duplex_stats[duplex_id]["neg_count"]) ])

    #             stats_file.write(output + "\n")

    #             del duplex_stats[duplex_id]
                    

    # def _write_out_original(
    #     self, adapter_class_to_alignments, forward_fastq, reverse_fastq):

    #     for data in _fetch_write_data(adapter_class_to_alignments):
            
    #         alignments = data[0]
    #         id = data[1]
    #         adapter = data[2]

    #         for a in alignments:

    #             if self.collection_type == "+":
    #                 for_seq = a.r1.seq
    #                 for_qual = a.r1.qual
    #                 rev_seq = nucleotide.reverseComplement(a.r2.seq)
    #                 rev_qual = a.r1.qual[::-1]

    #             else:
    #                 for_seq = nucleotide.reverseComplement(a.r2.seq)
    #                 for_qual = a.r2.qual[::-1]
    #                 rev_seq = a.r1.seq
    #                 rev_qual = a.r1.qual

    #             original_id = ':'.join([id, a.adapter, a.qname])
                
    #             forward_next = fastq.Fastq( original_id, for_seq, '+', for_qual )
    #             reverse_next = fastq.Fastq( original_id, rev_seq, '+', rev_qual )
                
    #             forward_fastq.next(forward_next)
    #             reverse_fastq.next(reverse_next)


    # def _write_out_strand_consensus(
    #     self, adapter_class_to_alignments, forward_fastq, reverse_fastq):


    # def _write_out_duplex_consensus(
    #     self, adapter_class_to_alignments, forward_fastq, reverse_fastq):





    def _write_out_consensus(
        self, adapter_class_to_alignments, forward_fastq, reverse_fastq,
        original_forward_fastq=None, original_reverse_fastq=None):
            
        write_original = original_forward_fastq != None \
                         and original_reverse_fastq != None
        # Here we are creating the consensus of each adapter class grouping
        for key in adapter_class_to_alignments:

            value = adapter_class_to_alignments[key]
            pos_id = value[0].id
            adapter = key

            if self.collection_type == "-":
                adapter = adapter_flip(key)

            id = ':'.join([
                ''.join(['@', adapter]),
                str(pos_id),
                str(len(value)),
                str(self.collection_type),
                str(duplex_ids[self.collection_type][pos_id][adapter])
                ])

            #
            forward_consensus = consensus([x.r1 for x in value], 'first')
            reverse_consensus = consensus([x.r2 for x in value], 'second')


            # Write to Fastq
            forward_next = None
            reverse_next = None

            #
            if self.collection_type == "+":
                
                forward_next = fastq.Fastq(
                    id , forward_consensus[0], 
                    '+', forward_consensus[1] )

                reverse_next = fastq.Fastq(
                    id , nucleotide.reverseComplement(reverse_consensus[0]), 
                    '+', reverse_consensus[1][::-1] )

                forward_fastq.next(forward_next)
                reverse_fastq.next(reverse_next)
                
                if write_original:
                    
                    for v in value:
                        
                        original_id = ':'.join([id, v.adapter, v.qname])
                        
                        original_forward_next = fastq.Fastq(
                            original_id, v.r1.seq, '+', v.r1.qual )
                        
                        original_reverse_next = fastq.Fastq(
                            original_id, 
                            nucleotide.reverseComplement(v.r2.seq),
                            '+', v.r2.qual[::-1] )
                        
                        original_forward_fastq.next(original_forward_next)
                        original_reverse_fastq.next(original_reverse_next)

            #
            elif self.collection_type == "-":
                
                forward_next = fastq.Fastq(
                    id , nucleotide.reverseComplement(reverse_consensus[0]), 
                    '+', reverse_consensus[1][::-1] )

                reverse_next = fastq.Fastq(
                    id , forward_consensus[0], 
                    '+', forward_consensus[1] )

                forward_fastq.next(forward_next)
                reverse_fastq.next(reverse_next)

                if write_original:

                    for v in value:
                        
                        original_id = ':'.join([id, adapter_flip(v.adapter), v.qname])
                        
                        original_forward_next = fastq.Fastq(
                            original_id, 
                            nucleotide.reverseComplement(v.r2.seq), 
                            '+', v.r2.qual[::-1] )
                        
                        original_reverse_next = fastq.Fastq(
                            original_id, v.r1.seq,
                            '+', v.r1.qual )
                        
                        original_forward_fastq.next(original_forward_next)
                        original_reverse_fastq.next(original_reverse_next)


    def adapter_table_average_consensus(
        self, forward_fastq, reverse_fastq, adapter=None, strand_mismatch=3,
        strand_indexes=None, duplex_mismatch=3, duplex_indexes=None, 
        original_forward_fastq=None, original_reverse_fastq=None, 
        stats_file= None):


        adapter_class_to_alignments = self._pair_strand_alignments(
            strand_mismatch, strand_indexes)

        self._pair_duplex_alignments(
            adapter_class_to_alignments, 
            duplex_mismatch, duplex_indexes)

        self._write_out_consensus(
            adapter_class_to_alignments, forward_fastq, reverse_fastq,
            original_forward_fastq, original_reverse_fastq)




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

                    collection_id = self.id

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
                self.misformed[id] = align
                continue

            # Ignore Alignments that contain a Split Alignment tag
            if align.is_split():
                self.split[id] = align
                continue

            # Ignore Alignments that map to different chromosomes
            if align.is_multi_ref():
                self.multi_ref[id] = align
                continue

            # If this ID is new, push it to process
            if not id in self.tmp_collections:
                self.next_id.append(id)

            # Process the first alignment if never done before
            if self.id == None:
                self.id = self.next_id.pop(0)

            # Enter this loop where we will begin processing alignments
            while self.id.lessthan(self.base_buffer, id):

                list_of_alignments = []
                list_of_positive_alignments = []
                list_of_negative_alignments = []

                collection_id = self.id

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

                self.id = self.next_id.pop(0)

            if not id in self.tmp_collections:
                self.tmp_collections[id] = {align.qname: align}
                continue

            else:
                self.tmp_collections[id][align.qname] = align
                continue
