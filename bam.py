import pysam
import sys
import itertools
import fastq
import collections
import nucleotide
import json

class AlignmentCollectionCreate:

    def __init__(self, bamfile, sense_matters=True):
        self.data = bamfile
        self.tmp_collections = {}
        self.collection_count_read_one = {}
        self.collection_count_read_two = {}
        self.ref = None
        self.start = None
        self.multi_ref_reads = {}
        self.qnames_to_id = {}
        self.read_queue = []
        self.next_collection_to_output = []
        self.sense_matters = sense_matters
        self.tmp_collection_sense = []

    def __str__(self):
       return json.dumps(self.qnames_to_id)

    def next(self):

        if not self.tmp_collection_sense == []:
            return self.tmp_collection_sense.pop(0)

        elif self.next_collection_to_output == []:
            while True:
                # Obtain a Read
                read = self.data.next()
                # Remove Reads that map to different chromosomes for Now
                while not read.reference_id == read.next_reference_id:
                    self.multi_ref_reads[read.qname] = read
                    read = self.data.next()

                # Process the first read if never done before
                if self.start == None :
                    self.ref = read.reference_id
                    self.start = read.reference_start
                    id = ':'.join([
                        str(read.reference_id), 
                        str(read.reference_start), 
                        str(read.next_reference_id), 
                        str(read.reference_start + abs(read.tlen))
                        ])
                    self.collection_count_read_one[id] = 1
                    self.collection_count_read_two[id] = 0
                    self.qnames_to_id[read.qname] = id
                    self.tmp_collections[id] = { read.qname : [ read ] }
                    read = self.data.next()

			    # Process Reads until the start position changes
                while True:

                    while not read.reference_id == read.next_reference_id:
                        self.multi_ref_reads[read.qname] = read
                        read = self.data.next()

                    if read.qname in self.qnames_to_id:
                        self.read_queue.append(read)

                    else:               
                        id = ':'.join([
                            str(read.reference_id), 
                            str(read.reference_start), 
                            str(read.next_reference_id), 
                            str(read.reference_start + abs(read.tlen))
                            ])

                        if not id in self.tmp_collections:
                            self.tmp_collections[id] = { read.qname : [ read ] }
                            self.collection_count_read_one[id] = 1
                            self.collection_count_read_two[id] = 0
                            self.qnames_to_id[read.qname] = id

                        else:
                            self.tmp_collections[id][read.qname] = [read]
                            self.collection_count_read_one[id] += 1
                            self.qnames_to_id[read.qname] = id

                    if self.start < read.reference_start or self.ref < read.reference_id:
                        self.start = read.reference_start
                        self.ref = read.reference_id
                        break

                    read = self.data.next()

                # Process Reads in Queue
                if not self.read_queue == []:
                    while True:
                        read = self.read_queue.pop(0)
                        id = self.qnames_to_id[read.qname]
                        self.tmp_collections[id][read.qname].append(read)
                        self.collection_count_read_two[id] += 1
                        if(self.collection_count_read_two[id] == self.collection_count_read_one[id]):
                            self.next_collection_to_output.append(id)

                        if self.read_queue == []:
                            break

                # Process Reads in Collection Queue
                if not self.next_collection_to_output == []:
                    break

        
        collection_id = self.next_collection_to_output.pop(0)
        list_of_alignments_positive = []
        list_of_alignments_negative = []
        for pairs in self.tmp_collections[collection_id]:
            tmp_align = Alignment(self.tmp_collections[collection_id][pairs][0], self.tmp_collections[collection_id][pairs][0])
            
            del self.qnames_to_id[tmp_align.read_one.qname]
            if tmp_align.is_positive():
                list_of_alignments_positive.append(tmp_align)
            else:
                list_of_alignments_negative.append(tmp_align)
        del self.tmp_collections[collection_id]
        del self.collection_count_read_one[collection_id]
        del self.collection_count_read_two[collection_id]
        if not list_of_alignments_positive == []:
            self.tmp_collection_sense.append(AlignmentCollection(list_of_alignments_positive, collection_type='+', collection_id=collection_id))
        if not list_of_alignments_negative == []:
            self.tmp_collection_sense.append(AlignmentCollection(list_of_alignments_negative, collection_type='-', collection_id=collection_id))
        return self.tmp_collection_sense.pop(0)

class AlignmentCollection:

    def __init__(self, list_of_alignments, is_adapters=True, collection_type='', collection_id=''):
        self.data = list_of_alignments
        self.is_adapters = is_adapters 
        self.length = len(list_of_alignments)
        self.collection_type = collection_type
        self.collection_id = collection_id

    def __len__(self):
        return self.length

    def __str__(self):
        tmp_str = '\t'.join([self.collection_type, self.collection_id])
        return tmp_str


    def consensus_average(self, max_mismatch, fastq_forward, fastq_reverse):
        
        if len(self.data) == 1:
            
            # Read Name, will be the adapter_class, followed by start and end points of the read
            id = ':'.join([
                ''.join(['@', self.data[0].read_one.qname.split("M00")[0]]), 
                str(self.data[0].read_one.reference_id), 
                str(self.data[0].read_one.reference_start),
                str(self.data[0].read_two.reference_id),
                str(self.data[0].read_two.reference_end),
                str(len(self.data))
                ])

            # Since there is only a single Pair in this adapter class group, just output the same seq and qual
            # May need to do some reverse complement magic here. We will see post alignment
            forward_seq = self.data[0].read_one.seq
            forward_qual = self.data[0].read_one.qual
            reverse_seq = self.data[0].read_two.seq
            reverse_qual = self.data[0].read_two.qual
            
            # Generate the Fastq entries and write them to file
            forward = fastq.Fastq(id, forward_seq, "+", forward_qual)
            reverse = fastq.Fastq(id, reverse_seq, "+", reverse_qual)
            fastq_forward.next(forward)
            fastq_reverse.next(reverse)

        elif self.is_adapters:
            
            # We are trying to identify all of the adapter classes in the collection of Reads
            adapter_to_adapter_class = {}
            adapter_class = []
            for key in collections.Counter([x.read_one.qname.split("M00")[0] for x in self.data]):
                
                # If the adapter_class list is empty, add the most frequent adapter as the first adapter class
                if len(adapter_class) == 0:
                    adapter_class.append(key)
                    adapter_to_adapter_class[key] = key

                # Otherwise we need to determine if this new adapter is its own class or already in a class
                elif not key in adapter_to_adapter_class:
                    
                    # Calculate the distance of the adapter to each known adapter classes
                    distance_to_adapter_class = [ nucleotide.distance(x, key) for x in adapter_class ]
                    # Determine which class (i.e. find the index) the adapter belongs to
                    index = [ n for n,i in enumerate(distance_to_adapter_class) if i <= max_mismatch ]
                    
                    # If no index was created, the adapter is distant from all current adapter classes
                    # Define this as a new adapter class
                    if index == []:
                        adapter_class.append(key)
                        adapter_to_adapter_class[key] = key
 
                    # Otherwise, update adapter to adapter class dict
                    else:
                        index = index[0]
                        adapter_to_adapter_class[key] = adapter_class[index]


						# Here we are grouping all adlignments to their respected adapter class
            adapter_class_to_reads = {}
            for alignment in self.data:
                adapter = alignment.read_one.qname.split("M00")[0]
                if adapter_to_adapter_class[adapter] in adapter_class_to_reads:
                    adapter_class_to_reads[adapter_to_adapter_class[adapter]].append(alignment)
                else:
                    adapter_class_to_reads[adapter_to_adapter_class[adapter]] = [alignment]

            # Here we are creating the consensus of each adapter class grouping
            for key in adapter_class_to_reads.iterkeys():
                value = adapter_class_to_reads[key]
                id = ':'.join([
                    ''.join(['@', key]), 
                    str(value[0].read_one.reference_id),
                    str(value[0].read_one.reference_start),
                    str(value[0].read_two.reference_id),
                    str(value[0].read_two.reference_end),
                    str(len(self.data))
                    ])

                # Construct consensus of Seq and Qual for forward reads
                seq_forward = []
                qual_forward = []
                for i in range(0, min([len(x.read_one.seq) for x in value])-1):
                    seq_forward.append(collections.Counter([ x.read_one.seq[i-1] for x in value]).most_common(1)[0][0])
                    qual_forward.append(collections.Counter([ x.read_one.qual[i-1] for x in value]).most_common(1)[0][0])
                
                # Construct consensus of Seq and Qual for reverse reads
                seq_reverse = []
                qual_reverse = []    
                for i in range(0, min([len(x.read_two.seq) for x in value])-1):
                    seq_reverse.append(collections.Counter([ x.read_two.seq[i-1] for x in value]).most_common(1)[0][0])
                    qual_reverse.append(collections.Counter([ x.read_two.qual[i-1] for x in value]).most_common(1)[0][0])

                # Write to Fastq
                fastq_forward.next(fastq.Fastq(str(id), ''.join(seq_forward), '+', ''.join(qual_forward)))
                fastq_reverse.next(fastq.Fastq(str(id), ''.join(seq_reverse), '+', ''.join(qual_reverse)))                

class Alignment:
    
    def __init__(self, read_one, read_two):
        self.read_one = read_one
        self.read_two = read_two

    def __str__(self):
        return "{0}:{1}:{2}:{3}:{4}:{5}".format(
            self.read_one.qname,
            str(self.read_one.reference_id),
            str(self.read_one.reference_start),
            self.read_two.qname,
            str(self.read_two.reference_id),
            str(self.read_two.reference_end))

    def coords(self):
        return '-'.join([
            ':'.join([
                str(self.read_one.reference_id), 
                str(self.read_one.reference_start)]),
            ':'.join([
                str(self.read_two.reference_id), 
                str(self.read_two.reference_end)])])

    def is_positive(self):
        return self.read_one.is_read1


class BamRead:    

    def __init__(self, file):
        self.fh = pysam.AlignmentFile(file, 'r')
        self.data = {}
        self.start = 0
        self.buffer = []

    def __iter__(self):
        return self

    def next(self):
        
        while True:
            read = self.fh.next()
            if ''.join(['1:', read.qname]) in self.data:
                tmp_read = self.data[''.join(['1:', read.qname])]
                del self.data[''.join(['1:', read.qname])]
                if self.buffer == []:
                     self.buffer.append(tmp_read)
                     self.buffer.append(read)
                else:
                     self.buffer.append(tmp_read)
                     self.buffer.append(read)
                     return Alignment(self.buffer.pop(0), self.buffer.pop(0))
            else:
                self.data[''.join(['1:', read.qname])] = read

    def __del__(self):
        self.close()

    def __exit__(self):
        self.close()
    
    def close(self):
        self.fh.close()

def BamOpen(file, mode):
    if mode == 'r':
        return BamRead(file)
    else:
        sys.exit(':'.join(['This mode is currently not supported by BamOpen', mode]))
    
