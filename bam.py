import pysam
import sys
import itertools
import fastq
import collections
import nucleotide

class AlignmentCollection:

    def __init__(self, list_of_alignments, is_adapters=True):
        self.data = list_of_alignments
        self.is_adapters = is_adapters 

    def consensus_average(self, max_mismatch, fastq_forward, fastq_reverse):
        
        if len(self.data) == 1:
            
            # Read Name, will be the adapter_class, followed by start and end points of the read
            id = ':'.join([
                ''.join(['@', self.data[0].read_one.qname.split(":")[0]]), 
                str(self.data[0].read_one.reference_id), 
                str(self.data[0].read_one.reference_start),
                str(self.data[0].read_two.reference_id),
                str(self.data[0].read_two.reference_end)
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
            for key in collections.Counter([x.read_one.qname.split(":")[0] for x in self.data]).most_common():
                
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
                adapter = alignment.read_one.qname.split(":")[0]
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
                    str(value[0].read_two.reference_end)
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
    
    def __init__(self, read_one, read_two, next_read_one, next_read_two):
        self.read_one = read_one
        self.read_two = read_two
        self.next_read_one = next_read_one
        self.next_read_two = next_read_two

    def __str__(self):
        return self.read_one.qname

    def coords(self):
        return '-'.join([
            ':'.join([
                str(self.read_one.reference_id), 
                str(self.read_one.reference_start)]),
            ':'.join([
                str(self.read_two.reference_id), 
                str(self.read_two.reference_end)])])

    def next_coords(self):
        return '-'.join([
            ':'.join([
                str(self.next_read_one.reference_id), 
                str(self.next_read_one.reference_start)]),
            ':'.join([
                str(self.next_read_two.reference_id), 
                str(self.next_read_two.reference_end)])])

    def is_positive(self):
        return self.read_one.is_read1


class BamRead:    

    def __init__(self, file):
        self.fh = pysam.AlignmentFile(file, 'r')
        self.data = {}
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
                     return Alignment(self.buffer.pop(0), self.buffer.pop(0), tmp_read, read)
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
    

def CreateCollection(bamfile, sense_matters=None):
    alignment = bamfile.next()
    if sense_matters == None:
        sense_matters = True
    if sense_matters:
        collections = {}
        collections['+'] = []
        collections['-'] = []    
        if alignment.is_positive():
            collections['+'].append(alignment)
        else:
            collections['-'].append(alignment)
        while alignment.coords() == alignment.next_coords():
            alignment = bamfile.next()
            if alignment.is_positive():
                collections['+'].append(alignment)
            else:
                collections['-'].append(alignment)
        collections['+'] = AlignmentCollection(collections['+'])
        collections['-'] = AlignmentCollection(collections['-'])
        return collections
    else:
        collections = [alignment]
        while alignment.coords() == alignment.next_coords():
            alignment = bamfile.next()
            collections.append(alignment)
        collections = AlignmentCollection(collections)
        return collections
