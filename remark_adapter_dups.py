import argparse
import printer
import pysam
import nucleotide
import collections

def get_adapter_classes(frequency):
    adapter_class = []
    for key in collections.Counter(frequency).most_common():
         if adapter_class = []:
             #This will be the most common, thus must be an adapter_class
             adapter_class.append(key);
         else:
             distances = [ nucleotide.distance(key, x) > args.max_mismatch for x in adapter_class ]
             if all(distances):
                 #If the next most common adapter differs by all current classes by the max_mismatch rate, it is an adapter class
                 adapter_class.append(key)

    return adapter_class

def get_remarked_sequences(sequences, adapters):
    keeper_sequence = {}
    for (adapt, read) in sequences:
        distances = [ nucleotide.distance(adapt, x) for x in adapters ]
        which = [ x <= args.max_mismatch for x in distances ].index(True)
        if adapters[which] in keeper_sequence:
            if keeper_sequence[adapters[which]].mapq < read.mapq:
                 keeper_sequence[adapters[which]] = read

    return keeper_sequence

def append_query_id(sequences, query_id):
    for read in sequences:
        query_id[read.qname] = "1"

def write_to_bam(sequences, bam):
    for read in sequences:
        bam.write(read)

def remark_dups(sequences, frequency, query_ids, output_bam):
    adapter_classes = get_adapter_classes(frequency)
    keeper_sequences = get_remarked_sequences(sequences, adapter_classes)
    append_query_id(keeper_sequences, query_ids)
    write_to_bam(keeper_sequences, output_bam)

if __name__ == '__main__':

    printer.general('REMARK')
    printer.general('Parsing Arguments for remark_adapter_dups.py')	
    desc = "Trim paired-end fastq files that contain an adapter sequence (paste this sequence in QNAME)"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-i", "--input",
        nargs=2,
        required=True,
        help="Takes the file locations of the forward and reverse fastq files from paired end sequencing"
        )
    parser.add_argument(
        "-o", "--output",
        nargs=2,
        required=True,
        help="Takes the file locations to store the trimmed forward and reverse fastq files"
        )
    parser.add_argument(
        "-as", "--adapter_sequence",
        type=str,
        default="WSWSWGACT",
        help="The adapter sequence following IUPAC DNA naming"
        )
    parser.add_argument(
        "-mm", "--max_mismatch",
        type=int,
        default=3,
        help="The maximum number of mismatches to accept in reads when trimming, set to the length of the adapter sequence to filter no reads"
        )
    args = parser.parse_args()

    input_bam = pysam.AlignmentFile(args.input, "r")
    output_bam = pysam.AlignmentFile(args.output, "wb", template=input_bam)

    chromosome_start = ''
    chromosome_end = ''
    start = 0
    end = 0
    sequences_forward = {}
    forward_frequency = {}
    sequences_reverse = {}
    reverse_frequency = {}
    query_ids = {}
    for read in input_bam:
        if read.qname in query_ids:
            output_bam.write(read)
            del query_ids[read.qname]

        elif read.start == start and read.end == end and read.chr == chromosome:
             tmp_adapter = read.qname.split(':')[0]
             if read.is_reverse:
                 if tmp_adapter in sequences_forward:
                     if sequences_forward[tmp_adapter].mapq < read.mapq:
                         sequences_forward[tmp_adapter] = read
                         forward_frequency[tmp_adapter] += 1
                 else:
                     sequences_forward[tmp_adapter] = read
                     forward_frequency[tmp_adapter] = 1
             else:
                 if tmp_adapter in sequences_reverse:
                     if sequences_reverse[tmp_adapter].mapq < read.mapq:
                         sequences_reverse[tmp_adapter] = read
                         reverse_frequency[tmp_adapter] += 1
                 else:
                     sequences_reverse[tmp_adapter] = read
                     reverse_frequecy[tmp_adapter] = 1

        else:
            if len(sequences_forward) == 0 and len(sequences_reverse) == 0:
                start = read.start
                end = read.end
                chromosome_start = read.chromosome
                chromosome_end = read.chromosome
                tmp_adapter = read.qname.split(':')[0]
                if read.is_reverse:
                    if tmp_adapter in sequences_forward:
                        if sequences_forward[tmp_adapter].mapq < read.mapq:
                            sequences_forward[tmp_adapter] = read
                else:
                    if tmp_adapter in sequences_reverse:
                        if sequences_reverse[tmp_adapter].mapq < read.mapq:
                            sequences_reverse[tmp_adapter] = read

            elif len(sequences_forward) == 1 or len(sequences_reverse) == 1:
                if len(sequences_forward) == 1:
								    (k, v), = sequences_forward.items()
                    query_ids[v.qname] = "1"
                    output_bam.write(v)
                if len(sequences_reverse) == 1:
                    (k, v), = sequences_reverse.items()
                    query_ids[v.qname] = "1"
                    output_bam.write(v)
                sequences_forward = {}
                sequences_reverse = {}
                forward_frequency = {}
                reverse_frequency = {}

            else:
                keeper_forward_seqs = remark_dups(sequences_forward, forward_frequency, query_ids, output_bam)
                keeper_reverse_seqs = remark_dups(sequences_reverse, reverse_frequency, query_ids, output_bam)
                sequences_forward = {}
                sequences_reverse = {}
                forward_frequency = {}
                reverse_frequency = {}
