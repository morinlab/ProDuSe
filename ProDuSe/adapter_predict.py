import argparse
import printer
import fastq
import nucleotide
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

desc = "Quality Control on the adapter sequences in a pair of fastq files"
parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--input",
    nargs=2,
    required=True,
    help="Takes the file locations of the forward and reverse fastq files from paired end sequencing"
    )
parser.add_argument(
    "-m", "--max_adapter_length",
    default=30,
    type=int
    )

def main(args=None):

    if args == None:
        args = parser.parse_args()
    
    # Check If Input is Gzip and call appropariate FastqOpen
    read = 'r'
    is_input_one_gzipped = not not re.search('.*\.gz', args.input[0])
    is_input_two_gzipped = not not re.search('.*\.gz', args.input[1])
    if is_input_one_gzipped and is_input_two_gzipped:
        read = ''.join([read, 'g']);
    elif is_input_one_gzipped or is_input_two_gzipped:
        print 'Input files must both be gzipped or both uncompressed\n'
        sys.exit()

    # Open Fastq files for reading
    forward_input = fastq.FastqOpen(args.input[0], read)
    reverse_input = fastq.FastqOpen(args.input[1], read)

    counts = [[0 for i in range(5)] for j in range(args.max_adapter_length)]

    for forward_read in forward_input:
        reverse_read = reverse_input.next()
        for i in range(args.max_adapter_length):
            counts[i][nucleotide.BASE_TO_INDEX[forward_read.seq[i]]] += 1
            counts[i][nucleotide.BASE_TO_INDEX[reverse_read.seq[i]]] += 1

    predicted_seq = ""

    for i in range(args.max_adapter_length):
        total = float(sum(counts[i][0:4]))
        props = [float(val) / total for val in counts[i][0:4] ]
        min_dist = 1
        min_base = ''
        for real_base, real_props in nucleotide.DIST.iteritems():
            new_dist = nucleotide.dist(real_props, props)
            if new_dist < min_dist:
                min_dist = new_dist
                min_base = real_base
        predicted_seq += min_base 

    print predicted_seq + "\n"


if __name__ == "__main__":
    main()
