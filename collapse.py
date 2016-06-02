import argparse
import sys
import alignment
import fastq
import printer
import pysam
import re

desc = "Trim paired-end fastq files that contain an adapter sequence (paste this sequence in QNAME)"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument(
    "-i", "--input",
    nargs=1,
    required=True,
    help="Takes the file locations of the forward and reverse fastq files from paired end sequencing"
    )
parser.add_argument(
    "-o", "--output",
    nargs=2,
    required=True,
    help="Takes the file location of the forward and reverse fastq files"
    )
parser.add_argument(
    "-mm", "--max_mismatch",
    type=int,
    default=3,
    help="The maximum mismatch acceptable for a new adapter to be considered in an adapter class"
    )
parser.add_argument(
    "-as", "--adapter_sequence",
    default="WSWSWGACT",
    type=str
    )

args = parser.parse_args()

bamfile = pysam.AlignmentFile(args.input[0], 'r')

# Check If Output is Gzip and call appropariate FastqOpen
write = 'w'
is_output_one_gzipped = not not re.search('.*\.gz', args.output[0])
is_output_two_gzipped = not not re.search('.*\.gz', args.output[1])
if is_output_one_gzipped and is_output_two_gzipped:
    write = ''.join([write, 'g']);
elif is_output_one_gzipped or is_output_two_gzipped:
    print 'Output files must both be gzipped or both uncompressed\n'
    sys.exit()

forward_fastq = fastq.FastqOpen(args.output[0], write)
reverse_fastq = fastq.FastqOpen(args.output[1], write)

collection_creator = alignment.AlignmentCollectionCreate(bamfile)
for collection in collection_creator:
    collection.adapter_table_average_consensus(forward_fastq, reverse_fastq, args.max_mismatch, args.adapter_sequence)