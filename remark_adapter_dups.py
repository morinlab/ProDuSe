import argparse
import sys
import bam
import fastq
import printer

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
    nargs=1,
    type=int,
    default=3,
    help="The maximum mismatch acceptable for a new adapter to be considered in an adapter class"
    )

args = parser.parse_args()

bamfile = bam.BamOpen(args.input[0], 'r')
forward_fastq = fastq.FastqOpen(args.output[0], "w")
reverse_fastq = fastq.FastqOpen(args.output[1], "w")

counter = 1
while True:
    printer.mark
    collection = bam.CreateCollection(bamfile)
    collection['+'].consensus_average(args.max_mismatch, forward_fastq, reverse_fastq)
    collection['-'].consensus_average(args.max_mismatch, forward_fastq, reverse_fastq)

