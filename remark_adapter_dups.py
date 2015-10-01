import argparse
import sys
import bam
import fastq
import printer
import pysam

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
parser.add_argument(
    "-as", "--adapter_sequence",
    default="WSWSWGACTWSWSWGACT",
    type=str
    )
#parser.add_argument(
#    "-a", "--algorithm",
#    choices = ['mapq', 'consensus'],


args = parser.parse_args()

bamfile = pysam.AlignmentFile(args.input[0], 'r')
forward_fastq = fastq.FastqOpen(args.output[0], "w")
reverse_fastq = fastq.FastqOpen(args.output[1], "w")
printer = printer.open('PROCESS_FASTQ', 'MARK', 'START', ['read'], 3)
collection_creator = bam.AlignmentCollectionCreate(bamfile, printer=printer)
for collection in collection_creator:
    collection.consensus_average(args.max_mismatch, forward_fastq, reverse_fastq, args.adapter_sequence)
