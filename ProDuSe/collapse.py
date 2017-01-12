#!/usr/bin/env python

# Imported Python Modules
import configargparse
import fastq
import nucleotide
import re
import sys
import os
import pysam
import alignment
import printer


# Command Line Arguments
desc = "Create strand specific consensus for independent molecules present in BAM files."

parser = configargparse.ArgParser( description=desc )

parser.add(
    "-c", "--config",
    required=False,
    is_config_file=True,
    help="Configuration file for any of the input arguments."
    )

parser.add(
    "-i", "--input",
    type=str,
    required=True,
    help="An input bam file for reading generated from trim.py and bwa.py"
    )

parser.add(
    "-o", "--output",
    required=True,
    action="append",
    help="A pair of empty fastq files for writing"
    )

parser.add(
    "-as", "--adapter_sequence",
    type=str,
    required=True,
    help="The randomized adapter sequence flanked in input fastq files described using IUPAC bases"
    )

parser.add(
    "-sp", "--strand_position",
    type=str,
    required=True,
    help="The positions in the adapter sequence to include in distance calculations, 0 for no, 1 for yes"
    )

parser.add(
    "-dp", "--duplex_position",
    type=str,
    required=True,
    help="The positions in the adapter sequence to include in distance calculations, 0 for no, 1 for yes"
    )

parser.add(
    "-smm", "--strand_max_mismatch",
    type=int,
    required=True,
    help="The maximum number of mismatches allowed between the expected and actual adapter sequences",
    )

parser.add(
    "-dmm", "--duplex_max_mismatch",
    type=int,
    required=True,
    help="The maximum number of mismatches allowed between the expected and actual adapter sequences",
    )

parser.add(
    "-mamt", "--max_alignment_mismatch_threshold",
    type=int,
    required=False,
    default=20,
    help="The maximum number of mismatches allowed in an alignment"
    )

parser.add(
    "-oo", "--original_output",
    type=str,
    required=False,
    action="append",
    default = None,
    help="A pair of empty fastq files to rewrite original fastqs with duplex information"
    )

parser.add(
    "-sf", "--stats_file",
    type=str,
    required=False,
    default = None,
    help="A single file to store stats on the run"
    )


def main(args=None):

    # Fetch command line arguments
    if args == None:
        args = parser.parse_args()

    # Check common errors
    if not len(args.strand_position) == len(args.adapter_sequence):
        printer.message("strand_position and adapter_sequence must have same length", error=True)

    if not len(args.duplex_position) == len(args.adapter_sequence):
        printer.message("duplex_position and adapter_sequence must have same length", error=True)

    if not os.path.isfile(args.input):
        printer.message("input bam file does not exist - " + args.input, error=True)

    if os.path.isfile(args.output[0]):
        printer.message("output fastq file(s) already exist - " + args.output[0], error=True)

    if os.path.isfile(args.output[1]):
        printer.message("output fastq file(s) already exist - " + args.output[1], error=True)

    print args.original_output

    if args.original_output != None:
        if os.path.isfile(args.original_output[0]):
            printer.message("original output fastq file(s) already exist - " + args.original_output[0], error=True)

        if os.path.isfile(args.original_output[1]):
            printer.message("original output fastq file(s) already exist - " + args.original_output[1], error=True)

    write = 'w'
    is_output_one_gzipped = not not re.search('.*\.gz', args.output[0])
    is_output_two_gzipped = not not re.search('.*\.gz', args.output[1])
    if is_output_one_gzipped and is_output_two_gzipped:
        write = ''.join([write, 'g']);
    elif is_output_one_gzipped or is_output_two_gzipped:
        printer.message("output files must be both gzipped or both uncompressed", error=True)

    original_write = 'w'
    if args.original_output != None:
        is_output_one_gzipped = not not re.search('.*\.gz', args.original_output[0])
        is_output_two_gzipped = not not re.search('.*\.gz', args.original_output[1])
        if is_output_one_gzipped and is_output_two_gzipped:
            original_write = ''.join([original_write, 'g']);
        elif is_output_one_gzipped or is_output_two_gzipped:
            printer.message("original output files must be both gzipped or both uncompressed", error=True)

    # Everything looks good to go
    printer.message("Starting...")

    bamfile = pysam.AlignmentFile(args.input, 'r')

    ref_adapter = ''.join([args.adapter_sequence, args.adapter_sequence])
    strand_indexes = list(''.join([args.strand_position, args.strand_position]))
    strand_indexes = [ i for i in range(len(strand_indexes)) if strand_indexes[i] == "1" ]
    duplex_indexes = list(''.join([args.duplex_position, args.duplex_position]))
    duplex_indexes = [ i for i in range(len(duplex_indexes)) if duplex_indexes[i] == "1" ]
    forward_fastq = fastq.FastqOpen(args.output[0], write)
    reverse_fastq = fastq.FastqOpen(args.output[1], write)


    original_forward_fastq = None
    original_reverse_fastq = None
    if args.original_output != None:
        original_forward_fastq = fastq.FastqOpen(args.original_output[0], original_write)
        original_reverse_fastq = fastq.FastqOpen(args.original_output[1], original_write)

    stats_file = None
    if not args.stats_file == None:
        stats_file = open(args.stats_file, 'w')

    collection_creator = alignment.AlignmentCollectionCreate(bamfile, max_alignment_mismatch_threshold=args.max_alignment_mismatch_threshold)
    counter = 0
    for collection in collection_creator:
        counter += 1
        collection.adapter_table_average_consensus(
            forward_fastq = forward_fastq, 
            reverse_fastq = reverse_fastq,
            strand_mismatch = args.strand_max_mismatch, 
            strand_indexes = strand_indexes, 
            duplex_mismatch = args.duplex_max_mismatch, 
            duplex_indexes = duplex_indexes,
            original_forward_fastq = original_forward_fastq, 
            original_reverse_fastq = original_reverse_fastq,
            stats_file = stats_file )

        if counter % 100000 == 0:
            printer.message("Positions Processed:" + str(counter))
    printer.message("Positions Processed:" + str(counter))

if __name__ == "__main__":
    main()
