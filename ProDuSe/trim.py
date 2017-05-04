#! /usr/bin/env python

# Usage:
#   See trim.py -h or get_args() for details
#
# Description:
#   Loops through each foward and reverse read pair in the fastq files, and compares their adapter
#   sequences with the expected adapter sequence. Reads exceeding (or below, if the option is
#   specified) the maximum allowed mismatch are discarded. Returns two new fastq files with the
#   reads that passed (or failed, if the option is specified) the filter
#
# AUTHORS:
#   Marco Albuquerque (Creator)
#   Christopher Rushton (ckrushto@sfu.ca)


# If not installed, or running in python2 this works fine
try:
    import fastq
    import nucleotide
except ImportError:
    # If installed and running in python3
    from ProDuSe import nucleotide, fastq

import configargparse
import configparser
import time
import sys
import os

"""
Processes command line arguments

    Returns:
        args: A namespace object containing parameters passed from the command line
    Raises:
        parser.error: If not two input or output files were specified
"""

desc = "Trim paired-end fastq files that contain an adapter sequence (paste this sequence in QNAME)"
parser = configargparse.ArgParser(description=desc)
parser.add(
    "-c", "--config",
    required=False,
    is_config_file=True,
    type=lambda x: is_valid_file(x, parser),
    help="An optional configuration file for any of the input arguments."
    )
parser.add(
    "-i", "--input",
    metavar="FASTQ",
    required=True,
    action="append",
    type=lambda x: is_valid_file(x, parser),
    help="A pair of fastq files containing flanked randomized adapter sequences"
    )
parser.add(
    "-o", "--output",
    required=True,
    action="append",
    type=str,
    help="Output file names"
    )
parser.add(
    "--adapter_sequence",
    type=str,
    required=True,
    help="The randomized adapter sequence flanked in input fastq files described using IUPAC bases"
    )
parser.add(
    "--adapter_position",
    type=str,
    required=True,
    help="The positions in the adapter sequence to include in distance calculations, 0 for no, 1 for yes"
    )
parser.add(
    "--max_mismatch",
    type=int,
    required=True,
    help="The maximum number of mismatches allowed between the expected and actual adapter sequences",
    )
parser.add(
    "-v",
    action="store_true",
    help="Instead, output entries that are distant from the adapter sequence"
    )
parser.add(
    "-u",
    action="store_true",
    help="Instead, output entries without trimming the adapter sequence"
    )


def is_valid_file(file, parser):
    """
    Checks to ensure the specified file exists, and throws an error if it does not

    Args:
        file: A filepath
        parser: An argparse.ArgumentParser() object. Used to throw an exception if the file does not exist

    Returns:
        type: The file itself

    Raises:
        parser.error(): An ArguParser.error() object, thrown if the file does not exist
    """

    if os.path.isfile(file):
        return file
    else:
        parser.error("Unable to find %s" % (file))


def main(args=None):

    """
        Determine the difference between the real and theoretical adapter sequences for each read in the file

        For each pair of reads in the forward and reverse fastq files, calculates the distance between
        the theoretical adapter sequence and the adapter sequence. If this distance is greater (or less
        than, if -v is specied) the maximum permitted mismatch (-mm), these reads are discarded.
        Outputs two new fastq  with the matching reads

        Args:
            args: Input arguments. See get_args() for what arguments are supported
    """
    if args is None:
        args = parser.parse_args()
    elif args.config:

        # Since configargparse does not parse commands from the config file if they are passed as argument here
        # They must be parsed manually
        cmdArgs = vars(args)
        config = configparser.ConfigParser()
        config.read(args.config)
        configOptions = config.options("config")
        for option in configOptions:
            param = config.get("config", option)
            # Convert arguments that are lists into an actual list
            if param[0] == "[" and param[-1] == "]":
                paramString = param[1:-1]
                param = paramString.split(",")

            # WARNING: Command line arguments will be SUPERSCEEDED BY CONFIG FILE ARGUMENTS
            cmdArgs[option] = param

        # This is gross, but set the flags to false
        if "v" not in cmdArgs:
            args.v = False
        if "u" not in cmdArgs:
            args.u = False

    # If input and output files were specified from the command line, ensures that a pair of files were provided
    if not len(args.input) == 2:
        parser.error('--input must be specified exactly twice (i.e. -i file1.fastq -i file2.fastq)')

    if not len(args.output) == 2:
        parser.error('--output must be specified exactly twice (i.e. -o file1.fastq -o file2.fastq)')

    print_prefix = "PRODUSE-TRIM\t"
    sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "Starting...\n"]))

    # Checks input and output files
    if not len(args.adapter_position) == len(args.adapter_sequence):
        sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "ERROR: adapter_position and adapter_sequence must have same length\n"]))
        sys.exit(1)

    if not os.path.isfile(args.input[0]):
        sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "ERROR: Input fastq files does not exist - " + args.input[0] + "\n"]))
        sys.exit(1)

    if not os.path.isfile(args.input[1]):
        sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "ERROR: Input fastq files does not exist - " + args.input[1] + "\n"]))
        sys.exit(1)

    if os.path.isfile(args.output[0]):
        sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "WARNING: Output file %s already exist, overwriting\n" % args.output[0]]))
        os.remove(args.output[0])
    if os.path.isfile(args.output[1]):
        sys.stdout.write("\t",join([print_prefix, time.strftime('%X'), "WARNING: Output file %s already exist, overwriting\n" % args.output[1]]))
        os.remove(args.output[1])

    # Determine Possible Adapters
    # TODO: Allow forward and reverse adapter sequences to be different?
    ref_adapter = ''.join([args.adapter_sequence, args.adapter_sequence])
    ref_indexes = list(''.join([args.adapter_position, args.adapter_position]))
    ref_indexes = [i for i in range(len(ref_indexes)) if ref_indexes[i] == "1"]

    # Checks if either input file is gzipped
    readOneType = 'r'
    readTwoType = 'r'
    input_one_gzipped = args.input[0].endswith(".gz")
    input_two_gzipped = args.input[1].endswith(".gz")
    if input_one_gzipped:
        readOneType = ''.join([readOneType, 'g'])
    if input_two_gzipped:
        readTwoType = ''.join([readTwoType, 'g'])

    # Checks if output files are gzipped
    writeOneType = 'w'
    writeTwoType = 'w'
    output_one_gzipped = args.output[0].endswith(".gz")
    output_two_gzipped = args.output[1].endswith(".gz")
    if output_one_gzipped:
        writeOneType = ''.join([writeOneType, 'g'])
    if output_two_gzipped:
        writeTwoType = ''.join([writeTwoType, 'g'])

    # Open Fastq files for reading and writing
    forward_input = fastq.FastqOpen(args.input[0], readOneType)
    reverse_input = fastq.FastqOpen(args.input[1], readTwoType)
    forward_output = fastq.FastqOpen(args.output[0], writeOneType)
    reverse_output = fastq.FastqOpen(args.output[1], writeTwoType)

    # For each read in the forward and reverse fastq files, trim the adapter, at it to the read id
    # and output this new fastq record to the temprary output fastq
    count, discard = 0, 0

    # Loop over every read in fastq file
    for forward_read in forward_input:

        count += 1

        if count % 100000 == 0:
            sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "Discard Rate:" + str(round(float(discard) / float(count), 3) * 100) + "%", "Count:" + str(count) + "\n"]))

        # Fetch associated reverse read
        reverse_read = reverse_input.next()

        # Get the actual adapter sequence from both fastq files
        # As well as trimming the seq and qual bases from these regions
        forward_adapter = forward_read.trim(len(args.adapter_sequence))
        reverse_adapter = reverse_read.trim(len(args.adapter_sequence))

        # Glue these together
        sequenced_adapter = ''.join([forward_adapter.seq, reverse_adapter.seq])

        # Determine the distance between the actual and expected adapter sequence
        actual_mismatch = nucleotide.distance(sequenced_adapter, ref_adapter, ref_indexes)

        # Throws out reads with less than the desired mismatch (if that option is specified)
        if args.v and actual_mismatch <= int(args.max_mismatch):

            # Do not include in output fastq files if distance is greater then max_mismatch
            discard += 1
            continue

        # Throws out reads with greater with the desired number of mismatches
        elif not args.v and actual_mismatch > int(args.max_mismatch):

            discard += 1
            continue

        if args.u:

            # Glue back the trimmed sequences and quality scores
            forward_read.prepend(forward_adapter)
            reverse_read.prepend(reverse_adapter)

        else:

            # Update read names to include adapter sequence prefixes
            forward_read.id = ''.join(['@', forward_adapter.seq, reverse_adapter.seq, ':', forward_read.id[1:]])
            reverse_read.id = ''.join(['@', forward_adapter.seq, reverse_adapter.seq, ':', reverse_read.id[1:]])

        # Write entries to fastq outputs
        forward_output.next(forward_read)
        reverse_output.next(reverse_read)

    if count % 100000 != 0:

        sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "Discard Rate:" + str(round(float(discard) / float(count), 3) * 100) + "% ", "Count:" + str(count) + "\n"]))

    # Close Fastq Files
    forward_input.close()
    reverse_input.close()
    forward_output.close()
    reverse_output.close()


if __name__ == '__main__':

    main()
