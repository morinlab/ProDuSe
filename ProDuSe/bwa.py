#! /usr/bin/env python

# USAGE:
#   See bwa.py -h for details
#
# DESCRIPTION:
#   Maps FASTA sequences to a reference genome using the Burrows-Wheeler algoririthm
#   And sorts the resulting BAM file
#   Designed for use in the ProDuSe pipeline
#
# AUTHORS
#   Marco Albuquerque (Creator)
#   Christopher Rushton (ckrushto@sfu.ca)

import configargparse
import configparser
import sys
import os
import subprocess
import time
import itertools
from distutils.version import StrictVersion


"""
Processes command line arguments

Returns:
    args: A namespace object containing parameters passed from the command line
Raises:
    parser.error: If not two input files were specified
"""
desc = "Aligns paired fastq files using BWA. For use with the PRODUSE Pipeline"
parser = configargparse.ArgParser(description=desc)
parser.add(
    "-c", "--config",
    required=False,
    is_config_file=True,
    help="An optional configuration file, can specify any of the input arguments."
    )
parser.add(
    "-t", "--threads",
    required=False,
    default=1,
    help="Number of threads BWA should use [Default: %(default)s] "
    )
parser.add(
    "-i", "--input",
    required=True,
    action="append",
    type=lambda x: is_valid_file(x, parser),
    help="Fastq files, coresponding to the forward and reverse reads"
    )
parser.add(
    "-o", "--output",
    required=True,
    type=str,
    help="Output BAM file name"
    )
parser.add(
    "-r", "--reference",
    required=True,
    type=lambda x: is_valid_file(x, parser),
    help="Reference genome file with BWA indexes"
    )
parser.add(
    "-R", "--readgroup",
    required=False,
    type=str,
    help="Fastq read group"
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
    Runs BWA on the supplied files, and sorts the output file into


    Args:
        args: A namespace object listing BWA paramters. See get_args() for supported paramters
    """

    if args is None:
        args = parser.parse_args()

        # Ensures bwa is installed on the system
        try:
            DEVNULL = open(os.devnull, "w")
            checkBWA = subprocess.Popen("bwa", stderr=DEVNULL, stdout=DEVNULL)
            checkBWA.wait()
        except OSError:
            sys.stderr.write("ERROR: Burrows-Wheeler Aligner (bwa) cannot be found\n")
            sys.stderr.write("Please install BWA and try again\n")
            sys.exit(1)

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

    # If input and output files were specified from the command line, ensures that a pair of files were provided
    if not len(args.input) == 2:
        parser.error('--input must be specified exactly twice (i.e. -i file1.fastq -i file2.fastq)')


    print_prefix = "PRODUSE-BWA\t"
    sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "Starting...\n"]))

    if not os.path.isfile(args.input[0]) or not os.path.isfile(args.input[1]):
        sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "ERROR: one or more of input fastq files does not exist - " + args.input[0] + args.input[1] + "\n"]))
        sys.exit(1)

    if os.path.isfile(args.output):
        sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "ERROR: output file already exists - " + args.output + "\n"]))
        sys.exit(1)

    # Determine samtools version and ensures it is V 1.3.1 or higher
    version_samtools_command1 = ['samtools']

    ps1 = subprocess.Popen(version_samtools_command1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps1_iterator = iter(ps1.stderr.readline, b"")
    for line in ps1_iterator:
        line = line.decode("utf-8")
        if line.startswith("Version"):
            if StrictVersion(line.split(" ")[1]) < StrictVersion("1.3.1"):
                sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "Error: samtools must be Version 1.3.1 or higher - " + line.split(" ")[1] + "\n"]))
                sys.exit(1)
            break


    # Sets up BWA args
    threads = "-t " + str(args.threads)
    bwaCom = ['bwa', 'mem', threads]
    if args.readgroup:
        bwaCom.append(args.readgroup)
    bwaCom.extend([args.reference, args.input[0], args.input[1]])

    # Samtools arguments
    samtoolsCom = ['samtools', 'sort', '--output-fmt', "BAM", '-']
    outFile = args.output

    # Runs processes
    with open(outFile, "w") as o:

        runBwa = subprocess.Popen(bwaCom, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        runSamtools = subprocess.Popen(samtoolsCom, stdout=o, stdin=runBwa.stdout, stderr=subprocess.PIPE)

        bwa_iterator = iter(runBwa.stderr.readline, b"")
        samtools_iterator = iter(runSamtools.stderr.readline, b"")

        counter = 0
        # Displays a status update at the command line
        for line in itertools.chain(bwa_iterator, samtools_iterator):
            line = line.decode("utf-8")
            if line.startswith("[M::mem_process_seqs]"):
                counter += int(line.split(" ")[2])
                sys.stdout.write("\t".join([print_prefix, time.strftime('%X'), "Reads Processed:" + str(counter) + "\n"]))
            elif line.startswith("[E::"):
                sys.stderr.write(line + "\n")

    runBwa.stdout.close()
    runBwa.wait()
    runSamtools.wait()
    bwaReturnCode = runBwa.returncode

    # Ensure that BWA completed sucessfully
    if bwaReturnCode != 0:
        sys.exit(1)


if __name__ == '__main__':

    main()
