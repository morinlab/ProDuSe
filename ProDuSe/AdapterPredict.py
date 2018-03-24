#! /usr/bin/env python

import argparse
import os
import sys
import gzip
import math


# Which nucleotides each IUPAC base represents, and their distribution
# A:0, C:1, G:2, T:3
DIST = {
    'A':(1,0,0,0),
    'C':(0,1,0,0),
    'G':(0,0,1,0),
    'T':(0,0,0,1),
    'W':(0.5,0,0,0.5),
    'S':(0,0.5,0.5,0),
    'M':(0.5,0.5,0,0),
    'K':(0,0,0.5,0.5),
    'R':(0.5,0,0.5,0),
    'Y':(0,0.5,0,0.5),
    'N':(0.25,0.25,0.25,0.25)
    }


def isValidFile(file, parser):
    """
    Checks to ensure the specified file exists

    :param file: A string containing a filepath
    :param parser: An argparse.ArgumentParser object
    :return: The input variable file, if it exists
    :raises: parser.error() if the specified file does not exist
    """

    if os.path.exists(file):
        return file
    else:
        raise parser.error("Unable to locate \'%s\'. Please ensure the file exists, and try again" % file)


def getLikelyBase(nucFreq):

    # So here's the deal. We need to figure out the consensus base at this position
    # Since it is degenerate, it could be one of several bases
    # However, if we assume that, in the case of a degenerate position, the nucleotides should be present at an equal proportion
    # then we can identify a most likely degenerate base based upon the proportions of each base at this position

    # First, determine the proportion of each nucleotide
    try:
        total = float(sum(nucFreq[i] for i in range(len(nucFreq))))
        prop = tuple(float(x) / total for x in nucFreq)
    except ZeroDivisionError:  # Aka, there are no nucleotides whatsoever at this position. The user likely entered an insane barcode length
        return None

    minDist = 10
    minBase = None

    # Check each combination of items
    for base, distances in DIST.items():
        # Determine the distance between the expected bases at this position, and the actual bases
        realDist = sum((float(distances[i]) - float(prop[i])) ** 2 for i in range(0, 4))
        realDist = math.sqrt(realDist)

        if realDist < minDist:
            minDist = realDist
            minBase = base

    return minBase


parser = argparse.ArgumentParser(description="Estimates the degenerate barcode sequence used in a given sample")
parser.add_argument("-i", "--input", metavar="FASTQ", nargs=2, required=True, type=lambda x: isValidFile(x, parser),
                        help="A set of paired-end FASTQ files")
parser.add_argument("-m", "--max_barcode_length", metavar="INT", default=16, type=int, help="Maximum barcode length to predict [Default: %(default)s]")


def main(args=None, sysStdin=None, supressOutput=False):

    if args is None:
        if sysStdin is None:
            args = parser.parse_args()
        else:
            args = parser.parse_args(sysStdin)

    # Sanity check: Is the maximum barcode length >1?
    if args.max_barcode_length <= 0:
        raise parser.error("-m/--max_adapter_length must be greater than 0")

    # Sanity check: Are the input FASTQ files the same?
    if os.path.samefile(args.input[0], args.input[1]):
        raise parser.error("\'%s\' and \'%s\' are the same file" % (args.input[0], args.input[1]))

    # Sanity check: Is the maximum barcode length something stupid?
    if args.max_barcode_length > 100:
        sys.stderr.write("WARNING: m/--max_adapter_length was set to %s, which is kind of insane.\n" % args.max_barcode_length)
        sys.stderr.write("We'll continue anyways, but you should really check that this is correct.\n")

    # Generate a tuple which will store the bases at each position in the adapter sequence
    baseCounts = []
    for i in range(0, args.max_barcode_length):
        baseCounts.append([0, 0, 0, 0, 0])
    baseCounts = tuple(baseCounts)

    # Open the input files for reading
    # Determine if they are gzipped based upon file extension
    # I have tried examining the magic number before, but that does not work
    if args.input[0].split(".")[-1] == "gz":  # Aka get the file extension
        f1in = gzip.open(args.input[0], "rt")
    else:
        f1in = open(args.input[0], "r")

    if args.input[1].split(".")[-1] == "gz":  # Aka get the file extension
        f2in = gzip.open(args.input[1], "rt")
    else:
        f2in = open(args.input[1], "r")

    lineNum = 0
    nucToIndex = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}

    for forwardRead, reverseRead in zip(f1in, f2in):
        lineNum += 1
        # Since the nucleotide sequence of each FASTQ record is only on line 2, ignore the other 3 lines
        if (lineNum - 2) % 4 != 0:
            continue

        # Cycle through the bases of this sequence which correspond to the barcode
        for i in range(0, args.max_barcode_length):
            try:
                baseCounts[i][nucToIndex[forwardRead[i]]] += 1  # aka obtain the base at the i'th position of the sequence, and add that
                baseCounts[i][nucToIndex[reverseRead[i]]] += 1  # count to the base count dictionary for both reads

            except KeyError:  # It looks like there are non-nucleotide sequences in one of the input files, or we have reached the end of the read
                if forwardRead[i] == os.linesep or reverseRead[i] == os.linesep:  # Aka, we have reached the end of this read pair. It is truncated
                    break
                else:
                    sys.stderr.write("ERROR: \'%s\' or \'%s\' do not appear to be nucleotides" % (forwardRead[i], reverseRead[i]))
                    exit(1)

    # Identify the base used at each position in the barcode
    barcode = ""
    for pos in baseCounts:

        newBase = getLikelyBase(pos[:4])
        if newBase:
            barcode += newBase
        else:
            break

    if not supressOutput:
        sys.stdout.write(barcode + "\n")
    return barcode

if __name__ == "__main__":
    main()
