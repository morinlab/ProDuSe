#! /usr/bin/env python

import argparse
import os
import gzip
import time
import sys
from configobj import ConfigObj


def isValidFile(file, parser):
    """
    Checks to ensure the provided file is exists, and throws an error if it is not.

    :param file: A string containing a filepath to the file of interest
    :param parser: An argparse.ArgumentParser() object.

    :returns: The "file" variable, if the file is valid
    :raises parser.error: If the file is not valid
    """

    if os.path.exists(file):
        return file
    else:
        raise parser.error("Unable to locate %s. Please ensure the file exists, and try again." % (file))


def validateArgs(args):
    """
    Validates that the specified set of arguments are valid
    :param args: A dictionary listing {argument: parameter}
    :return: A dictionary listing {argument: parameter} that have been validatd
    """

    # Convert the dictionary into a list, to allow the arguments to be re-parsed
    listArgs = []
    for argument, parameter in args.items():

        if parameter is None or parameter is False:
            continue
        # Something was provided as an argument
        listArgs.append("--" + argument)

        # If the argument is a flag, ignore the boolean, as it will be re-added once the arguments re-parsed
        if isinstance(parameter, bool):
            continue

        # If the parameter is a list, we need to add each element seperately
        if isinstance(parameter, list):
            for p in parameter:
                listArgs.append(str(p))
        else:
            listArgs.append(str(parameter))

    parser = argparse.ArgumentParser(description="Trims degenerate barcodes")
    parser.add_argument("-c", "--config", metavar="INI", type=lambda x: isValidFile(x, parser),
                        help="An optional configuration file, which can provide one or more arguments")
    parser.add_argument("-i", "--input", metavar="FASTQ", required=True, nargs=2, type=lambda x: isValidFile(x, parser),
                        help="A pair of FASTQ files generated from paired-end sequencing. May be gzipped")
    parser.add_argument("-o", "--output", metavar="FASTQ", required=True, nargs=2,
                        help="A pair of output FASTQ files, to which the trimmed reads will be written to. May be gzipped")
    parser.add_argument("-mm", "--max_mismatch", metavar="INT", required=True, type=int,
                        help="Maximum mismatch allowed between expected and actual adapter sequences")
    parser.add_argument("-b", "--barcode_sequence", metavar="NNNWSMRWSYWKMWWT", required=True, type=str,
                        help="Degenerate barcode sequence, represented in IUPAC bases")
    parser.add_argument("-p", "--barcode_position", metavar="0001111111111110", required=True, type=str,
                        help="Positions to consider when comparing expected and actual barcode sequences (1=Yes, 0=No)")
    parser.add_argument("--reverse", action="store_true",
                        help="Instead, output reads which fall outside the mismatch threshold")
    parser.add_argument("--no_trim", action="store_true", help="Instead, output entries without trimming the adapter sequence")
    parser.add_argument("--trim_other_end", action="store_true",
                        help="In addition, examine the other end of the read for a barcode. Will not remove partial barcodes")
    validatedargs = parser.parse_args(listArgs)
    return vars(validatedargs)


parser = argparse.ArgumentParser(description="Trims degenerate barcodes")
parser.add_argument("-c", "--config", metavar="INI", type=lambda x: isValidFile(x, parser),
                    help="An optional configuration file, which can provide one or more arguments")
parser.add_argument("-i", "--input", metavar="FASTQ", nargs=2, type=lambda x: isValidFile(x, parser),
                    help="A pair of FASTQ files generated from paired-end sequencing. May be gzipped")
parser.add_argument("-o", "--output", metavar="FASTQ", nargs=2,
                    help="A pair of output FASTQ files, to which the trimmed reads will be written to. May be gzipped")
parser.add_argument("-mm", "--max_mismatch", metavar="INT", type=int,
                    help="Maximum mismatch allowed between expected and actual adapter sequences")
parser.add_argument("-b", "--barcode_sequence", metavar="NNNWSMRWSYWKMWWT", type=str,
                    help="Degenerate barcode sequence, represented in IUPAC bases")
parser.add_argument("-p", "--barcode_position", metavar="0001111111111110", type=str,
                    help="Positions to consider when comparing expected and actual barcode sequences (1=Yes, 0=No)")
parser.add_argument("--reverse", action="store_true", help="Instead, output reads which fall outside the mismatch threshold")
parser.add_argument("--no_trim", action="store_true", help="Instead, output entries without trimming the adapter sequence")
parser.add_argument("--trim_other_end", action="store_true", help="In addition, examine the other end of the read for a barcode. Will not remove partial barcodes")


def main(args=None, sysStdin=None, printPrefix="PRODUSE-TRIM\t"):
    IUPACCodeDict = {
        'A': {'A'},  # Adenine
        'C': {'C'},  # Cytosine
        'G': {'G'},  # Guanine
        'T': {'T', 'U'},  # Thymine
        'U': {'T', 'U'},  # Uracil
        'R': {'A', 'G'},  # A or G
        'Y': {'C', 'T'},  # C or T
        'S': {'G', 'C'},  # G or C
        'W': {'A', 'T'},  # A or T
        'K': {'G', 'T'},  # G or T
        'M': {'A', 'C'},  # A or C
        'B': {'C', 'G', 'T'},  # C or G or T
        'D': {'A', 'G', 'T'},  # A or G or T
        'H': {'A', 'C', 'T'},  # A or C or T
        'V': {'A', 'C', 'G'},  # A or C or G
        'N': {'A', 'C', 'G', 'T'},  # any base
        '.': {'.', '-'},  # gap
        '-': {'.', '-'}  # gap
    }

    compliment = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'U': 'A',
        'R': 'Y',
        'Y': 'R',
        'S': 'S',
        'W': 'W',
        'K': 'M',
        'M': 'K',
        'B': 'V',
        'D': 'H',
        'H': 'D',
        'V': 'B',
        'N': 'N'
    }

    if args is None:
        if sysStdin is None:
            args = parser.parse_args()
            args = vars(args)
        else:
            args = parser.parse_args(sysStdin)
            args = vars(args)

    # If a config file is specified, parse the arguments from that
    if args["config"] is not None:
        config = ConfigObj(args["config"])
        try:
            for argument, parameter in config["trim"].items():
                if argument in args and args[argument] is None:  # Aka this argument is used by trim, and a parameter was not provided at the command line
                    args[argument] = parameter
        except KeyError:  # i.e. there is no section named "trim" in the config file
            sys.stderr.write(
                "ERROR: Unable to locate a section named \'trim\' in the config file \'%s\'\n" % (args["config"]))
            exit(1)

    # Finally, re-parse the arguments to ensure are required arguments are specified
    args = validateArgs(args)

    # Sanity check to ensure that the input FASTQ files are not the same
    if os.path.realpath(args["input"][0]) == os.path.realpath(args["input"][1]):
        parser.error("The input files specified are the same file!")

    # Check that the barcode sequence and mask are the same length
    if len(args["barcode_sequence"]) != len(args["barcode_position"]):
        parser.error("The barcode sequence and barcode mask must be the same length")

    # Just in case someone used lowercase
    args["barcode_sequence"] = args["barcode_sequence"].upper()

    # Make sure the user used IUPAC bases in the barcode sequence
    for base in args["barcode_sequence"]:
        if base not in IUPACCodeDict:
            parser.error("Unrecognized IUPAC base: %s" % base)

    barcodeLength = len(args["barcode_position"])
    barcodeIndexes = tuple(i for i in range(0, barcodeLength) if args["barcode_position"][i] == "1")
    barcodeBases = tuple(args["barcode_sequence"][i] for i in barcodeIndexes)

    # Obtain the reverse compliment of the barcode, if the other end of the read is to be examined
    if args["trim_other_end"]:
        reverseIndexes = tuple(barcodeLength - i - 1 for i in barcodeIndexes)
        reverseIndexes = reverseIndexes[::-1]
        reverseBarcode = barcodeBases[::-1]
        reverseBarcode = tuple(compliment[x] for x in reverseBarcode)

    # Open the input and output fastq files for reading/writing
    # Determine if the input files are gzipped, and open appropriately
    if args["input"][0].split(".")[-1] == "gz":
        f1in = gzip.open(args["input"][0], "rt")
    else:
        f1in = open(args["input"][0])

    if args["input"][1].split(".")[-1] == "gz":
        f2in = gzip.open(args["input"][1], "rt")
    else:
        f2in = open(args["input"][1])

    # Open outputs
    if args["output"][0].endswith(".gz"):
        f1out = gzip.open(args["output"][0], "wt")
    else:
        f1out = open(args["output"][0], "w")
    if args["output"][1].endswith(".gz"):
        f2out = gzip.open(args["output"][1], "wt")
    else:
        f2out = open(args["output"][1], "w")

    # Begin processing the FASTQ files
    # Read 1
    r1name = None
    r1seq = None
    r1strand = None
    r1qual = None

    # Read 2
    r2name = None
    r2seq = None
    r2strand = None
    r2qual = None

    # Status messages
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Starting...\n"]))
    count = 0
    discard = 0
    discardWarning = False

    for line1, line2 in zip(f1in, f2in):

        # Load the entire record before processing it
        if r1name is None:  # i.e. we have not recieved a record name yet, lets add them here
            r1name = line1
            r2name = line2
        elif r1seq is None:
            r1seq = line1
            r2seq = line2
        elif r1strand is None:  # i.e. we have not recieved the record's strand info yet
            r1strand = line1
            r2strand = line2
        else:  # These lines contain the quality scores. Once we recieve those, we can begin processing the read pair for real
            r1qual = line1
            r2qual = line2

            # First, obtain the barcode sequences
            barcode1 = r1seq[:barcodeLength]
            barcode2 = r2seq[:barcodeLength]

            # Append the family barcode to the read name
            # We will specify to BWA that this is a tag later
            familyBarcode = barcode1 + barcode2

            # Remove the Illumina tag from the reads, since BWA's -C option will add that tag to the BAM file as well
            # which will corrupt the output SAM
            r1name = r1name.split(" ")[0].rstrip() + "\tOX:Z:" + familyBarcode + os.linesep
            r2name = r2name.split(" ")[0].rstrip() + "\tOX:Z:" + familyBarcode + os.linesep

            # Second, determine if the barcode sequences fall within the mismatch theshold
            # Obtain the positions that will actually be compared
            b1pos = tuple(barcode1[i] for i in barcodeIndexes)
            b2pos = tuple(barcode2[i] for i in barcodeIndexes)

            # Next,compare these positions to the barcode sequence
            mismatch = 0
            for b1, b2, bRef in zip(b1pos, b2pos, barcodeBases):
                if b1 not in IUPACCodeDict[bRef]:
                    mismatch += 1
                if b2 not in IUPACCodeDict[bRef]:
                    mismatch += 1

            # If the mismatch is greater than the mismatch tolerated, "discard" this read by emptying the sequence
            # Unless the option is specified to keep reads which fall outside the mismatch threshold. In which case,
            # do the opposite
            if (not args["reverse"] and mismatch > args["max_mismatch"]) or (args["reverse"] and mismatch < args["max_mismatch"]):
                discard += 2
            else:
                # Otherwise, the barcode of this read pair is within thresholds
                # Trim the sequence and quality scores, unless the user has specified otherwise
                if not args["no_trim"]:
                    r1seq = r1seq[barcodeLength:]
                    r2seq = r2seq[barcodeLength:]
                    r1qual = r1qual[barcodeLength:]
                    r2qual = r2qual[barcodeLength:]


                # If we are checking the other end of the read for the presence of a barcode, do so now
                # Note that we do NOT discard reads here
                if args["trim_other_end"] and not args["no_trim"]:
                    # Yeah yeah, I know I'm duplicating code. But it's not that much code
                    # Obtain the candidate reverse barcode sequence
                    barcode1 =r1seq.rstrip()[-barcodeLength:]  # Ignore the newline
                    barcode2 = r1seq.rstrip()[-barcodeLength:]
                    familyBarcode = barcode1 + barcode2
                    b1pos = tuple(barcode1[i] for i in reverseIndexes)
                    b2pos = tuple(barcode2[i] for i in reverseIndexes)
                    # Next,compare these positions to the barcode sequence
                    mismatch = 0
                    for b1, b2, bRef in zip(b1pos, b2pos, reverseBarcode):
                        if b1 not in IUPACCodeDict[bRef]:
                            mismatch += 1
                        if b2 not in IUPACCodeDict[bRef]:
                            mismatch += 1

                    # If the trailing sequence is within the mismatch, assume that is is a barcode, and trim it
                    if mismatch < args["max_mismatch"]:
                        r1seq = r1seq.rstrip()[:-barcodeLength] + os.linesep
                        r2seq = r2seq.rstrip()[:-barcodeLength] + os.linesep
                        r1qual = r1qual.rstrip()[:-barcodeLength] + os.linesep
                        r2qual = r2qual.rstrip()[:-barcodeLength] + os.linesep

                # Finally, output these reads
                f1out.write(r1name)
                f1out.write(r1seq)
                f1out.write(r1strand)
                f1out.write(r1qual)
                f2out.write(r2name)
                f2out.write(r2seq)
                f2out.write(r2strand)
                f2out.write(r2qual)

            # Finally, clear the previous read pair (the one that was just written out) from the buffer
            r1name = None
            r1strand = None
            r1seq = None
            r1qual = None
            r2name = None
            r2strand = None
            r2seq = None
            r2qual = None

            count += 2
            if count % 100000 == 0:

                sys.stderr.write("\t".join([printPrefix, time.strftime('%X'),
                                            "Discard Rate:" + str(float(discard) / float(count) * 100)[:3] + "%",
                                            "Count:" + str(count) + "\n"]))
                # Check to see if the discard rate is excessive (>70%). If so, the barcode is probably incorrect
                if float(discard) / float(count) > 0.70 and not discardWarning:
                    sys.stderr.write("WARNING: The read discard rate is excessively high. Are you sure the barcode sequence is correct?\n")
                    sys.stderr.write("You can check the barcode using \'adapter_predict\'\n")
                    discardWarning = True

    f1in.close()
    f2in.close()
    f1out.close()
    f2out.close()

    # Final messages
    if count % 100000 != 0:
        sys.stderr.write("\t".join([printPrefix, time.strftime('%X'),
                                    "Discard Rate:" + str(float(discard) / float(count) * 100)[:3] + "%",
                                    "Count:" + str(count) + "\n"]))
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Trimming Complete\n"]))


if __name__ == "__main__":
    main()
