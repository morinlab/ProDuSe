#!/usr/bin/env python

# Import modules
import pickle
import argparse
import os
import sys
import time
import random
import seaborn
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from configobj import ConfigObj
try:
    from ProDuSe import Call
except ImportError:  # i.e. ProDuSe is not installed
    import Call


def isValidFile(file, parser, allowNone=False):
    """
    Checks to ensure the specified file path is valid

    :param file: A string containing a filepath
    :param parser: An argparse.ArgumentParser() object
    :param allowNone: If "None" should be treated as a valid file
    :return: file, if the specified file exists
    :raises parser.error: If the specified file does not exist
    """

    if os.path.exists(file):
        return file
    elif allowNone and file.upper() == "NONE":
        return file
    parser.error("Unable to locate \'%s\': No such file or directory." % file)


def loadVariants(inFile):
    """
    Load variants from the specified VCF file, and store them in to dictionaries based upon their validation status
    :param inFile: A string containing a filepath to a VCF file listing mutation status
    :return:
    """

    trueVariants = {}
    unclassifiedVariants = {}
    skippedVar = 0
    validatedVar = 0

    with open(inFile) as f:
        # Sanity check that the input file is a VCF file
        line = f.readline()
        if "fileformat=VCF" not in line:
            raise TypeError("Input file \'%s\' does not meet VCF specifications" % inFile)

        for line in f:
            line = line.rstrip()
            # Ignore header lines. We will assume that the remainder of the file complies with VCF specifications
            if line.startswith("#"):
                continue
            # Parse the chromosome, position, alternate allele, and validation status from this variant entry
            cols = line.split("\t")
            chrom = cols[0]
            position = int(cols[1]) - 1  # Since pysam uses 0-based indexing
            alt = cols[4]
            info = cols[7]
            validated = None
            infoFields = info.split(";")
            for field in infoFields:
                name, value = field.split("=")
                if name == "VALIDATED":
                    validated = None if value.upper() == "UNKNOWN" else value.upper() == "TRUE"
                    break
            if validated is None:  # No validation status was specified for this variant
                skippedVar += 1
                continue
            # Add this variant to the coresponding validations
            validatedVar += 1
            if validated is True:
                if chrom not in trueVariants:
                    trueVariants[chrom] = {}
                if position not in trueVariants[chrom]:
                    trueVariants[chrom][position] = set()
                trueVariants[chrom][position].add(alt)
            elif validated is None:
                if chrom not in unclassifiedVariants:
                    unclassifiedVariants[chrom] = {}
                if position not in unclassifiedVariants[chrom]:
                    unclassifiedVariants[chrom][position] = set()
                unclassifiedVariants[chrom][position].add(alt)

        if skippedVar > 0:
            sys.stderr.write("WARNING: Validation status was not specified for %s variants. These will be ignored.\n" % skippedVar)

    return trueVariants, unclassifiedVariants


def validateArgs(args):
    """
    Validates that the specified set of arguments are valid, and that all required arguments are set

    This is done here, as we can't check ahead of time if the arguments specified in a config file
    are valid until they have been processed

    """

    # Convert the dictionary into a list, to allow the arguments to be re-parsed
    listArgs = []
    for argument, parameter in args.items():

        # If this argument was not set, ignore it
        if parameter is None or parameter is False or parameter == "None" or parameter == "False":
            continue
        # Something was provided as an argument
        listArgs.append("--" + argument)

        # Ignore booleans, as we will re-add them when the arguments are re-parsed
        if parameter == "True" or parameter is True:
            continue
        # If the parameter is a list, we need to add each element separately
        if isinstance(parameter, list):
            for p in parameter:
                listArgs.append(str(p))
        else:
            listArgs.append(str(parameter))
    parser = argparse.ArgumentParser(
        description="Creates a new random forest filtering model based upon the specified validated variants")
    parser.add_argument("-c", "--config", metavar="INI", type=lambda x: isValidFile(x, parser),
                        help="An optional configuration file which can provide one or more arguments")
    parser.add_argument("-b", "--bam", metavar="BAM", required=True, nargs="+", type=lambda x: isValidFile(x, parser),
                        help="Input post-collapse or post-clipoverlap BAM file")
    parser.add_argument("-v", "--validations", metavar="VCF", required=True, nargs="+", type=lambda x: isValidFile(x, parser),
                        help="Input VCF file listing validated variants. Validation status must be specified using the INFO field \'VALIDATED\'")
    parser.add_argument("-o", "--output", metavar="PICKLE", required=True, help="Output pickle file, containing the filtering classifier")
    parser.add_argument("-r", "--reference", metavar="FASTA", required=True, type=lambda x: isValidFile(x, parser),
                        help="Reference Genome, in FASTA format")
    parser.add_argument("-t", "--targets", metavar="BED", nargs="+", type=lambda x: isValidFile(x, parser),
                        help="A BED file containing regions in which to restrict variant calling")
    parser.add_argument("--true_stats", metavar="TSV", type=lambda x: isValidFile(x, parser),
                        help="An optional file to dump stats coresponding to validated variants")
    parser.add_argument("--false_stats", metavar="TSV", type=lambda x: isValidFile(x, parser),
                        help="An optional file to dump stats coresponding to false variants that were chosen to train the filter")

    args = parser.parse_args(listArgs)
    return vars(args)


parser = argparse.ArgumentParser(description="Creates a new random forest filtering model based upon the specified validated variants")
parser.add_argument("-c", "--config", metavar="INI", type=lambda x: isValidFile(x, parser), help="An optional configuration file which can provide one or more arguments")
parser.add_argument("-b", "--bam", metavar="BAM", type=lambda x: isValidFile(x, parser), nargs="+", help="One or more post-collapse or post-clipoverlap BAM files. BAM files must be specified in the same order as \'--validations\'")
parser.add_argument("-v", "--validations", metavar="VCF", type=lambda x: isValidFile(x, parser), nargs="+", help="One or more VCF file listing validated variants. Validation status must be specified using the INFO field \'VALIDATED\'")
parser.add_argument("-o", "--output", metavar="PICKLE", help="Output pickle file, containing the filtering classifier")
parser.add_argument("-r", "--reference", metavar="FASTA", type=lambda x: isValidFile(x, parser), help="Reference Genome, in FASTA format")
parser.add_argument("-t", "--targets", metavar="BED", type=lambda x: isValidFile(x, parser, allowNone=True), nargs="+", help="One or more BED file(s) containing regions in which to restrict variant calling. Must be specified in the same order as \'--bam\' (use \'None\' for no BED file)")
parser.add_argument("--true_stats", metavar="TSV", type=lambda x: isValidFile(x, parser), help="An optional file to dump stats coresponding to validated variants")
parser.add_argument("--false_stats", metavar="TSV", type=lambda x: isValidFile(x, parser), help="An optional file to dump stats coresponding to false variants that were chosen to train the filter")

def main(args=None, sysStdin=None, printPrefix="PRODUSE-TRAIN\t"):
    if args is None:
        args = parser.parse_args(args=sysStdin)
        args = vars(args)

    # If a config file was specified, parse the arguments from
    if args["config"] is not None:
        config = ConfigObj(args["config"])
        try:
            for argument, parameter in config["train"].items():
                if argument in args and not args[argument]:  # Aka this argument is used by call, and a parameter was not provided at the command line
                    args[argument] = parameter
        except KeyError:  # i.e. there is no section named "call" in the config file
            sys.stderr.write(
                "ERROR: Unable to locate a section named \'train\' in the config file \'%s\'\n" % (args["config"]))
            exit(1)
    # Check that the args generated from the config file are validated
    args = validateArgs(args)

    # Sanity check input arguments
    if len(args["bam"]) != len(args["validations"]):
        sys.stderr.write("ERROR: The number of \'-v/--validations\' files must match the number of \'-b/--bam\' files")
        exit(1)
    if args["targets"] is None:
        args["targets"] = [None] * len(args["bam"])

    baseToIndex = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4}
    # Aggregate the statistics to generate the filter
    trueVarStats = []
    falseVarStats = []

    # Load validated variants
    i = 0
    for validations in args["validations"]:
        trueVariants, unclassifiedVariants = loadVariants(validations)

        # Run the pileup to identify candidate variants
        pileup = Call.PileupEngine(args["bam"][i], args["reference"], args["targets"][i], printPrefix=printPrefix)

        # Find candidate variants (same as Call)
        pileup.generatePileup()

        sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Classifying Candidate Variants\n"]))
        for chrom, positions in pileup.candidateVar.items():

            for location, pos in positions.items():

                try:
                    posStats, altBases = pos.collapsePosition()
                except KeyError:  # i.e. The reference based at this position is degenerate
                    continue
                # Is this a true or false variant
                for posStat, base in zip(posStats, altBases):

                    if chrom in unclassifiedVariants and location in unclassifiedVariants[chrom] and base in unclassifiedVariants[chrom][location]:
                        continue  # We can't train the filter based on this variant
                    if chrom in trueVariants and location in trueVariants[chrom] and base in trueVariants[chrom][location]:
                        trueVarStats.append(posStat)
                    else:
                        falseVarStats.append(posStat)
        i += 1
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Training Random Forest\n"]))

    # Subset the false variants so that there are an equal number of false and true variants
    # If there are more validated variants than false variants, do not subset
    if len(falseVarStats) > len(trueVarStats):
        falseVarStats = random.sample(falseVarStats, len(trueVarStats))

    # Dump the stats used to train the filters to the specified output files
    with open(args["true_stats"] if args["true_stats"] is not None else os.devnull, "w") as t,\
            open(args["false_stats"] if args["false_stats"] is not None else os.devnull, "w") as f:
        # Write the file headers
        f.write("\t".join(["quality", "family_size", "distance_to_read_end", "number_of_read_mismatches", "mapping_qual", "alt_count", "strand_bias_p", "family_in_duplex"]))
        for line in trueVarStats:
            t.write("\t".join(list(str(x) for x in line)) + os.linesep)
        for line in falseVarStats:
            f.write("\t".join(list(str(x) for x in line)) + os.linesep)

    # Merge the false and true variant stats
    varStats = trueVarStats
    realOrNot = [0] * len(trueVarStats)
    varStats.extend(falseVarStats)
    realOrNot.extend([1]*len(falseVarStats))

    filter = RandomForestClassifier(n_estimators=30, max_features="auto")
    filter.fit(varStats, realOrNot)

    # Write the trained classifier out to the file
    with open(args["output"], "w+b") as o:
        pickle._dump(filter, o)

    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Filter Trained\n"]))

if __name__ == "__main__":
    main()
