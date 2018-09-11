#!/usr/bin/env python

# Import modules
import pickle
import argparse
import os
import sys
import time
import random
from sklearn.ensemble import RandomForestClassifier
from configobj import ConfigObj
from pyfaidx import Fasta
import multiprocessing
try:
    import Call
    import ProDuSeExceptions as pe
except ImportError:
    from ProDuSe import Call
    from ProDuSe import ProDuSeExceptions as pe

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
    Load variants from the specified VCF file, and store them in to dictionaries
    :param inFile: A string containing a filepath to a VCF file listing mutation status
    :return:
    """

    variants = {}
    varCount = 0

    if inFile is None:  # i.e. no VCF file was specified
        return {}

    with open(inFile) as f:
        # Sanity check that the input file is a VCF file
        line = f.readline()
        if "fileformat=VCF" not in line:
            raise pe.InvalidInputException("Input file \'%s\' does not meet VCF specifications" % inFile)

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

            # Add this variant to the coresponding validations
            varCount += 1
            if chrom not in variants:
                variants[chrom] = {}
            if position not in variants[chrom]:
                variants[chrom][position] = set()
            variants[chrom][position].add(alt)

    return variants


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
                        help="Input VCF file listing validated variants.")
    parser.add_argument("--ignore_vcfs", metavar="VCF", type=lambda x: isValidFile(x, parser), nargs="+",
                        help="One or more optional VCF files listing variants that are to be excluded from filter training")
    parser.add_argument("-o", "--output", metavar="PICKLE", required=True, help="Output pickle file, containing the filtering classifier")
    parser.add_argument("-r", "--reference", metavar="FASTA", required=True, type=lambda x: isValidFile(x, parser),
                        help="Reference Genome, in FASTA format")
    parser.add_argument("-t", "--targets", metavar="BED", nargs="+", type=lambda x: isValidFile(x, parser),
                        help="A BED file containing regions in which to restrict variant calling")
    parser.add_argument("-j", "--jobs", metavar="INT", type=int, default=1, help="Number of chromosomes to process simultaneously")
    parser.add_argument("--true_stats", metavar="TSV", type=lambda x: isValidFile(x, parser),
                        help="An optional output file containing all true variant status used to train the filter")
    parser.add_argument("--false_stats", metavar="TSV", type=lambda x: isValidFile(x, parser),
                        help="An optional output file containing all false variant status used to train the filter")
    parser.add_argument("-d", "--plot_dir", metavar="DIR",
                        help="An output directory for density plots visualizing the various characteristics")

    args = parser.parse_args(listArgs)
    return vars(args)

def processSample(bamFile, refGenome, targets, contig, trueVariants, ignoreVariants, printPrefix):

    # Run the pileup to identify candidate variants
    pileup = Call.PileupEngine(bamFile, refGenome, targets, printPrefix=printPrefix)

    trueVarStats = []
    falseVarStats = []

    #if "_" not in contig:  # Don't print status messages for minor contigs
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Processing " + contig + "\n"]))

    # Find candidate variants (same as Call)
    pileup.generatePileup(chrom=contig)

    #if "_" not in contig:  # Don't print status messages for minor contigs
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), contig + ": Classifying Candidate Variants\n"]))
    chromToDelete = []
    for chrom, positions in pileup.candidateVar.items():

        for location, pos in positions.items():
            if pos.summarizeVariant():  # i.e. this variant passes the minimum depth filter

                # Only examine the major allele
                maxWeight = 0
                maxAltAllele = None
                for allele in pos.altAlleles.keys():
                    weight = pos.strandCounts[allele]
                    if maxAltAllele is None or weight > maxWeight:
                        maxAltAllele = allele
                        maxWeight = weight
                if chrom in ignoreVariants and location in ignoreVariants[chrom] and maxAltAllele in \
                        ignoreVariants[chrom][location]:
                    continue  # Skip this variant

                # Skip variants where the reference allele is supported by too few reads
                if pos.strandCounts[pos.ref] < 4:
                    continue
                # Obtain the stats that we wish to filter upon

                posStats = pileup.varToFilteringStats(pos, maxAltAllele)

                if chrom in trueVariants and location in trueVariants[chrom] and maxAltAllele in trueVariants[chrom][
                    location]:
                    trueVarStats.append(posStats)
                else:
                    falseVarStats.append(posStats)
        chromToDelete.append(chrom)
    for chrom in chromToDelete:
        del pileup.candidateVar[chrom]

    return trueVarStats, falseVarStats


parser = argparse.ArgumentParser(description="Creates a new random forest filtering model trained using the supplied variants")
parser.add_argument("-c", "--config", metavar="INI", type=lambda x: isValidFile(x, parser), help="An optional configuration file which can provide one or more arguments")
parser.add_argument("-b", "--bam", metavar="BAM", type=lambda x: isValidFile(x, parser), nargs="+", help="One or more post-collapse BAM files. BAM files must be specified in the same order as \'--validations\'")
parser.add_argument("-v", "--validations", metavar="VCF", type=lambda x: isValidFile(x, parser), nargs="+", help="One or more VCF files listing validated variants.")
parser.add_argument("--ignore_vcfs", metavar="VCF", type=lambda x: isValidFile(x, parser, allowNone=True), nargs="+", help="One or more optional VCF files listing variants that are to be excluded from filter training")
parser.add_argument("-o", "--output", metavar="PICKLE", help="Output pickle file, containing the trained filter")
parser.add_argument("-r", "--reference", metavar="FASTA", type=lambda x: isValidFile(x, parser), help="Reference Genome, in FASTA format")
parser.add_argument("-t", "--targets", metavar="BED", type=lambda x: isValidFile(x, parser, allowNone=True), nargs="+", help="One or more BED files containing regions in which to restrict variant calling. Must be specified in the same order as \'--bam\' (use \'None\' for no BED file)")
parser.add_argument("-j", "--jobs", metavar="INT", type=int, help="Number of chromosomes to process simultaneously")
parser.add_argument("--true_stats", metavar="TSV", type=lambda x: isValidFile(x, parser), help="An output file containing all true variant status used to train the filter")
parser.add_argument("--false_stats", metavar="TSV", type=lambda x: isValidFile(x, parser), help="An optional output file containing all false variant status used to train the filter")
parser.add_argument("-d", "--plot_dir", metavar="DIR", help="An output directory for density plots visualizing the various characteristics")

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
    if args["targets"] is None:
        args["targets"] = [None] * len(args["bam"])
    if args["ignore_vcfs"] is None:
        args["ignore_vcfs"] = [None] * len(args["bam"])

    # Parse contig names
    contigNames = Fasta(args["reference"]).records.keys()
    # Aggregate the statistics to generate the filter
    trueVarStats = []
    falseVarStats = []

    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Starting...\n"]))

    # Load validated variants
    i = 0
    for validations in args["validations"]:
        trueVariants = loadVariants(validations)
        ignoreVariants = loadVariants(args["ignore_vcfs"][i])

        sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Processing Sample \'%s\'\n" % args["bam"][i].split("/")[-1]]))

        # Run single threaded
        if args["jobs"] == 1:

            # Process each chromosome seperately to reduce memory footprint
            for contig in contigNames:

                sTrueVarStats, sFalseVarStats = processSample(args["bam"][i], args["reference"], args["targets"][i], contig, trueVariants, ignoreVariants, printPrefix)
                trueVarStats.extend(sTrueVarStats)
                falseVarStats.extend(sFalseVarStats)
        else:  # Run Multi threaded
            # Parallelize by chromosome
            if args["jobs"] == 0:  # i.e. use as many threads as possible
                threads = None  # The multiprocessing module will chose the number of threads
            else:
                threads = min(len(contigNames), args["jobs"], os.cpu_count())

            # Prepare arguments for pool
            multiprocArgs = list([args["bam"][i], args["reference"], args["targets"][i], x, trueVariants, ignoreVariants, printPrefix] for x in contigNames)
            processPool = multiprocessing.Pool(processes=threads)
            try:
                varStats = processPool.starmap_async(processSample, multiprocArgs).get()
                processPool.close()
                processPool.join()
            except KeyboardInterrupt as e:
                processPool.terminate()
                processPool.join()
                raise e

            # Extract true and false variant stats
            for vStats in varStats:
                trueVarStats.extend(vStats[0])
                falseVarStats.extend(vStats[1])

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

    filter = RandomForestClassifier(n_estimators=50, max_features="auto")
    filter.fit(varStats, realOrNot)

    # Write the trained classifier out to the file
    with open(args["output"], "w+b") as o:
        pickle._dump(filter, o)

    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Filter trained and saved in \'%s\'\n" % args["output"]]))

if __name__ == "__main__":
    main()
