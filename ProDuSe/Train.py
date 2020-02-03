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
import matplotlib
matplotlib.use('Agg')
import seaborn
import matplotlib.pyplot as plt
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
    parser.add_argument("-b", "--bam", metavar="BAM", nargs="+", type=lambda x: isValidFile(x, parser),
                        help="Input post-collapse or post-clipoverlap BAM file")
    parser.add_argument("-v", "--validations", metavar="VCF", nargs="+", type=lambda x: isValidFile(x, parser),
                        help="Input VCF file listing validated variants.")
#    parser.add_argument("-m", "--mappability", metavar="WIG", type=lambda x: isValidFile(x, parser), required=True,
#                        help="UCSC genome mappability file in wiggle format")
    parser.add_argument("-o", "--output", metavar="PICKLE", required=True, help="Output pickle file, containing the filtering classifier")
    parser.add_argument("-r", "--reference", metavar="FASTA", required=True, type=lambda x: isValidFile(x, parser),
                        help="Reference Genome, in FASTA format")
    parser.add_argument("-t", "--targets", metavar="BED", nargs="+", type=lambda x: isValidFile(x, parser),
                        help="A BED file containing regions in which to restrict variant calling")
    miscArgs = parser.add_argument_group(description="Miscellaneous arguments")
    miscArgs.add_argument("--ignore_vcf", metavar="VCF", type=lambda x: isValidFile(x, parser),
                        help="An optional VCF file listing variants that are to be excluded from filter training")
    miscArgs.add_argument("-j", "--jobs", metavar="INT", type=int, default=1,
                          help="Number of chromosomes to process simultaneously")
    miscArgs.add_argument("--input_files", metavar="TSV", type=lambda x: isValidFile(x, parser),
                          help="An input tab-deliniated file listing \'validations1.vcf input1.bam\'. May be used instead of \'-v\' and \'-b\'. May also specify \'-t\'.")
    miscArgs.add_argument("--true_stats", metavar="TSV", type=str,
                          help="An output file containing all true variant status used to train the filter")
    miscArgs.add_argument("--false_stats", metavar="TSV", type=str,
                          help="An optional output file containing all false variant status used to train the filter")
    miscArgs.add_argument("-d", "--plot_dir", metavar="DIR",
                          help="An output directory for density plots visualizing the various characteristics")

    args = parser.parse_args(listArgs)
    args = vars(args)

    # Check input
    if (not args["bam"] or not args["validations"]) and not args["input_files"]:
        raise parser.error("The following arguments are required: -b/--bam AND -v/--validations OR --input_files")

    return args

def processSampleMulti(args):
    """
    Unpacks processSample arguments

    :param args: A list of arguments to be unpacked
    :return: Twi tuples which contain variant stats for true and false variants
    """

    return processSample(*args)

def processSample(bamFile, refGenome, targets, contig, trueVariants, ignoreVariants, printPrefix):

    # Parse a temporary sample name
    baseSName = os.path.basename(bamFile)
    sampleName = baseSName.split(".")[0]

    # Run the pileup to identify candidate variants
    pileup = Call.PileupEngine(bamFile, refGenome, targets, printPrefix=printPrefix)

    trueVarStats = []
    falseVarStats = []

    #if "_" not in contig:  # Don't print status messages for minor contigs
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName, "Processing " + contig + os.linesep]))

    # Find candidate variants (same as Call)
    pileup.generatePileup(chrom=contig)

    #if "_" not in contig:  # Don't print status messages for minor contigs
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName, contig + ": Classifying Candidate Variants" + os.linesep]))
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

                posStats = list(pileup.varToFilteringStats(pos, maxAltAllele, chrom, location))
                # Add the location and sample nameof this variant to the end of the position stats. We will remove
                # these before before filter training
                posStats.append(chrom + ":" + str(location))
                posStats.append(sampleName)
                if chrom in trueVariants and location in trueVariants[chrom] and maxAltAllele in trueVariants[chrom][
                    location]:
                    trueVarStats.append(posStats)
                else:
                    falseVarStats.append(posStats)

        chromToDelete.append(chrom)
    for chrom in chromToDelete:
        del pileup.candidateVar[chrom]

    # If there are more false variants than true variants (which is usually the case), subset them
    # Take double the number of false variants as true variants, since there are far more artifact variants
    # than real variants, and this will help us measure all the possible different artifact variants
    if len(falseVarStats) > len(trueVarStats) * 2:
        falseVarStats = random.sample(falseVarStats, len(trueVarStats) * 2)

    return trueVarStats, falseVarStats

def parseInputFiles(args):
    """
    Parses input files from a specified config file

    Examines the config file specified by args["input_files"], and parses out the input BAM files, VCF files, and
    capture space files (if specified)
    The input file format is assumed to be in the format of (tab-deliniated):
    sample1.bam    sample1.vcf
    :param args: A dictionary listing input parameters
    :return: args: A dictionary storing the input parameters
    """
    bamFiles = []
    vcfFiles = []
    targetFiles = []

    with open(args["input_files"]) as f:
        for line in f:
            # Remove line endings
            line = line.rstrip("\n").rstrip("\r")
            cols = line.split("\t")
            try:
                vcfFile = cols[0]
                bamFile = cols[1]
                try:
                    targetFile = cols[2]
                    # Check that this file exists
                    if not os.path.exists(targetFile):
                        raise FileNotFoundError("Unable to locate \'%s\'" % targetFile)
                except IndexError:  # i.e. No target space was specified for this sample. That's ok, just set it to None
                    targetFile = None
                # Check that the remaining input files exist
                if not os.path.exists(bamFile):
                    raise FileNotFoundError("Unable to locate \'%s\'" % bamFile)
                if not os.path.exists(vcfFile):
                    raise FileNotFoundError("Unable to locate \'%s\'" % vcfFile)
                # Now that we know these files actuall exist, save them
                bamFiles.append(bamFile)
                vcfFiles.append(vcfFile)
                targetFiles.append(targetFile)
            except IndexError as e:  # i.e we are missing one or more columns in the input file
                raise TypeError("--input_files must be a tab-deliniated file listing \'variants.vcf    bamFile.bam\'") from e

    # Now that we have parsed all the input files, modify the arguments dictionary to include them
    args["bam"] = bamFiles
    args["validations"] = vcfFiles
    args["targets"] = targetFiles

    return args


def generatePlots(outDir, trueVarStats, falseVarStats, featureNames):

    # As position and sample are stored at the end of the variant stats, by iterating over the feature names
    # we will automatically skip plotting those two elements
    for i in range(0, len(featureNames)):
        # What is this feature?
        name = featureNames[i]
        # Don't plot artifact variants, since that is a boolean
        if name == "C->A mutation":
            continue
        # Parse the statistics of true variants and false variants for this feature
        tStats = list(x[i] for x in trueVarStats)
        fStats = list(x[i] for x in falseVarStats)
        outFileName = outDir + os.sep + name + ".png"
        seaborn.distplot(tStats, hist=False, color="b", label="True Variants")
        hist = seaborn.distplot(fStats, hist=False, color="r", label="False Variants", axlabel=name)
        plot = hist.get_figure()
        plot.savefig(outFileName)
        plt.clf()


parser = argparse.ArgumentParser(description="Creates a new random forest filtering model trained using the supplied variants")
parser.add_argument("-c", "--config", metavar="INI", type=lambda x: isValidFile(x, parser), help="An optional configuration file which can provide one or more arguments")
parser.add_argument("-b", "--bam", metavar="BAM", type=lambda x: isValidFile(x, parser), nargs="+", help="One or more post-collapse BAM files. BAM files must be specified in the same order as \'--validations\'")
parser.add_argument("-o", "--output", metavar="PICKLE", help="Output pickle file, containing the trained filter")
parser.add_argument("-r", "--reference", metavar="FASTA", type=lambda x: isValidFile(x, parser), help="Reference Genome, in FASTA format")
parser.add_argument("-t", "--targets", metavar="BED", type=lambda x: isValidFile(x, parser, allowNone=True), nargs="+", help="One or more BED files containing regions in which to restrict variant calling. Must be specified in the same order as \'--bam\' (use \'None\' for no BED file)")
parser.add_argument("-v", "--validations", metavar="VCF", type=lambda x: isValidFile(x, parser), nargs="+", help="One or more VCF files listing validated variants.")
# parser.add_argument("-m", "--mappability", metavar="WIG", type=lambda x: isValidFile(x, parser), help="UCSC genome mappability file in wiggle format")
miscArgs = parser.add_argument_group(description="Miscellaneous arguments")
miscArgs.add_argument("--ignore_vcf", metavar="VCF", type=lambda x: isValidFile(x, parser, allowNone=True), help="An optional VCF file listing variants that are to be excluded from filter training")
miscArgs.add_argument("-j", "--jobs", metavar="INT", type=int, help="Number of chromosomes to process simultaneously")
miscArgs.add_argument("--input_files", metavar="TSV", type=lambda x:isValidFile(x, parser), help="An input tab-deliniated file listing \'validations1.vcf input1.bam\'. May be used instead of \'-v\' and \'-b\'. May also specify targets.")
miscArgs.add_argument("--true_stats", metavar="TSV", type=str, help="An output file containing all true variant status used to train the filter")
miscArgs.add_argument("--false_stats", metavar="TSV", type=str, help="An optional output file containing all false variant status used to train the filter")
miscArgs.add_argument("-d", "--plot_dir", metavar="DIR", help="An output directory for density plots visualizing the various characteristics")

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
    # Check that the args generated from the config file are valid
    args = validateArgs(args)
    # Was a input config specified? If so, lets parse the input files from that
    if args["input_files"]:
        args = parseInputFiles(args)
    if args["targets"] is None:
        args["targets"] = [None] * len(args["bam"])
    elif len(args["targets"]) == 1:  # i.e. all samples have the same capture space
        args["targets"] = args["targets"] * len(args["bam"])

    # Parse contig names
    contigNames = Fasta(args["reference"]).records.keys()
    # Aggregate the statistics to generate the filter
    trueVarStats = []
    falseVarStats = []

    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Starting...\n"]))

    # Load validated variants
    i = 0
    ignoreVariants = loadVariants(args["ignore_vcf"])
    for validations in args["validations"]:
        trueVariants = loadVariants(validations)

        sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Processing Sample \'%s\'\n" % os.path.basename(args["bam"][i])]))

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
                varStats = processPool.imap_unordered(processSampleMulti, multiprocArgs)
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
    varFeatures = ["Total_Supporting_Molecules", "Strand_Bias", "Mean_Base_Quality", "Base_Quality_Bias",
                        "Mean_Mapping_Quality", "Mapping_Quality_Bias", "Total_Depth", "Homopolymer", "Nearby_Mismatch_Proportion",
                        "Mean_Read_Mismatches", "Read_Mismatch_Bias", "Mean_Family_Size", "Family_Size_Bias",
                        "Duplex_Counts", "C->A mutation"]

    # Dump the stats used to train the filters to the specified output files
    with open(args["true_stats"] if args["true_stats"] is not None else os.devnull, "w") as t,\
            open(args["false_stats"] if args["false_stats"] is not None else os.devnull, "w") as f:

        # Write the file headers
        outHeader = []
        outHeader.extend(varFeatures)
        outHeader.extend(["Position", "Sample"])
        f.write("\t".join(outHeader))
        t.write("\t".join(outHeader))
        f.write(os.linesep)
        t.write(os.linesep)
        for line in trueVarStats:
            t.write("\t".join(list(str(x) for x in line)) + os.linesep)
        for line in falseVarStats:
            f.write("\t".join(list(str(x) for x in line)) + os.linesep)

    # Merge the false and true variant stats
    varStats = list(x[:-2] for x in trueVarStats)  # Strip off position and sample name
    realOrNot = [0] * len(trueVarStats)
    varStats.extend(list(x[:-2] for x in falseVarStats))
    realOrNot.extend([1]*len(falseVarStats))

    filter = RandomForestClassifier(n_estimators=100, max_features=None)
    filter.fit(varStats, realOrNot)

    # Write the trained classifier out to the file
    with open(args["output"], "w+b") as o:
        pickle._dump(filter, o)

    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Filter trained and saved in \'%s\'\n" % args["output"]]))

    # Print a summary of the feature weights
    sys.stdout.write("Filter feature importances" + os.linesep)
    for feature, importance in zip(varFeatures, filter.feature_importances_):
        sys.stdout.write(feature + "\t" + str(importance) + os.linesep)

    # Generate plots of these characteristics
    if args["plot_dir"]:
        sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Generating Figures..." + os.linesep]))
        generatePlots(args["plot_dir"], trueVarStats, falseVarStats, varFeatures)

if __name__ == "__main__":
    main()
