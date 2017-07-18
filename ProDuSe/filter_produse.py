#! /usr/bin/env python

# USAGE:
# 	See filter_produse.py -h for details
#
# Description:
# 	Filters ProDuSe variant calls
# 	TODO: Expand this
#
# AUTHOR:
# 	Christopher Rushton (ckrushto@sfu.ca)

import configargparse
import configparser
import os
import numpy as np
import time
import sys
import pysam
import subprocess

"""
Processes command line arguments
"""
parser = configargparse.ArgumentParser(description="Filters ProDuSe variant calls")
parser.add_argument("-c", "--config", is_config_file=True, type=lambda x: isValidFile(x, parser), help="Optional configuration file, which can provide any of the input arguments.")
parser.add_argument("-i", "--input", type=lambda x: isValidFile(x, parser), required=True, help="Input ProDuSe vcf file, output of \'snv.py\'")
parser.add_argument("-m", "--molecule_stats", type=lambda x: isValidFile(x, parser), required=True, help="Input molecule stats file, output of \'snv.py\'")
parser.add_argument("-o", "--output", required=True, help="Output VCF file name")
parser.add_argument("-sb", "--strand_bias_threshold", type=float, default=0.05, help="Strand bias p-value threshold, below which vairants will be discarded")
parser.add_argument("-ss", "--strong_singleton_threshold", default=1, type=int, help="Base threshold for strong singleton bases (SP, SN) [Default: %(default)s]")
parser.add_argument("-sd", "--strong_duplex_threshold", default=1, type=int, help="Base threshold for strong duplex bases (DPN, DPn, DpN) [Default: %(default)s]")
parser.add_argument("-wt", "--weak_base_threshold", default=2, type=int, help="Weak supported base count theshold [Default: %(default)s]")
parser.add_argument("-sv", "--allow_single_stranded", action="store_true", default=False, help="Allow variants with only single stranded support")
parser.add_argument("-md", "--min_depth", type=int, default=2, help="Minimum depth threshold [Default: %(default)s]")
parser.add_argument("-nb", "--normal_bam", type=lambda x: isValidFile(x, parser), default=None, metavar="BAM", help="A BAM file coresponding to a matched normal sample")
parser.add_argument("-nv", "--normal_vaf", default=0.05, type=float, metavar="FLOAT", help="VAF threshold for the normal sample, above which variants will be called as germline [Default: %(default)s]")
parser.add_argument("-r", "--reference", type=lambda x: isValidFile(x, parser), metavar="FASTA", help="Reference genome, in FASTA format. Required if a normal BAM file is supplied")
parser.add_argument("-fl", "--filter_log", metavar="FILE", help="A log file to explain the thresholds used for each variant, and why variants failed filters")
parser.add_argument("-g", "--germline_output", metavar="FILE", help="If a matched normal was supplied, an output file for germline variants")


def isValidFile(file, parser):
    """
    Checks to ensure the specified file exists, and throws an error if it does not

    Args:
        file: A filepath
        parser: An argparse.ArgumentParser() object. Used to throw an exception if the file does not exist

    Returns:
        file: The file itself

    Raises:
        parser.error(): An ArgumentParser.error() object, thrown if the file does not exist
    """
    if os.path.isfile(file):
        return file
    else:
        parser.error("Unable to locate %s" % (file))


def getMoleculeIndex(header):
    """
    Identifies the column of each molecule type from the header

    Molecule types: DPN, DPn, DpN, Dpn, SN, SP, Sn, Sp

    Args:
        header: A tab-delinated string coresponding to the file header line
    Returns:
        moleculeIndexes: A dictionary listing molecule:Index
    """
    moleculeIndexes = {}

    headerColumns = header.split()
    i = 0
    for column in headerColumns:
        moleculeIndexes[column] = i
        i += 1

    return moleculeIndexes


def setThresholds(statsFile, strongSingletonBaseline, strongDuplexBaseline, weakBaseline):
    """
    Sets filtering thresolds based upon molecule counts

    TODO: Expand explination

    Args:
        statsFile: Path to a file listing molecule abundance and total counts for each position
        strongBaseline: An int representing the minimum number of strong molecules required to call a variant as real
        weakBaseline: An int representing the minimum number of strong molecules required to call a variant as real

    Return:
        molecThresholds: A dictionary containing the slope and offset (intersect) for each molecule

    """

    # Lists to store the counts and total depth for each type of molecule
    moleculeCounts = {"DPN": [], "DPn": [], "DpN": [], "Dpn": [], "SN": [], "SP": [], "Sn": [], "Sp": []}

    # Obtain molecule abundance and total depth at each locus
    with open(statsFile) as f:
        for line in f:
            # Skip header lines
            if line[0] == "#":
                if line[0:2] == "##":
                    continue
                columnIndex = getMoleculeIndex(line)
                continue

            molCounts = line.split()
            try:
                # Obtain the total molecule abundance at this locus
                for molecule in moleculeCounts:
                    totalMol = int(molCounts[columnIndex[molecule + "_TOTAL"]])
                    altMol = int(molCounts[columnIndex[molecule + "_ALT"]])

                    # Weight this locus based upon overall depth
                    for i in range(0, totalMol // 8):
                        moleculeCounts[molecule].append([totalMol, altMol])
            except KeyError as x:
                raise TypeError("The header of %s is malformed or missing" % (statsFile))

    molecThreshold = {}

    # Here, determine the minimum TOTAL number of molecules required to accurately filter something with this type of molecule
    for molecule in moleculeCounts:

        lineCoord = sorted(moleculeCounts[molecule], key=lambda x: x[0])
        x = list(lineCoord[x][0] for x in range(0, len(lineCoord)))
        y = list(lineCoord[x][1] for x in range(0, len(lineCoord)))
        # If there are no molecules of this type across the entire capture space, they will not be used for filtering
        if (len(x) == 0 or np.mean(x) == 0) and (len(y) == 0 or np.mean(y) == 0):
            slope = 0
            offset = 0
        else:
            slope, offset = np.polyfit(x, y, 1)
        # If the offset is below 0 (which doesn't really make sense biologically), set it to 0
        if offset < 0:
            offset = 0
        # Determine the type of molecule, and thus the minimum molecule threshold
        if "P" in molecule or "N" in molecule:  # Strong molecule

            if "D" in molecule:
                molecThreshold[molecule] = {"slope": slope, "offset": offset + strongDuplexBaseline}
            else:
                molecThreshold[molecule] = {"slope": slope, "offset": offset + strongSingletonBaseline}
        else:
            molecThreshold[molecule] = {"slope": slope, "offset": offset + weakBaseline}

    return molecThreshold


def processFields(line):
    """
    Creates a dictionary from each field in the supplied string

    Args:
        line: A string containing semi-colon deliminated fields, in the format of NAME=E1,E2,E3,E4,E5 etc
    Returns:
        fieldDict: A dictionary listing the fields, in the format of fieldDict[NAME] = [E1, E2, E3, E4, E5]
    """

    fieldDict = {}
    lineFields = line.split(";")
    for field in lineFields:
        fieldName, fieldElements = field.split("=")
        elementList = fieldElements.split(",")
        fieldDict[fieldName] = elementList
    return fieldDict


def vafFromPileup(chrom, locus, bamFile, refAllele, altAlleles):
    """
    Calculates the VAF at a given locus

    Using samtools mpileup, calculate the number of reference and alternate alleles at this position.

    Args:
        chrom: A string representing the chromsome name
        locus: An integer representing the locus of the variant (1-indexed)
        bamFile: A pysam.AlignmentFile object coresponding to the BAM file of interest. Must be indexed
        refAllele: A string containing the reference allele at this position (i.e. A,C,G,T)
        altAlleles: A list containing the alternate alleles at this position (i.e. ["A"])
    """

    """
    for locus in bamFile.pileup(chrom, locus-1, locus):
        for sequence in locus.pileups:
            print(sequence.alignment.query_sequence[sequence.query_position])
    """
    alleleCounts = {"A":0, "C":0, "G":0, "T":0, "I":0, "D":0}
    vafs = {}

    DEVNULL = open(os.devnull, "w")
    # Generate a pileup of this region
    region = chrom + ":" + str(locus) + "-" + str(locus)
    samtoolsCom = ["samtools", "mpileup", "-r", region, bamFile]
    pileup = subprocess.Popen(samtoolsCom, stdout=subprocess.PIPE, stderr=DEVNULL)  # DANGER ZONE: Supressing error messages

    for line in pileup.stdout:
        line = line.decode("utf-8")
        # Determine allele counts at this position
        alleles = line.split("\t")[4]
        indel = False
        for allele in alleles:  # Loop through each letter in the pileup
            # If we previously dealt with an indel, loop until we reach the next allele
            # Format here will be +[0-9][ACGTactg], so skip past the letter
            if indel:
                if isinstance(allele, int):
                    continue
                else:
                    indel = False
                    continue
            # Ignore info on the start and end of a read, as this is not useful at this point
            if allele == "$" or allele == "^":
                continue
            # If there is an insertion or deletion at this position
            if allele == "+":
                alleleCounts["I"] += 1
                indel = True
            elif allele == "-":
                alleleCounts["D"] += 1
                indel = True
            # Process nucleotide alleles
            elif allele == "T" or allele == "t":
                alleleCounts["T"] += 1
            elif allele == "A" or allele == "a":
                alleleCounts["A"] += 1
            elif allele == "G" or allele == "g":
                alleleCounts["G"] += 1
            elif allele == "C" or allele == "c":
                alleleCounts["C"] += 1
            # Ignore everything else

    refCount = float(alleleCounts[refAllele])
    for altAllele in altAlleles:
        altCount = float(alleleCounts[altAllele])
        if refCount == 0 and altCount > 0:
            vafs[altAllele] = 1
        elif refCount == 0 and altCount == 0:
            vafs[altAllele] = 0
        else:
            vafs[altAllele] = altCount / (refCount + altCount)

    return vafs


def runFilter(vcfFile, thresholds, outFile, minDepth=0, strandBiasThresh=0.05, allowSingle=False, logFile=None, matchedNormal=None, vafThreshold=0.05, germOut=None, noReformat=False):
    """
    Filters variants at each locus based upon molecule counts at that locus

    TODO: Expand this

    Args:
        vcfFile: Path to ProDuSe raw variants file, output of produse snv
        threholds: A dictionary listing slope and offset for each molecule type
        outFile: Path to use for the output VCF
        minDepth: Minimum depth required at a locus to call a variant
        strandBias: Strand bias p-value, blow which variants will be discarded
        allowSingle: Allow variants without duplex support to pass the filter
        logFile: A filepath string to explain why each variant passed or failed filters
        matchedNormal: A filepath string to a BAM file coresponding to the matched normal
        vafThresold: A float indicating the germline VAF threshold, above which variants called in the matched normal will be considered germline
        germOut: A string to a filepath for germline variants
        noReformat: Boolean indicating if the output VCF format should be identical to the input VCF format
    """

    printPrefix = "PRODUSE-FILTER\t"

    # Don't create a germline vcf file if no matched normal was supplied
    if not matchedNormal:
        germOut = None

    # Counters for the number of alternate alleles passing and failing filters
    passedVar = 0
    failedVar = 0

    baseToIndex = {"A": 0, "C": 1, "G": 2, "T": 3}
    molTypes = ["DPN", "DPn", "DpN", "SN", "SP", "Dpn", "Sn", "Sp"]
    with open(vcfFile) as f, open(outFile, "w") as o, (open(logFile, "w") if logFile is not None else open(os.devnull, "w")) as log, (open(germOut, "w") if germOut is not None else open(os.devnull, "w")) as germVCF:
        for line in f:
            # Ignore header lines
            if line[0] == "#":
                o.write(line)
                germVCF.write(line)
                continue

            log.write("\n")
            chrom, start, ID, refAllele, altAlleles = line.split()[0:5]

            passingPosAlleles = []
            passingNegAlleles = []

            log.write("Variant " + chrom + ":" + start + "\n")

            if "," in altAlleles:
                altAlleles = altAlleles.split(",")
            else:
                altAlleles = [altAlleles]

            # Filter this variant out if it is germline (if a matched normal was provided)
            if matchedNormal:
                normVaf = vafFromPileup(chrom, int(start), matchedNormal, refAllele, altAlleles)
            else:
                normVaf = None

            infoCol = line.split()[7]
            infoFields = processFields(infoCol)

            # Filter for minimum depth
            totalDepth = sum(int(infoFields["MC"][x]) for x in range(0, 4))
            if totalDepth < minDepth:
                log.write("Alternate Allele \'%s\': FAILED\n" % (",".join(altAlleles)))
                log.write("Cause: Variant depth %i is less than minimum depth %i\n" % (totalDepth, minDepth))
                continue

            for altAllele in altAlleles:

                log.write("Processing Allele \'%s\'\n" % (altAllele))

                # Remove germline variants
                if normVaf and normVaf[altAllele] > vafThreshold:
                    log.write("Alternate Allele \'%s\': FAILED\n" % (altAllele))
                    log.write("Cause: Allele is germline: Germline VAF: %f, VAF Threshold: %f\n" % (normVaf[altAllele], vafThreshold))
                    germVCF.write(line)  # Rare case: May output duplicates
                    continue

                # Here, I am implementing a "3-strike rule". If there is a lot (>100) of molecules at this locus, and none of them indicate a variant,
                # there is likely no variant at this position
                # This should remove a lot of noisy variants in high-duplicate samples
                strikes = []

                # Ignore alleles with significant strand bias
                strandBias = float(infoFields["StrBiasP"][baseToIndex[altAllele]])
                if strandBias <= strandBiasThresh:
                    log.write("Alternate Allele \'%s\': FAILED\n" % (altAllele))
                    log.write("Cause: Allele strand bias %f is below threshold %f\n" % (strandBias, strandBiasThresh))
                    continue

                for molecule in molTypes:

                    supportsVariant = False
                    if len(strikes) >= 3:
                        log.write("Alternate Allele \'%s\': FAILED\n" % (altAllele))
                        log.write("Cause: Insufficent variant support in molecules %s\n" % (",".join(strikes)))
                        break

                    # Set the threshold based upon total molecule depth
                    totalMolecules = float(infoFields[molecule][baseToIndex[refAllele]])
                    for x in altAlleles:
                        totalMolecules += float(infoFields[molecule][baseToIndex[x]])

                    threshold = int(round(thresholds[molecule]["offset"])) + thresholds[molecule]["slope"] * totalMolecules
                    log.write("%s Molecule Threshold: %i: " % (molecule, threshold))

                    altMolecules = float(infoFields[molecule][baseToIndex[altAllele]])
                    log.write("Alternate Allele Molecules: %f: Status: " % (altMolecules))

                    if altMolecules > 0:
                        supportsVariant = True

                    if altMolecules > threshold:

                        log.write("PASS")
                        # What strand supports this variant?
                        if "p" in molecule or "P" in molecule:
                            if altAllele not in passingPosAlleles:
                                passingPosAlleles.append(altAllele)
                        if "n" in molecule or "N" in molecule:
                            if altAllele not in passingNegAlleles:
                                passingNegAlleles.append(altAllele)
                    else:
                        log.write("FAIL")

                    if not supportsVariant and sum(int(infoFields[molecule][x]) for x in range(0, 4)) > threshold * 2:
                        # Do not strike if the same molecule class has already been given a strike
                        if molecule == "DpN" and "DPn" in strikes:
                            log.write("\n")
                            continue
                        elif molecule == "SP" and "SN" in strikes:
                            log.write("\n")
                            continue
                        elif molecule == "Sp" and "Sn" in strikes:  # Not really necessary
                            log.write("\n")
                            continue
                        else:
                            strikes.append(molecule)
                            log.write(": Strike added\n")  # This should be documented somewhere
                    else:
                        log.write("\n")

                if len(passingPosAlleles) == 0 and len(passingNegAlleles) == 0:
                    log.write("Alternate Allele \'%s\': FAILED\n" % (altAllele))
                    log.write("Cause: No base types passed filters\n")

            # Store alternate alleles which pass filters
            passingAltAlleles = []
            if not allowSingle:
                for allele in altAlleles:
                    if allele in passingPosAlleles and allele in passingNegAlleles:
                        log.write("Alternate Allele \'%s\': PASSED\n" % (altAllele))
                        passingAltAlleles.append(allele)
                    elif allele in passingPosAlleles or allele in passingNegAlleles:
                        log.write("Alternate Allele \'%s\': FAILED\n" % (altAllele))
                        log.write("Cause: Lack of Duplex support\n")
            else:
                for allele in altAlleles:
                    if allele in passingPosAlleles or allele in passingNegAlleles:
                        passingAltAlleles.append(allele)
                        log.write("Alternate Allele \'%s\': PASSED\n" % (altAllele))

            if passingAltAlleles:
                # Replace the existing VAF with the VAF of alt alleles which passed filters
                altMolecules = 0
                for allele in passingAltAlleles:
                    altMolecules += float(infoFields["MC"][baseToIndex[allele]])
                totalMolecules = sum(float(infoFields["MC"][x]) for x in range(0, 4))
                vaf = altMolecules / totalMolecules
                newInfoCol = infoCol.split(";")[:-1]
                newInfoCol.append("VAF=" + str(vaf))
                outInfo = ";".join(newInfoCol)
                outLine = line.replace(infoCol, outInfo)
                origAltAlleles = line.split()[4]
                newAltAlleles = ",".join(passingAltAlleles)
                outLine = outLine.replace(origAltAlleles, newAltAlleles)

                passedVar += 1

                o.write(outLine)
            else:
                failedVar += 1

        sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Filtering Complete\n"]))

        # Correct Plurals!
        if passedVar == 1:
            passedStatusLine = "1 variant passed filters"
        else:
            passedStatusLine = "\f variants passed filters" % (passedVar)
        sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), passedStatusLine + "\n"]))

        if failedVar == 1:
            failedStatusLine = "1 variant failed filters"
        else:
            failedStatusLine = "\f variants failed filters" % (failedVar)
        sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), failedStatusLine + "\n"]))


def checkRef(bamFile, vcfFile):
    """
    Checks to ensure that the reference genome path supplied was used to align the supplied BAM file and for the vcf
    variants. If it was not, the user is given a warning
    
    Args:
        bamFile: A file path to the BAM file of interest
        vcfFile: A file path to the VCF file of interest
        refGenome: A file path to the reference genome
    """

    # Check the VCF file
    vcfContigs = []
    with open(vcfFile) as f:
        line = f.readline()
        while line[0] == "#":
            if "contig" in line:
                vcfContig = line.split("=")[2].split(",")[0]
                vcfContigs.append(vcfContig)

            line = f.readline()

    # Check BAM file:
    openBAM = pysam.AlignmentFile(bamFile)
    for type, value in openBAM.header.items():

        if type == "SQ":  # Contig names

            for contig in value:
                contigName = contig["SN"]
                if contigName not in (vcfContigs):
                    sys.stderr.write("WARNING: Contigs used by BAM file %s do not match contigs used by VCF file %s\n" % (
                    bamFile, vcfFile))
                    sys.stderr.write("Ensure the normal sample was aligned to the same reference genome that ProDuSe snv was run with\n")
                    break


def main(args=None):
    """

    """
    # Obtain arguments
    if args is None:
        args = parser.parse_args()
    elif args.config:
        # Since configargparse does not parse commands from the config file if they are passed as argument here
        # They must be parsed manually
        cmdArgs = vars(args)
        cmdArgs["normal_bam"] = None
        config = configparser.ConfigParser()
        config.read(args.config)
        configOptions = config.options("config")
        for option in configOptions:
            param = config.get("config", option)
            # Convert arguments that are lists into an actual list
            if param[0] == "[" and param[-1] == "]":
                param = param[1:-1].split(",")

            # WARNING: Command line arguments will be SUPERSCEEDED BY CONFIG FILE ARGUMENTS
            cmdArgs[option] = param

    printPrefix = "PRODUSE-FILTER\t"
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Starting...\n"]))

    # If a normal sample was supplied, process it
    if args.normal_bam:

        # Check to ensure that a reference genome was provided

        # If the supplied normal files are fastq files
        if args.normal_bam.endswith(".bam") or args.normal_bam.endswith(".sam"):
            # Create an index if one does not already exist
            if not os.path.exists(args.normal_bam + ".bai"):
                sys.stdout.write("\t".join([printPrefix, time.strftime('%X'), "Indexing Normal\n"]))
                pysam.index(args.normal_bam)
            checkRef(args.normal_bam, args.input)
        else:
            sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "ERROR: Normal File must be in either BAM or SAM format\n"]))
            sys.exit(1)

    thresholds = setThresholds(args.molecule_stats, args.strong_singleton_threshold, args.strong_duplex_threshold, args.weak_base_threshold)
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Thresholds Set\n"]))
    runFilter(args.input, thresholds, args.output, int(args.min_depth), float(args.strand_bias_threshold), args.allow_single_stranded, args.filter_log, args.normal_bam, float(args.normal_vaf), args.germline_output, args.no_reformat)


if __name__ == "__main__":
    main()
