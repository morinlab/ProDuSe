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

"""
Processes command line arguments
"""
parser = configargparse.ArgumentParser(description="Filters ProDuSe variant calls")
parser.add_argument("-c", "--config", is_config_file=True, type=lambda x: isValidFile(x, parser), help="Optional configuration file, which can provide any of the input arguments.")
parser.add_argument("-i", "--input", type=lambda x: isValidFile(x, parser), required=True, help="Input ProDuSe vcf file, output of \'snv.py\'")
parser.add_argument("-m", "--molecule_stats", type=lambda x: isValidFile(x, parser), required=True, help="Input molecule stats file, output of \'snv.py\'")
parser.add_argument("-o", "--output", required=True, help="Output VCF file name")
parser.add_argument("-sb", "--strand_bias_threshold", default=0.05, help="Strand bias p-value threshold, below which vairants will be discarded")
parser.add_argument("-st", "--strong_base_threshold", default=1, type=int, help="Strong supported base count threshold [Default: %(default)s]")
parser.add_argument("-wt", "--weak_base_threshold", default=2, type=int, help="Weak supported base count theshold [Default: %(default)s]")
parser.add_argument("-ss", "--allow_single_stranded", action="store_true", default=False, help="Allow variants with only single stranded support")


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


def setThresholds(statsFile, strongBaseline, weakBaseline):
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
			molecThreshold[molecule] = {"slope": slope, "offset": offset + strongBaseline}
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


def runFilter(vcfFile, thresholds, outFile, strandBiasThresh=0.05, allow_single=False):
	"""
	Filters variants at each locus based upon molecule counts at that locus

	TODO: Expand this

	Args:
		vcfFile: Path to ProDuSe raw variants file, output of produse snv
		threholds: A dictionary listing slope and offset for each molecule type
		outFile: Path to use for the output VCF
		require_dual: Require variant support on both strands to call a variant as real
	"""
	printPrefix = "PRODUSE-FILTER\t"
	sys.stdout.write("\t".join([printPrefix, time.strftime('%X'), "Starting...\n"]))

	baseToIndex = {"A": 0, "C": 1, "G": 2, "T": 3}
	molTypes = ["DPN", "DPn", "DpN", "SN", "SP", "Dpn", "Sn", "Sp"]
	with open(vcfFile) as f, open(outFile, "w") as o:
		for line in f:
			# Ignore header lines
			if line[0] == "#":
				o.write(line)
				continue

			passingPosAlleles = []
			passingNegAlleles = []

			refAllele, altAlleles = line.split()[3:5]
			if "," in altAlleles:
				altAlleles = altAlleles.split(",")
			else:
				altAlleles = [altAlleles]

			infoCol = line.split()[7]
			infoFields = processFields(infoCol)

			# Here, I am implementing a "3-strike rule". If there is a lot (>100) of molecules at this locus, and none of them indicate a variant,
			# there is likely no variant at this position
			# This should remove a lot of noisy variants in high-duplicate samples
			strikes = []
			for molecule in molTypes:

				supportsVariant = False
				if len(strikes) >= 3:
					break
				# Count the number of ref and alt molecules at this locus
				totalMolecules = float(infoFields[molecule][baseToIndex[refAllele]])
				for altAllele in altAlleles:
					totalMolecules += float(infoFields[molecule][baseToIndex[altAllele]])

				# Set the threshold based upon total molecule depth
				threshold = int(round(thresholds[molecule]["offset"])) + thresholds[molecule]["slope"] * totalMolecules

				# Process each possible alternate allele individually
				for altAllele in altAlleles:

					# Throw ignore alleles with significant strand bias
					if float(infoFields["StrBiasP"][baseToIndex[altAllele]]) <= strandBiasThresh:
						continue

					altMolecules = float(infoFields[molecule][baseToIndex[altAllele]])
					if altMolecules > threshold:

						supportsVariant = True
						# What strand supports this variant?
						if "p" in molecule or "P" in molecule:
							if altAllele not in passingPosAlleles:
								passingPosAlleles.append(altAllele)
						if "n" in molecule or "N" in molecule:
							if altAllele not in passingNegAlleles:
								passingNegAlleles.append(altAllele)
				if not supportsVariant and sum(int(infoFields[molecule][x]) for x in range(0, 4)) > threshold * 2:
					# Do not strike if the same molecule class has already been given a strike
					if molecule == "DpN" and "DPn" in strikes:
						continue
					elif molecule == "SP" and "SN" in strikes:
						continue
					elif molecule == "Sp" and "Sn" in strikes:
						continue
					strikes.append(molecule)

			passingAltAlleles = []
			if not allow_single:
				for allele in ["A", "C", "G", "T"]:
					if allele in passingPosAlleles and allele in passingNegAlleles:
						passingAltAlleles.append(allele)
			else:
				for allele in ["A", "C", "G", "T"]:
					if allele in passingPosAlleles or allele in passingNegAlleles:
						passingAltAlleles.append(allele)

			if passingAltAlleles:
				# Calculate the VAF of these alleles
				altMolecules = 0
				for allele in passingAltAlleles:
					altMolecules += float(infoFields["MC"][baseToIndex[allele]])
				totalMolecules = sum(float(infoFields["MC"][x]) for x in range(0, 4))
				vaf = altMolecules / totalMolecules
				outInfo = ";".join([infoCol, "VAF=" + str(vaf)])
				outLine = line.replace(infoCol, outInfo)
				origAltAlleles = line.split()[4]
				newAltAlleles = ",".join(passingAltAlleles)
				outLine = outLine.replace(origAltAlleles, newAltAlleles)
				o.write(outLine)


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

    thresholds = setThresholds(args.molecule_stats, args.strong_base_threshold, args.weak_base_threshold)
    runFilter(args.input, thresholds, args.output, args.strand_bias_threshold, args.allow_single_stranded)


if __name__ == "__main__":
	main()
