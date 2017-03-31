#! /usr/bin/env python


import configargparse
import os

"""
Processes command line arguments
"""
parser = configargparse.ArgumentParser(description="Filters ProDuSe variants calls")
parser.add_argument("-c", "--config", is_config_file=True, type=lambda x: isValidFile(x, parser), help="Optional configuration file, which can provide any of the input arguments.")
parser.add_argument("-dsv", "--totalvaf", type=float, default=0.05, help="Dual-strand VAF threshold [Default: %(default)s]")
parser.add_argument("-md", "--min_duplex", type=int, default=3, help="Minimum duplex support required to call a variant [Default: %(default)s]")
parser.add_argument("-mp", "--min_pos_strand", type=int, default=1, help="Minimum positive strand support required to call a varaint [Default: %(default)s]")
parser.add_argument("-mn", "--min_neg_strand", type=int, default=1, help="Minimum negative strand support required to call a variant [Default: %(default)s]")
parser.add_argument("-ms", "--min_singleton", type=int, default=3, help="Minimum singleton support required to call a variant [Default: %(default)s]")
parser.add_argument("-ww", "--weak_base_weight", type=lambda x: isValidWeight(x, parser), default=0.1, help="Weight of weak bases, relative to strong supported bases [Default: %(default)s]")
parser.add_argument("-i", "--input", type=lambda x: isValidFile(x, parser), required=True, help="Input ProDuSe vcf file, the output of \'snv.py\'")
parser.add_argument("-sb", "--strand_bias", type=float, default=0.05, help="Strand bias p-value threshold, below which variants will be ignored")
parser.add_argument("-o", "--output", required=True, help="Output VCF file name")


def isValidFile(file, parser):
	"""
    Checks to ensure the specified file exists, and throws an error if it does not

    Args:
        file: A filepath
        parser: An argparse.ArgumentParser() object. Used to throw an exception if the file does not exist

    Returns:
        type: The file itself

    Raises:
        parser.error(): An ArgumentParser.error() object, thrown if the file does not exist
    """
	if os.path.isfile(file):
		return file
	else:
		parser.error("Unable to locate %s" % (file))


def isValidWeight(weight, parser):
	"""
	Confirms that the weak supported base weight is a sensible value

	Args:
		weight: The user-provided value
		paser: A configargparse object
	Returns:
		weight: The user-provided value, if it is value
	Throws:
		parser.error(): A ConfigArgparser error object, if the value is not a float or is not reasonable
	"""

	try:
		weight = float(weight)
		if weight < 0 or weight > 1:
			raise parser.error("Weak base weight must be a float between 0 and 1")
		return weight
	except ValueError:
		raise parser.error("Weak base weight must be a float between 0 and 1")


def findFields(infoCol):
	"""
	Identifies the index of fields in the supplied string

	Args:
		infoCol: A list generated from the VCF info column, seperated by ";"
	Returns:
		fieldCols: A dictionary listing the name of each field and its index
	Raises
		TypeError: If one or more of the fields cannot be identified
	"""
	fieldCols = {"DPN":None, "DPn":None, "DpN":None, "Dpn":None, "SP":None, "Sp":None, "Sn":None, "SN":None, "StrBiasP":None, "MC":None}

	# Find the index of each field
	i = 0
	for field in infoCol:
		fieldName = field.split("=")[0]
		if fieldName in fieldCols:
			fieldCols[fieldName] = i

		i += 1

	# Check to ensure each field was found
	for name, index in fieldCols.items():
		if index is None:
			raise TypeError("Unable to locate %s in the info column of the input file" % (name))

	return fieldCols


def main(args=None):

	if args is None:
		args = parser.parse_args()

	with open(args.input) as f, open(args.output, "w") as o:

		headerFound = False
		fieldsFound = False

		for line in f:

			# Ignore header and info lines, just print them out
			if line[0] == "#":

				# Quick sanity check: For the header line, check to ensure the file meets VCF requirements
				if line[0:6] == "#CHROM":
					headerFound=True
					if line != "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n":
						raise TypeError("The input file is not in VCF format")
				o.write(line)

			else:
				# If no header was found, this file may not be in VCF format
				if not headerFound:
					raise TypeError("The input file is not in VCF format")

				# Identify normal allele and variant(s)
				refAllele = line.split("\t")[3]
				variants = line.split("\t")[4]
				if "," in variants:
					variants = variants.split(",")
				else:
					variants = list(variants)

				# Grab the VCF info column
				infoCol = line.split("\t")[7].split(";")

				# If this is the first content line, determine the locations of the relevant fiends in the INFO column
				if not fieldsFound:
					fieldIndex = findFields(infoCol)
					fieldsFound = True

				# Parse the fields of interest
				fieldValue = {}
				for field, index in fieldIndex.items():
					fieldValue[field] = infoCol[index].split("=")[1]

				# Time to actual filter variants
				valAlleles, vaf = filterVariant(
					refAllele,
					variants,
					fieldValue,
					args.totalvaf,
					args.min_duplex,
					args.min_pos_strand,
					args.min_neg_strand,
					args.min_singleton,
					args.strand_bias,
					args.weak_base_weight)

				if valAlleles is not None: # i.e. the variant passed filters
					vafOutput = "VAF=" + str(vaf) + "\n"
					# If a locus has multiple variants, and not all of those variants passed filters, adjust the output line
					if valAlleles != variants:
						outList = []
						outList.extend(line.split("\t")[0:4])
						outList.append(",".join(valAlleles))
						outList.extend(line.split("\t")[5:])
						outLine = "\t".join(outList)
					else:
						outLine = line
					o.write(",".join([outLine[:-1], vafOutput]))


def filterVariant(refAllele, altAllele, desc, minVaf, minDuplex, minPosStrand, minNegStrand, minSingleton, strandBiasThreshold, weakBaseWeight):
	"""
	Filters supplied variant in accordance with specified thresholds

	This variant fails filtering if:
	1) There is significant indication of strand bias (falls below the specified p-value threshold)
	2) The variant occurs exists below the specified VAF theshold
	3) Variant has too little duplex support AND the variant has too little singleton support on each strand and overall singleton support

	Note that locus properties are stored in the format of Name:A, C, G, T

	TODO: Treat a locus with multiple alternate alleles better
	Args:
		refAllele: Reference allele
		altAllele: A list of one or more alternate alleles
		desc: A dictionary of locus properies, including allele support, strand bias, and molecule counts
		minVaf: Minimum VAF threshold (See 2)
		minDuplex: Minimum number of duplex molecules supporting the variant (See 3)
		minPosStrand: Minimum number of singleton molecules from the positive strand supporting the variant (See 3)
		minNegStrand: Minimum number of singleton molecules from the negative strand supporting the variant (See 3)
		minSingleton: Minimum total number of singleton molecules supporting the variant (See 3)
		strandBiasThreshold: P-value threshold for strand bias (See 1)
		weakBaseWeight: Weight that should be given to weakly supporting bases relative to strong supporting bases (See 3)

	Returns:
		passedAlt: A list of alternate alleles that passed filters
		vaf: Variant allele fraction of the variant, if it passes filters. If it fails filters, 'None' is returned

	"""
	nucToIndex = {"A":0, "C":1, "G":2, "T":3}
	indexToNuc = {0:"A", 1:"C", 2:"G", 3:"T"}
	refIndex = nucToIndex[refAllele]
	altIndex = list(nucToIndex[alt] for alt in altAllele)

	# Step 1: Check Strand Bias
	# If multiple alternative alleles are present, each are checked for strand bias.
	# If an alternative allele exhibits strand bias, it is excluded from further analysis

	# I'm just going to play it safe and create a copy of the array
	altPassedStrandB = []

	for index in altIndex:
		strandP = float(desc["StrBiasP"].split(",")[index])
		if strandP > strandBiasThreshold:
			# Passes filter
			altPassedStrandB.append(index)

	# If no alleles passed strand bias filters, stop here
	# Otherwise, continue filtering with the remaining alleles
	if len(altPassedStrandB) == 0:
		return None, None
	else:
		altIndex = altPassedStrandB
		passedAlleles = list(indexToNuc[allele] for allele in altIndex)

	# Step 2: Check VAF
	# Calculate VAF
	refCount = float(desc["MC"].split(",")[refIndex])
	altCount = sum(float(desc["MC"].split(",")[i]) for i in altIndex)
	vaf = altCount / (altCount + refCount)

	# Run vaf filter
	if vaf < minVaf:
		return None, None

	# Step 3: Check for sufficient duplex or singleton support
	dupStrSupport = sum(int(desc["DPN"].split(",")[i]) for i in altIndex)
	dupWeakSupport = sum(int(desc["DPn"].split(",")[i]) for i in altIndex)
	dupWeakSupport += sum(int(desc["DpN"].split(",")[i]) for i in altIndex)
	dupWeakSupport += sum(int(desc["Dpn"].split(",")[i]) for i in altIndex)
	dupWeakSupport = dupWeakSupport * weakBaseWeight

	# If there is sufficient duplex support, stop here. The variant has passed all filters
	if dupStrSupport + dupWeakSupport >= minDuplex:
		return passedAlleles, vaf
	# If not, check singleton support
	else:
		negSupport = sum(int(desc["SP"].split(",")[i]) for i in altIndex)
		negSupport += sum(int(desc["Sp"].split(",")[i]) for i in altIndex) * weakBaseWeight
		posSupport = sum(int(desc["SN"].split(",")[i]) for i in altIndex)
		posSupport += sum(int(desc["Sn"].split(",")[i]) for i in altIndex) * weakBaseWeight
		if negSupport >= minNegStrand and posSupport >= minPosStrand and negSupport + posSupport >= minSingleton:
			return passedAlleles, vaf
		else:
			return None, None

if __name__ == "__main__":

	main()
