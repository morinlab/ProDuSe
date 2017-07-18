#!/usr/bin/env python

# USAGE:
# 	See python SplitMerge.py -h for details
#
# DESCRIPTION:
#  	Seperates paired-ends reads which have been merged together by Stitcher
#
# 	If a read was not merged, it is simply written to the output BAM file
# 	If a read was merged (stitched) together, the read tag written by Stitcher (i.e. 'XD')
# 	is examined to identify which bases in the sequence originated from which read.
# 	The merged (stitched) sequence is then divided into apropriate forward and reverse reads.
# 	Note that stitching is simplified, so in the case of complicated stitching (ex. 12F2R3F16S20R)
# 	some bases may be placed on one read when they originated from the opposite read. In addition,
# 	changes in read alignment caused by stitching are detected, and the cigar string of output reads
# 	is modified to include an INDEL at the apropriate location to restore proper read alignment
#
# AUTHOR:
# 	Christopher Rushton (ckrushto@sfu.ca)


import configargparse
import configparser
import sys
import os
import pysam
import copy
import re
import time

# Processes command line arguments
parser = configargparse.ArgumentParser(description="Splits reads merged by stitcher")
parser.add_argument("-c", "--config", required=False, is_config_file=True, type=lambda x: isValidFile(x, parser), help="Optional configuration file, which can provide any of the input arguments.")
parser.add_argument("-i", "--input", metavar="BAM", required=True, type=lambda x: isValidFile(x, parser), help="Input sorted BAM file, containing merged reads to be split")
parser.add_argument("-u", "--unstitched_input", metavar="BAM", required=True, type=lambda x: isValidFile(x, parser), help="Input sorted BAM file, coresponding to the original unstitched reads")
parser.add_argument("-o", "--output", required=True, metavar="BAM", help="Output BAM file name")


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
		parser.error("The file %s does not exist" % (file))


def checkInput(bamFile):
	"""
	Checks to ensure the provided BAM file contains reads stitched by stitcher

	Cycles through the BAM file until it finds a read with an "XD" tag. If no tag is found,
	it indicates that no reads have been stitched. This could be due to a lack of overlapping reads
	(very, very, very unlikely), or user error. In this case, an error is thrown

	Args:
		bamFile: An open pysam.AlignmentFile() object coresponding to the BAM file of interest

	"""
	# I included a counter to allow for very small BAM files (ex. test files) to be run sucessfully
	i = 0
	foundStitchedRead = False
	for read in bamFile.fetch(until_eof=True, multiple_iterators=True):
		if read.has_tag("XD"):
			foundStitchedRead = True
			break
		else:
			i += 1

	if not foundStitchedRead and i > 100:
		sys.stderr.write("ERROR: No stitched reads were found in the supplied --input BAM file\n")
		sys.stderr.write("Ensure Stitcher was run on the --input BAM file, and try again\n")
		sys.exit(1)


def checkSort(bamFile):
	"""
	Checks the provided BAM file to ensure it is sorted by read name

	Only checks the first 50 reads

	Prints an error message if the first 50 reads are not sorted

	Args:
		bamFile: An open pysam.AlignmentFile() object coresponding to the BAM file of interest
	"""

	head = bamFile.head(50)
	readList = list(read.query_name.split(":")[0] for read in head)
	if readList != sorted(readList):
		sys.stderr.write("ERROR: Input BAM files must be sorted by read name\n")
		sys.stderr.write("Try \'samtools sort -n file.bam > file.n_sorted.bam\'\n")
		exit(1)


def cigarToList(cigar):
	"""
	Converts a cigar string to a list, with each index indicating the base type

	Args:
		cigar: A cigar string (ex. 14M2S)
	Returns:
		cigarArray: A list, with each index coresponding to a cigar base
	"""
	cigarRegex = '[0-9]+[A-Z]'
	cigarComponents = re.findall(cigarRegex, cigar)
	cigarArray = []
	for component in cigarComponents:
		baseType = component[-1]
		baseLength = int(component[:-1])
		for i in range(0, baseLength):
			cigarArray.append(baseType)

	return cigarArray


def listToCigar(array):
	"""
	Converts a list to a cigar string

	I'm calling it an array, even though I know it's a list

	Args:
		array: A list coresponding to the cigar of interest
	Returns:
		cigar: A string representing the cigar sequence
	"""
	if array is None or len(array) == 0:
		return None
	cigar = ""
	previousChar = array[0]
	previousCount = 0
	for character in array:
		if previousChar != character:
			cigar += str(previousCount) + previousChar
			previousChar = character
			previousCount = 1
		else:
			previousCount += 1
	cigar += str(previousCount) + previousChar
	return cigar


def divideCigar(iCigar):
	"""

	"""

	# Find the sections coresponding to the forward and reverse reads
	softClippedRegex = '[0-9]+S'
	allClippedRegex = '[0-9]+[A-Z]'
	softClippedCigars = re.findall(softClippedRegex, iCigar)
	allCigars = re.findall(allClippedRegex, iCigar)
	forwardCigar = None
	reverseCigar = None
	sharedCigar = None

	if softClippedCigars[0] == iCigar:
		sharedCigar = iCigar
		return forwardCigar, sharedCigar, reverseCigar

	solutionFound = False
	# Finds the overlapping
	for cigar in softClippedCigars:
		for i in range(0, len(allCigars)):
			if allCigars[i] == cigar:
				if i + 1 >= len(allCigars):
					cigar2 = None
				else:
					cigar2 = "".join(allCigars[i + 1:])
				if i - 1 < 0:
					cigar1 = None
				else:
					cigar1 = "".join(allCigars[:i])

				if cigar1 is not None and cigar2 is not None:

					# The cigar has been split into discrete sections for each read
					if "F" in cigar1 and "R" not in cigar1 and "R" in cigar2 and "F" not in cigar2:
						forwardCigar = cigar1
						reverseCigar = cigar2
						sharedCigar = cigar
						solutionFound = True
						break
					elif "R" in cigar1 and "F" not in cigar1 and "F" in cigar2 and "R" not in cigar2:
						forwardCigar = cigar2
						reverseCigar = cigar1
						sharedCigar = cigar
						solutionFound = True
						break
					# The reverse read is completely encliped by the forward read:
					elif "F" in cigar2 and "R" not in cigar2 and "F" in cigar1 and "R" not in cigar1:
						forwardCigar = iCigar
						return forwardCigar, sharedCigar, reverseCigar
					# The forward read is completely encliped by the reverse read
					elif "F" not in cigar2 and "R" in cigar2 and "F" not in cigar1 and "R" in cigar1:
						reverseCigar = iCigar
						return forwardCigar, sharedCigar, reverseCigar
				elif cigar1 is None and cigar2 is not None:
					# The reverse read is a subset of the forward read
					if "F" in cigar2 and "R" not in cigar2:
						forwardCigar = cigar2
						sharedCigar = cigar
						solutionFound = True
						break
					# The forward read is a subset of the reverse read
					elif "R" in cigar2 and "F" not in cigar2:
						reverseCigar = cigar2
						sharedCigar = cigar
						solutionFound = True
						break
				elif cigar2 is None and cigar1 is not None:
					# The reverse read is a subset of the forward read
					if "F" in cigar1 and "R" not in cigar1:
						forwardCigar = cigar1
						sharedCigar = cigar
						solutionFound = True
						break
					# The forward read is a subset of the reverse read
					elif "R" in cigar1 and "F" not in cigar1:
						reverseCigar = cigar1
						sharedCigar = cigar
						solutionFound = True
						break
		if solutionFound:
			break

	if not solutionFound:

		# Simpist case: If the cigar starts and ends with bases from a single read, then the other read is fully eclipsed
		# In this case, just use that read
		iCigarList = cigarToList(iCigar)
		if iCigarList[0] == "F" and iCigarList[-1] == "F":
			forwardCigar = iCigar
			return forwardCigar, None, None
		elif iCigarList[0] == "R"and iCigarList[-1] == "R":
			reverseCigar = iCigar
			return None, None, reverseCigar

		# Find the largest soft-clipped region, and assume that coresponds to the overlap
		maxSize = 0
		maxCigar = None
		for cigar in softClippedCigars:
			cigarLen = len(cigar[:-1])
			if cigarLen > maxSize:
				maxSize = cigarLen
				maxCigar = cigar
		# Split into unique forward and reverse sections
		cigar1 = iCigar.split(maxCigar)[0]
		cigar2 = iCigar.split(maxCigar)[1]
		sharedCigar = maxCigar

		# Next, determine which base the starting and ending cigar start and end with, and use that as te cigar
		cigar1List = cigarToList(cigar1)
		if len(cigar1List) > 0:
			if cigar1List[0] == "F":
				forwardCigar = cigar1
			elif cigar1List[0] == "R":
				reverseCigar = cigar1
		cigar2List = cigarToList(cigar2)
		if len(cigar2List) > 0:
			if cigar2List[-1] == "F":
				forwardCigar = cigar1
			elif cigar2List[-1] == "R":
				reverseCigar = cigar2

	return forwardCigar, sharedCigar, reverseCigar


def effectiveLengthFromCigar(cigarList):
	"""
	Determines the actual length of the cigar sequence relative to the reference
	Ignores leading soft-clipped bases

	Args:
		cigarList: Cigar sequence
	Returns:
		length: Length of the sequence
	"""
	if cigarList is None:
		return 0
	length = 0
	leadingSoftClipped = True
	for char in cigarList:
		if char == "S" and leadingSoftClipped:
			continue
		if char != "I":
			length += 1
			leadingSoftClipped = False
	return length


def overlayCigars(iCigar, template, sequence):
	"""
	Determines the portion of the stitched cigar originating from each element in origCigar

	Args:
		stitchedCigar: The cigar to be split
		origCigar: The template of elemets to be used for splitting
	Returns:
		fowardCigar: A cigar string coresponding to the unique forward sequence
		forwardOffset: Offset of this subsequence from the start of the sequence
		sharedCigar: A cigar string coresponding to the shared sequence
		reverseCigar: A cigar string coresponding to the unique reverse sequence
		reverseOffset: Offset of this subsequence from the start of the seqeuence
	Raises:
		TypeError: If the length of the candidate and template cigars do not match
	"""

	iCigarList = cigarToList(iCigar)
	templateList = cigarToList(template)
	if len(iCigarList) != len(templateList):
		raise TypeError("The two cigar strings cannot be divided, as they are of different lengths")
	elif "S" not in templateList:
		raise TypeError("Cigar %s cannot be as a template, as no soft-clipped bases are indicated" % (templateList))

	# Package output in a dictionary
	sepSeqCig = {"fSeq":"", "fOffset":None, "fCigarList":None, "rSeq":"", "rOffset":None, "rCigarList":None, "sSeq":"", "sOffset":None, "sCigarList":None}

	forwardCigar, sharedCigar, reverseCigar = divideCigar(template)
	if not sharedCigar:
		if not forwardCigar:
			sepSeqCig["rOffset"] = 0
			sepSeqCig["rCigarList"] = iCigar
			sepSeqCig["rSeq"] = sequence
		elif not reverseCigar:
			sepSeqCig["fOffset"] = 0
			sepSeqCig["fCigarList"] = iCigar
			sepSeqCig["fSeq"] = sequence
		return sepSeqCig
	forwardOffset, fCigarList, fSeq = getCigarSeq(forwardCigar, template, iCigarList, sequence)
	reverseOffset, rCigarList, rSeq = getCigarSeq(reverseCigar, template, iCigarList, sequence)
	sharedOffset, sharedCigarList, sharedSeq = getCigarSeq(sharedCigar, template, iCigarList, sequence)

	# Edge case: If the bases unique to a read are entirely soft-clipped, append the sequence to another, as it is not informative alone
	if fCigarList is not None and "M" not in fCigarList and "I" not in fCigarList:
		if sharedCigarList is None:
			# If the forward cigar was second, append it to the reverse cigar
			if forwardOffset > reverseOffset:
				rCigarList.extend(rCigarList)
				rSeq = rSeq + fSeq
			else:
				fCigarList.extend(rCigarList)
				rCigarList = fCigarList
				reverseOffset = forwardOffset
				rSeq = fSeq + rSeq
		else:
			# If the forward cigar is after the shared cigar, append it to the merged cigar. Otherwise, prepend the shared cigar
			if forwardOffset > sharedOffset:
				sharedCigarList.extend(fCigarList)
				sharedSeq = sharedSeq + fSeq
			else:
				fCigarList.extend(sharedCigarList)
				sharedCigarList = fCigarList
				sharedOffset = forwardOffset
				sharedSeq = fSeq + sharedSeq
		fCigarList = None
		fSeq = ""
	if rCigarList is not None and "M" not in rCigarList and "I" not in rCigarList:
		if sharedCigarList is None:
			# If the reverse cigar was second, append it to the forward cigar
			if reverseOffset > forwardOffset:
				fCigarList.extend(rCigarList)
				fSeq = fSeq + rSeq
			else:
				rCigarList.extend(fCigarList)
				fCigarList = rCigarList
				forwardOffset = reverseOffset
				fSeq = rSeq + fSeq
		else:
			# If the reverse cigar is after the shared cigar, append it to the merged cigar. Otherwise, prepend the shared cigar
			if reverseOffset > sharedOffset:
				sharedCigarList.extend(rCigarList)
				sharedSeq = sharedSeq + rSeq
			else:
				rCigarList.extend(sharedCigarList)
				sharedCigarList = rCigarList
				sharedOffset = reverseOffset
				sharedSeq = rSeq + sharedSeq
		rCigarList = None
		rSeq = ""
	if sharedCigarList is not None and "M" not in sharedCigarList and "I" not in sharedCigarList:
		if fCigarList is not None:
			if sharedOffset > forwardOffset:
				fCigarList.extend(sharedCigarList)
				fSeq = fSeq + sharedSeq
			else:
				sharedCigarList.extend(fCigarList)
				fCigarList = sharedCigarList
				forwardOffset = sharedOffset
				fSeq = sharedSeq + fSeq
		else:
			if sharedOffset > reverseOffset:
				rCigarList.extend(sharedCigarList)
				rSeq = rSeq + sharedSeq
			else:
				sharedCigarList.extend(rCigarList)
				rCigarList = sharedCigarList
				reverseOffset = sharedOffset
				rSeq = sharedSeq + rSeq
		sharedCigarList = None
		sharedSeq = ""
	sepSeqCig["fCigarList"] = fCigarList
	sepSeqCig["rCigarList"] = rCigarList
	sepSeqCig["sCigarList"] = sharedCigarList
	sepSeqCig["fSeq"] = fSeq
	sepSeqCig["rSeq"] = rSeq
	sepSeqCig["sSeq"] = sharedSeq
	sepSeqCig["fOffset"] = forwardOffset
	sepSeqCig["rOffset"] = reverseOffset
	sepSeqCig["sOffset"] = sharedOffset

	return sepSeqCig


def findCigar(subCigar, origCigar):
	"""
	Determines the start and end location of the sub-cigar in the originalCigar
	"""

	startIndex = origCigar.find(subCigar)
	# Checks to ensure the cigar is a descrete unit (i.e. 3S, not 23S)
	# Regexes would be much slower
	while startIndex > 0 and origCigar[startIndex - 1].isdigit():
		startIndex = origCigar.find(subCigar, startIndex + 1)
	endIndex = startIndex + len(subCigar)
	startLocus = len(cigarToList(origCigar[:startIndex]))
	endLocus = len(cigarToList(origCigar[:endIndex]))

	return startLocus, endLocus


def getCigarSeq(subCig, template, origCigList, sequence):
	"""
	Finds the indexes of the sub-cigar string, and uses that to divide the sequence

	Determines the start and end index of the sub-cigar in the template cigar. Using this,
	the offset from the start of the sequence, and the sequence represented by this cigar,
	is determined. The proportion of the original cigar presented by this sub-cigar is also
	determined, and returned.

	Args:
		subCigar: A string coresponding to a sub-cigar of template cigar
		template: The original cigar string (string)
		origCigList: A list representing the read's actual cigar sequence
		seqnece: A nucleotide sequence (string)
	Returns:
		offset: The offset of this sub-cigar from the start of the cigar (int)
		cigOutList: The sequence of the original cigar represented by the sub-cigar (list)
		outSequence: The nucleotide sequence represented by the sub-cigar sequence (string)
	"""
	if subCig is None:
		return None, None, ""
	start, end = findCigar(subCig, template)
	cigOutList = origCigList[start:end]
	precCigarList = origCigList[:start]

	outsequence = ""
	offset = 0
	for char in precCigarList:
		if char == "D":
			continue
		offset += 1

	i = 0
	for char in cigOutList:
		if char == "D":
			continue
		outsequence += sequence[offset + i]
		i += 1
	return offset, cigOutList, outsequence


def getStartIndex(start, cigarList, origCigar=None):
	"""
	Calculates the start location of the next read, based upon leading soft-clipped bases and indels in the previous read

	Args:
		start: Expected start locus
		cigarList: Cigar sequence of the read in question, as a list
		previousCigar: Cigar string of the previous read
	Returns:
		start: Corrected start locus
	"""
	# Check the cigar sequence for leading soft-clipped bases, and offset as necessary
	if cigarList[0] == "S":
		if origCigar:
			origList = cigarToList(origCigar)
			i = 0
			try:
				while cigarList[i] == "S":
					if origList[i] != "S":
						start += 1
					i += 1
			except IndexError as e:
				print(cigar)
				raise e
		else:
			i = 0
			try:
				while cigarList[i] == "S":
					start += 1
					i += 1
			except IndexError as e:
				print(cigar)
				raise e
	return start


def createRead(iRead, flag, divStart, divEnd, tempLen, cigarSeq, mateStart, start):
	"""
	Creates a pysam.AlignmentSequence() object with the listed characteristics

	Creates a copy of iRead, and divides it accoriding to the specified paramters

	Args:
		iRead: A pysam.AlignmentSequence() object coresponding to the stitched (merged) read
		flag: An int corespoding to the Read flag (see https://samtools.github.io/hts-specs/SAMv1.pdf for details)
		divStart: Start locus of the divided sequence (int)
		divEnd: End locus of the divided sequence (int)
		tempLen: An int coresponding to template length or insert size
		cigarSeq: A string coresponding to the cigar string (see https://samtools.github.io/hts-specs/SAMv1.pdf for details)
		mateStart: An int coresponding to the mate read's start index
		start: An int coresponding to this read's start index
	Returns:
		oRead: A pysam.AlignmentSequence() object representing the output read
	"""

	oRead = copy.deepcopy(iRead)
	oRead.flag = flag
	oRead.reference_start = start
	oRead.query_sequence = iRead.query_sequence[divStart:divEnd]
	oRead.query_qualities = iRead.query_qualities[divStart:divEnd]
	oRead.template_length = tempLen
	oRead.cigarstring = cigarSeq
	oRead.next_reference_start = mateStart
	oRead.next_reference_id = oRead.reference_id
	return oRead


def findIndel(origForward, origReverse, fStitchedCigar, rStitchedCigar, offset):
	"""

	Determines the location of the indel in the stitched reads

	The predicted start and end locus for the original and stitched reads is offset,
	indicating that there is an indel in the original reads. This is likely represented
	by the locus where a read becomes soft-clipped. Here, the read containing relevant
	soft-clipping is identified, and an insertion or deletion is added to the stitcehd
	cigar string.

	Args:
		origForward: A pysam.AlignmentSequence() object coresponding to the forward read
		origReverse: A pysam.AlignmentSequence() object coresponding to the reverse read
		fStitchedCigar: A list coresponding to the cigar string of the forward de-stitched read
		rStitchedCigar: A list coresponding to the cigar string of the reverse de-stitched read
		offset: Difference between the end of the stitched cigar minus the the end of the original reads
	Returns:
		fOutCigar: A list coresponding to the cigar string of the forward de-stitched read,
					containing the INDEL
		rOutCigar: A list coresponding to the cigar string of the reverse de-stitched read,
					containing the INDEL
	TODO: Handle insertions
	"""
	if offset < 0:
		offsetType = "D"
	else:
		offsetType = "I"

	origFCigar = cigarToList(origForward.cigarstring)
	origRCigar = cigarToList(origReverse.cigarstring)

	# Determines which reads contain soft-clipped bases
	fSoftClipped = origFCigar[-1] == "S"
	rSoftClipped = origRCigar[0] == "S"
	fOutCigar = fStitchedCigar
	rOutCigar = rStitchedCigar

	# If the indel is on the forward read, add it to the de-stitched cigar
	if fSoftClipped and not rSoftClipped:
		for i in range(len(origFCigar) - 1, 0, -1):
			if origFCigar[i] != "S":
				i += 1
				break

		# Modifify the de-stitched cigar as apropriate
		if offsetType == "D":
			for j in range(0, abs(offset)):
				fOutCigar.insert(i, offsetType)
		else:
			pass
			# TODO: Handle insertions
			"""
			for j in range(0, abs(offset)):
				# Here, the insertion is too big to be explained by soft-clipping alone. Rescue failed
				if i + j >= len(fStitchedCigar):
					break
				fOutCigar[i + j] = offsetType
			"""
	# If the index is on the reverse read, add it to the de-stitched cigar
	elif not fSoftClipped and rSoftClipped:
		for i in range(0, len(origRCigar)):
			if origRCigar[i] != "S":
				i -= 1
				break
			offsetLength = len(origRCigar) - i
			offsetLength = len(rStitchedCigar) - offsetLength
		# Modifify the de-stitched cigar as apropriate
		if offsetType == "D":
			for j in range(0, abs(offset)):
				rOutCigar.insert(0, offsetType)
		else:
			pass
			# TODO: Handle insertions
			"""
			for j in range(0, abs(offset)):
				# Here, the insertion is too big to be explained by soft-clipping alone. Rescue failed
				if i + j >= len(rStitchedCigar):
					break
				rOutCigar[i + j] = offsetType
			"""
	# There are soft-clipped bases on the forward and reverse read. In this case,
	else:
		pass

	return fOutCigar, rOutCigar


def processRead(read, origReads, counter):
	"""
	Seperates Stitched Reads

	If a read is passed which is not stitched, it is simply returned.
	If a N-read is passed, nothing is returned
	If a stitched read was passed, the bases originating from each original read are identified from the stitched tag
	One or two reads are returned depending on if the original reads completely overlap, or if one is soft-clipped
	Also identifies reads where the stitched alignment is different than the original alignment, and either corrects the
	alignment (if the difference is likely due to a real INDEL), or discards them

	Args:
		read: A pysam.AlignmentSequence() object coresponding to a read from the merged (stitched) BAM file
		origReads: A list containing one or two pysam.AlignmentSequence(), coresponding to the original reads for the stitched read
		counter: An int listing the number of reads which have been processed so far
	Returns:
		outReads: A list containing zero, one, or two pysam.AlignmentSequence() objects, coresponding to the de-stitched (if-applicable) reads
	"""

	# If the read is simply a bunch of Ns, ignore this read
	NRead = True
	for base in read.query_sequence:
		if base != "N":
			NRead = False
			break
	if NRead:
		return []

	# Stitcher sets a tag saving the original cigar strings of each read. Lets parse that out, and see what type of read we are dealing with here
	# Easiest case: The read was NOT merged by stitcher. In this case, just save the read, uneddited.
	if not read.has_tag("XD"):
		return [read]

	sepReads = overlayCigars(read.cigarstring, read.get_tag("XD"), read.query_sequence)

	# =========================================================================================== #
	# This may just look like needless blocks of code, but there are subtle differences between each
	# Beware when making changes!
	# =========================================================================================== #

	# One read completely encapulates the other. In this case, simply output the larger read alone
	if not sepReads["sCigarList"]:
		if not sepReads["rCigarList"]:
			fRead = copy.deepcopy(read)
			fRead.flag = origReads[0].flag
			fRead.template_length = origReads[0].template_length
			fRead.next_reference_id = read.reference_id
			fRead.next_reference_start = read.reference_start
			return [fRead]
		else:
			rRead = copy.deepcopy(read)
			rRead.flag = origReads[1].flag
			rRead.template_length = origReads[1].template_length
			rRead.next_reference_id = read.reference_id
			rRead.next_reference_start = read.reference_start
			return [rRead]
	# Both reads overlap completely. In this case, just duplicate the read, and set the flags apropriately
	elif not sepReads["fCigarList"] and not sepReads["rCigarList"]:
		oRead = copy.deepcopy(read)
		oRead.flag = 67 # Reads are paired, first in pair, but no strand info
		oRead.template_length = origReads[0].template_length
		oRead.next_reference_id = read.reference_id
		oRead.next_reference_start = read.reference_start
		return [oRead]
	# There are no bases unique to the foward read, but some are unique to the reverse read
	elif not sepReads["fCigarList"] and sepReads["rCigarList"]:

		if sepReads["sOffset"] > sepReads["rOffset"]:
			stitchedOffset = read.reference_start + effectiveLengthFromCigar(cigarToList(read.cigarstring)) - origReads[0].reference_start - effectiveLengthFromCigar(cigarToList(origReads[0].cigarstring))
			if stitchedOffset != 0:
				# I'm going to throw out these reads, as they seem to be errors by stitcher more commonly than not
				return []

			rStart = getStartIndex(read.reference_start, sepReads["rCigarList"], origReads[1].cigarstring)
			fStart = getStartIndex(read.reference_start + effectiveLengthFromCigar(sepReads["rCigarList"]), sepReads["sCigarList"], origReads[0].cigarstring)
			fRead = createRead(read, origReads[0].flag, sepReads["sOffset"], sepReads["sOffset"] + len(sepReads["sSeq"]), origReads[0].template_length, listToCigar(sepReads["sCigarList"]), rStart, fStart)
			rRead = createRead(read, origReads[1].flag, 0, len(sepReads["rSeq"]), origReads[1].template_length, listToCigar(sepReads["rCigarList"]), fStart, rStart)
		else:
			stitchedOffset = read.reference_start + effectiveLengthFromCigar(cigarToList(read.cigarstring)) - origReads[1].reference_start - effectiveLengthFromCigar(cigarToList(origReads[1].cigarstring))
			if stitchedOffset != 0:
				# I'm going to throw out these reads, as they seem to be errors by stitcher more commonly than not
				return []
			rStart = getStartIndex(read.reference_start + effectiveLengthFromCigar(sepReads["sCigarList"]), sepReads["rCigarList"], origReads[1].cigarstring)
			fStart = getStartIndex(read.reference_start, sepReads["sCigarList"], origReads[0].cigarstring)
			fRead = createRead(read, origReads[0].flag, 0, len(sepReads["sSeq"]), origReads[0].template_length, listToCigar(sepReads["sCigarList"]), rStart, fStart)
			rRead = createRead(read, origReads[1].flag, sepReads["rOffset"], sepReads["rOffset"] + len(sepReads["rSeq"]), origReads[1].template_length, listToCigar(sepReads["rCigarList"]), fStart, rStart)
		return [fRead, rRead]

	# There are no bases unique to the foward read, but some are unique to the reverse read
	elif sepReads["fCigarList"] and not sepReads["rCigarList"]:

		if sepReads["sOffset"] > sepReads["fOffset"]:
			stitchedOffset = read.reference_start + effectiveLengthFromCigar(cigarToList(read.cigarstring)) - origReads[1].reference_start - effectiveLengthFromCigar(cigarToList(origReads[1].cigarstring))
			if stitchedOffset != 0:
				# I'm going to throw out these reads, as they seem to be errors by stitcher more commonly than not
				return []
			rStart = getStartIndex(read.reference_start + effectiveLengthFromCigar(sepReads["fCigarList"]), sepReads["sCigarList"], origReads[1].cigarstring)
			fStart = getStartIndex(read.reference_start, sepReads["sCigarList"], origReads[0].cigarstring)
			fRead = createRead(read, origReads[0].flag, sepReads["fOffset"], sepReads["fOffset"] + len(sepReads["fSeq"]), origReads[0].template_length, listToCigar(sepReads["fCigarList"]), rStart, fStart)
			rRead = createRead(read, origReads[1].flag, sepReads["sOffset"], sepReads["sOffset"] + len(sepReads["sSeq"]), origReads[1].template_length, listToCigar(sepReads["sCigarList"]), fStart, rStart)
		else:
			stitchedOffset = read.reference_start + effectiveLengthFromCigar(cigarToList(read.cigarstring)) - origReads[0].reference_start - effectiveLengthFromCigar(cigarToList(origReads[0].cigarstring))
			if stitchedOffset != 0:
				# I'm going to throw out these reads, as they seem to be errors by stitcher more commonly than not
				return []
			rStart = getStartIndex(read.reference_start, sepReads["sCigarList"], origReads[1].cigarstring)
			fStart = getStartIndex(read.reference_start + effectiveLengthFromCigar(sepReads["sCigarList"]), sepReads["fCigarList"], origReads[0].cigarstring)
			fRead = createRead(read, origReads[0].flag, sepReads["fOffset"], sepReads["fOffset"] + len(sepReads["fSeq"]), origReads[0].template_length, listToCigar(sepReads["fCigarList"]), rStart, fStart)
			rRead = createRead(read, origReads[1].flag, sepReads["sOffset"], sepReads["sOffset"] + len(sepReads["sSeq"]), origReads[1].template_length, listToCigar(sepReads["sCigarList"]), fStart, rStart)
		return [fRead, rRead]

	elif sepReads["fCigarList"] and sepReads["rCigarList"]:

		if sepReads["fOffset"] > sepReads["rOffset"]:

			stitchedOffset = read.reference_start + effectiveLengthFromCigar(cigarToList(read.cigarstring)) - origReads[0].reference_start - effectiveLengthFromCigar(cigarToList(origReads[0].cigarstring))
			if stitchedOffset != 0:
				sepReads["rCigarList"], sepReads["fCigarList"] = findIndel(origReads[1], origReads[0], sepReads["rCigarList"], sepReads["fCigarList"], stitchedOffset)

			# Randomly assign the shared bases to one of the two reads
			# If an indel was added to one of the two reads, add the shared bases to that read
			if sepReads["fCigarList"][0] == "D" or sepReads["fCigarList"][0] == "I":
				sepReads["fOffset"] = sepReads["sOffset"]
				sepReads["fSeq"] = sepReads["sSeq"] + sepReads["fSeq"]
				sepReads["sCigarList"].extend(sepReads["fCigarList"])
				sepReads["fCigarList"] = sepReads["sCigarList"]
			elif sepReads["rCigarList"][-1] == "D" or sepReads["rCigarList"][-1] == "I":
				sepReads["rSeq"] = sepReads["rSeq"] + sepReads["sSeq"]
				sepReads["rCigarList"].extend(sepReads["sCigarList"])
			elif counter % 2 == 0:
				sepReads["fOffset"] = sepReads["sOffset"]
				sepReads["fSeq"] = sepReads["sSeq"] + sepReads["fSeq"]
				sepReads["sCigarList"].extend(sepReads["fCigarList"])
				sepReads["fCigarList"] = sepReads["sCigarList"]
			else:
				sepReads["rSeq"] = sepReads["rSeq"] + sepReads["sSeq"]
				sepReads["rCigarList"].extend(sepReads["sCigarList"])
			rStart = getStartIndex(read.reference_start, sepReads["rCigarList"], origReads[1].cigarstring)
			fStart = getStartIndex(read.reference_start + effectiveLengthFromCigar(sepReads["rCigarList"]), sepReads["fCigarList"], origReads[0].cigarstring)
			rRead = createRead(read, origReads[1].flag, sepReads["rOffset"], sepReads["rOffset"] + len(sepReads["rSeq"]), origReads[1].template_length, listToCigar(sepReads["rCigarList"]), rStart, fStart)
			fRead = createRead(read, origReads[0].flag, sepReads["fOffset"], sepReads["fOffset"] + len(sepReads["fSeq"]), origReads[0].template_length, listToCigar(sepReads["fCigarList"]), fStart, rStart)
		elif sepReads["rOffset"] > sepReads["fOffset"]:

			stitchedOffset = read.reference_start + effectiveLengthFromCigar(cigarToList(read.cigarstring)) - origReads[1].reference_start - effectiveLengthFromCigar(cigarToList(origReads[1].cigarstring))
			if stitchedOffset != 0:
				sepReads["fCigarList"], sepReads["rCigarList"] = findIndel(origReads[0], origReads[1], sepReads["fCigarList"], sepReads["rCigarList"], stitchedOffset)

			# If an indel was added to one of the two reads, add the shared bases to that read
			if sepReads["fCigarList"][-1] == "D" or sepReads["fCigarList"][-1] == "I":
				sepReads["fSeq"] = sepReads["fSeq"] + sepReads["sSeq"]
				sepReads["fCigarList"].extend(sepReads["sCigarList"])
			elif sepReads["rCigarList"][0] == "D" or sepReads["rCigarList"][0] == "I":
				sepReads["rOffset"] = sepReads["sOffset"]
				sepReads["rSeq"] = sepReads["sSeq"] + sepReads["rSeq"]
				sepReads["sCigarList"].extend(sepReads["rCigarList"])
				sepReads["rCigarList"] = sepReads["sCigarList"]
			# Randomly assign the shared bases to one of the two reads
			elif counter % 2 == 0:
				sepReads["fSeq"] = sepReads["fSeq"] + sepReads["sSeq"]
				sepReads["fCigarList"].extend(sepReads["sCigarList"])
			else:
				sepReads["rOffset"] = sepReads["sOffset"]
				sepReads["rSeq"] = sepReads["sSeq"] + sepReads["rSeq"]
				sepReads["sCigarList"].extend(sepReads["rCigarList"])
				sepReads["rCigarList"] = sepReads["sCigarList"]

			rStart = getStartIndex(read.reference_start + effectiveLengthFromCigar(sepReads["fCigarList"]), sepReads["rCigarList"], origReads[1].cigarstring)
			fStart = getStartIndex(read.reference_start, sepReads["fCigarList"], origReads[0].cigarstring)
			rRead = createRead(read, origReads[1].flag, sepReads["rOffset"], sepReads["rOffset"] + len(sepReads["rSeq"]), origReads[1].template_length, listToCigar(sepReads["rCigarList"]), fStart, rStart)
			fRead = createRead(read, origReads[0].flag, 0, len(sepReads["fSeq"]), origReads[0].template_length, listToCigar(sepReads["fCigarList"]), rStart, fStart)
		else:
			# Something awful has happened here. I'm not even sure how this is possible, but it does occur *very* rarely in some samples
			# TODO: Figure out what is going on here. For now, I am just going to discard these reads
			return []

		return [fRead, rRead]


def main(args=None):

	# Argument parsing
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
				param = param[1:-1].split(",")

			# WARNING: Command line arguments will be SUPERSCEEDED BY CONFIG FILE ARGUMENTS
			cmdArgs[option] = param

	# Opens BAM files for reading
	mergedBAM = pysam.AlignmentFile(args.input, "rb")
	originalBAM = pysam.AlignmentFile(args.unstitched_input, "rb")

	# Check to ensure the BAM files are sorted by read name
	checkSort(mergedBAM)
	checkSort(originalBAM)

	# Checks to ensure the BAM files were given in the correct order
	checkInput(mergedBAM)

	printPrefix = "PRODUSE-SPLITMERGE"
	sys.stdout.write("\t".join([printPrefix, time.strftime('%X'), "Starting...\n"]))
	counter = 0

	with pysam.AlignmentFile(args.output, "wb", header=mergedBAM.header) as outputBAM:

		# Read from both files simultaneously
		originalReads = originalBAM.fetch(until_eof=True)
		for read in mergedBAM.fetch(until_eof=True):

			# Pull coresponding original reads from the original BAM file
			matchedOrig = [originalReads.__next__()]
			# Skip over reads missed by stitcher
			while matchedOrig[0].query_name != read.query_name:
				matchedOrig = [originalReads.__next__()]

			if read.has_tag("XD"):
				# This read has been stitched. Pull an additional read from the BAM file
				matchedOrig.append(originalReads.__next__())
				# Sanity check: If this does not have the same name, then the BAMs are likely unsorted
				if matchedOrig[1].query_name != read.query_name:
					sys.stderr.write("ERROR: The reads pulled from the pre-stitched BAM do not corespond to the reads from the post-stitched BAM\n")
					sys.stderr.write("Ensure both BAMs are sorted by READ NAME (try 'samtools sort -n <BAMfile>\n\n")
					sys.stderr.write("Stitched read name: %s\n" % (read.query_name))
					sys.stderr.write("Pre-stitched read name: %s\n" % (matchedOrig[0].query_name))
					sys.exit(1)

				# Ensure the forward read is in position one
				if matchedOrig[0].flag != 99 and matchedOrig[0].flag != 163:
					matchedOrig.reverse()

			# De-stitch the read
			outputReads = processRead(read, matchedOrig, counter)

			for outRead in outputReads:
				outputBAM.write(outRead)
			counter += 1
			if counter % 100000 == 0:
				sys.stdout.write("\t".join([printPrefix, time.strftime('%X'), "Reads Processed: " + str(counter) + "\n"]))

	sys.stdout.write("\t".join([printPrefix, time.strftime('%X'), "Reads Processed: " + str(counter) + "\n"]))
	mergedBAM.close()
	originalBAM.close()


if __name__ == "__main__":
	main()
