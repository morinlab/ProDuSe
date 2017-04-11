#! /usr/bin/env python

import configargparse
import configparser
import sys
import os
import pysam
import copy
import re
import time

parser = configargparse.ArgumentParser(description="Splits reads merged ")
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

	Args:
		array: A list coresponding to the cigar of interest
	Returns:
		cigar: A string representing the cigar sequence
	"""
	if len(array) == 0:
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


def findIndel(fRead, rRead, predEnd, actualEnd, iForwardCigar, iReverseCigar):
	"""
	Find the location of insertions or deletions in the paired reads

	If soft-clipped bases are unique to either the forward or reverse read, assume the indel is in that read
	If there are soft-clipped bases on both read, perform realignment of the soft-clipped sections to find the
	locations of the indels

	TODO: This may fail if the cigar ends with a insertion or deletion; however, considering there is an large,
	earlier indel, this is extremely unlikely

	Args:
		fRead: A pysam.AlignmentSegment() object coresponding to the foward read
		rRead: A pysam.AlignmentSegment() object coresponding to the reverse read
		predEnd: The predicted end locus of the next read, based upon stitcher's alignment
		actualEnd: The actual end locus of the next read, based upon the unmerged reads
		iForwardCigar: The forward cigar string obtained from splitting the stitched read
		iReverseCigar: The reverse cigar string obtained from splitting the stitched read
	Returns:
		oForwardCigar, oReverseCigar: Cigar strings modified with the indel
	"""

	# Determine the type of event (Insertion or deletion)
	if predEnd < actualEnd:
		event = "D"
	else:
		event = "I"
	estEventSize = abs(predEnd - actualEnd)

	iFList = cigarToList(iForwardCigar)
	iRList = cigarToList(iReverseCigar)
	fCigar = cigarToList(fRead.cigarstring)
	rCigar = cigarToList(rRead.cigarstring)
	# If soft-clipped bases are unique to the forward read, add the INDEL there
	if fCigar[-1] == "S" and rCigar[0] != "S":
		# Find the location of the indel
		i = len(fCigar) - 1
		while fCigar[i] == "S" and i > 0:
			i -= 1
		# i += 1
		# Here, lets determine the amount of soft-clipped bases on the forward read, and assume those are matched to the overlapping
		eventSize = estEventSize
		# Add the event
		if event == "D":
			for j in range(0, eventSize):
				iFList.insert(len(iFList) - 1, event)
		else:
			for j in range(0, eventSize):
				iFList[len(iFList) - 1 - j] = event

	# If the soft-clipped bases are unique to the reverse read, add the indel there
	elif fCigar[-1] != "S" and rCigar[0] == "S":
		# Find the location of the indel
		i = 0
		while rCigar[i] == "S" and i < len(rCigar):
			i += 1

		# Same case here: lets determine the amount of soft-clipped bases on the forward read, and assume those are matched to the overlapping
		eventSize = estEventSize
		# Add the event
		if event == "D":
			for j in range(0, eventSize):
				iRList.insert(0, event)
		else:
			for j in range(0, eventSize):
				iRList[j] = event
	# Hardest case: There are two indels, one unique to the forward and reverse read
	# TODO: Finish this, since this is unlikely
	# This is going to require some local realignment to determine the exact alignment of the overlapping reads
	else:
		return None, None
	oForwardCigar = listToCigar(iFList)
	oReverseCigar = listToCigar(iRList)
	return oForwardCigar, oReverseCigar


def mergeCigars(*cigars):
	"""
	Combines multiple cigar strings, in the order they are recieved

	Args:
		cigars: A list containing cigar strings
	Returns:
		outCigar: The merged cigar string
	"""
	outCigarList = []
	for cigar in cigars:
		cigarList = cigarToList(cigar)
		outCigarList.extend(cigarList)
	outCigar = listToCigar(outCigarList)

	return outCigar


def overlayCigars(stitchedCigar, origCigar):
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

	# Find the sections coresponding to the forward and reverse reads
	softClippedRegex = '[0-9]+S'
	allClippedRegex = '[0-9]+[A-Z]'
	softClippedCigars = re.findall(softClippedRegex, origCigar)
	allCigars = re.findall(allClippedRegex, origCigar)
	forwardCigar = None
	forwardOutCigar = None
	forwardStart = 0
	reverseCigar = None
	reverseOutCigar = None
	reverseStart = 0
	sharedCigar = None
	sharedStart = 0
	stitchedList = cigarToList(stitchedCigar)
	origList = cigarToList(origCigar)
	if len(stitchedList) != len(origList):
		raise TypeError("The two cigar strings cannot be divided, as they are of different lengths")
	elif "S" not in origCigar:
		raise TypeError("Cigar %s cannot be used to divide %s, as no soft-clipped bases are indicated" % (origCigar, stitchedCigar))

	solutionFound = False
	for cigar in softClippedCigars:
		for i in range(0, len(allCigars)):
			if allCigars[i] == cigar:
				if i + 1 > len(allCigars):
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
						forwardCigar = stitchedCigar
						return forwardCigar, 0, sharedCigar, 0, reverseCigar, 0
					# The forward read is completely encliped by the reverse read
					elif "F" not in cigar2 and "R" in cigar2 and "F" not in cigar1 and "R" in cigar1:
						reverseCigar = stitchedCigar
						return forwardCigar, 0, sharedCigar, 0, reverseCigar, 0
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

		# Find the largest soft-clipped region, and assume that coresponds to the overlap
		maxSize = 0
		maxCigar = None
		for cigar in softClippedCigars:
			cigarLen = len(cigar[:-1])
			if cigarLen > maxSize:
				maxSize = cigarLen
				maxCigar = cigar
		# Split into unique forward and reverse sections
		cigar1, cigar2 = origCigar.split(maxCigar)
		cigar = maxCigar
		# Next, calculate the number of forward and reverse bases coresponding to each read for each cigar
		cigar1F = 1
		cigar1R = 1
		cigar1List = cigarToList(cigar1)
		for char in cigar1List:
			if char == "F":
				cigar1F += 1
			elif char == "R":
				cigar1R += 1

		cigar2F = 1
		cigar2R = 1
		cigar2List = cigarToList(cigar2)
		for char in cigar2List:
			if char == "F":
				cigar2F += 1
			elif char == "R":
				cigar2R += 1

		if float(cigar1F) / float(cigar1R) > float(cigar2F) / float(cigar2R):
			forwardCigar = cigar1
			reverseCigar = cigar2
		else:
			forwardCigar = cigar2
			reverseCigar = cigar1

	# Perform the overlapping
	sharedOutList = stitchedList
	if forwardCigar:
		forwardStart = len(cigarToList(origCigar[:origCigar.find(forwardCigar)]))
		forwardEnd = len(cigarToList(origCigar[:origCigar.find(forwardCigar) + len(forwardCigar)]))
		forwardOutList = stitchedList[forwardStart:forwardEnd]
		forwardOutCigar = listToCigar(forwardOutList)
		for i in range(forwardStart, forwardEnd):
			sharedOutList[i] = "N"
	if reverseCigar:
		reverseStart = len(cigarToList(origCigar[:origCigar.find(reverseCigar)]))
		reverseEnd = len(cigarToList(origCigar[:origCigar.find(reverseCigar) + len(reverseCigar)]))
		reverseOutList = stitchedList[reverseStart:reverseEnd]
		reverseOutCigar = listToCigar(reverseOutList)
		for i in range(reverseStart, reverseEnd):
			sharedOutList[i] = "N"
	tmpList = []
	sharedStart = 0
	leadingNs = True
	for char in sharedOutList:
		if char != "N":
			leadingNs = False
			tmpList.append(char)
		elif leadingNs:
			sharedStart += 1
	sharedOutList = tmpList
	sharedOutCigar = listToCigar(sharedOutList)

	# Edge case: If the bases unique to a read are entirely soft-clipped, append the sequence to another, as it is not informative alone
	if forwardOutCigar is not None and "M" not in forwardOutCigar and "I" not in forwardOutCigar and "D" not in forwardOutCigar:
		if sharedOutCigar is None:
			# If the forward cigar was second, append it to the reverse cigar
			if forwardStart > reverseStart:
				reverseOutCigar = mergeCigars(reverseOutCigar, forwardOutCigar)
			else:
				reverseOutCigar = mergeCigars(forwardOutCigar, reverseOutCigar)
				reverseStart = forwardStart
		else:
			# If the forward cigar is after the shared cigar, append it to the merged cigar. Otherwise, prepend the shared cigar
			if forwardStart > sharedStart:
				sharedOutCigar = mergeCigars(sharedOutCigar, forwardOutCigar)
			else:
				sharedOutCigar = mergeCigars(forwardOutCigar, sharedOutCigar)
				sharedStart = forwardStart
		forwardOutCigar = None
	if reverseOutCigar is not None and "M" not in reverseOutCigar and "I" not in reverseOutCigar and "D" not in reverseOutCigar:
		if sharedOutCigar is None:
			# If the reverse cigar was second, append it to the forward cigar
			if reverseStart > forwardStart:
				forwardOutCigar = mergeCigars(forwardOutCigar, reverseOutCigar)
			else:
				forwardOutCigar = mergeCigars(reverseOutCigar, forwardOutCigar)
				forwardStart = reverseStart
		else:
			# If the reverse cigar is after the shared cigar, append it to the merged cigar. Otherwise, prepend the shared cigar
			if reverseStart > sharedStart:
				sharedOutCigar = mergeCigars(sharedOutCigar, reverseOutCigar)
			else:
				sharedOutCigar = mergeCigars(reverseOutCigar, sharedOutCigar)
				sharedStart = reverseStart
		reverseOutCigar = None
	if sharedOutCigar is not None and "M" not in sharedOutCigar and "I" not in sharedOutCigar and "D" not in sharedOutCigar:
		if forwardOutCigar is not None:
			if sharedStart > forwardStart:
				forwardOutCigar = mergeCigars(forwardOutCigar, sharedOutCigar)
			else:
				forwardOutCigar = mergeCigars(sharedOutCigar, forwardOutCigar)
				forwardStart = sharedStart
		else:
			if sharedStart > reverseStart:
				reverseOutCigar = mergeCigars(reverseOutCigar, sharedOutCigar)
			else:
				reverseOutCigar = mergeCigars(sharedOutCigar, reverseOutCigar)
				reverseStart = sharedStart
		sharedOutCigar = None
	return forwardOutCigar, forwardStart, sharedOutCigar, sharedStart, reverseOutCigar, reverseStart


def getStartIndex(start, cigar, origCigar=None, previousCigar=None):
	"""
	Calculates the start location of the next read, based upon leading soft-clipped bases and indels in the previous read

	Args:
		start: Expected start locus
		cigar: Cigar sequence of the read in question
		previousCigar: Cigar string of the previous read
	Returns:
		start: Corrected start locus
	"""
	# Check the cigar sequence for leading soft-clipped bases, and offset as necessary
	cigarList = cigarToList(cigar)
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

	# Correct for insertions and deletions in the previous read
	if previousCigar:
		pCigarList = cigarToList(previousCigar)
		i = 0
		for char in pCigarList:
			if i == 0 and char == "S":
				start -= 1
			elif char == "D":
				i += 1
				start += 1
			elif char == "I":
				i += 1
				start -= 1
	return start


def getEndIndex(start, cigar):
	"""
	Determines the end locus of the mapped read

	Args:
		start: Start locus
		cigar: Cigar sequence of the read
	Returns:
		end: End locus
	"""
	cigarList = cigarToList(cigar)
	end = start
	start = True
	for char in cigarList:
		if start and char == "S":
			continue
		elif char != "I":
			end += 1
		start = False
	return end


def effectiveLengthFromCigar(cigar):
	"""
	Determines the actual length of the cigar sequence relative to the reference

	Args:
		cigar: Cigar sequence
	Returns:
		length: Length of the sequence
	"""
	length = 0
	cigarList = cigarToList(cigar)
	for char in cigarList:
		if char != "D":
			length += 1
	return length


def createRead(iRead, flag, divStart, divEnd, tempLen, cigarSeq, mateStart, start):
	"""
	TODO
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


def stealSomeBases(forwardCigar, reverseCigar, basesBorrowed=1):
	"""
	Caps an INDEL with bases from the other read

	Idenfied the read which is ends or begins with an indel. This indel is capped
	by "borrowing" a handful bases from the other read

	TODO: Add support for "split" events.

	Args:
		forwardCigar: Earlier Cigar String (starts first)
		reverseCigar: Later Cigar string (starts second)
		basesBorrowed: Max number of bases to "borrow"
	Returns:
		oForwardCigar: Forward Cigar String
		oReverseCigar: Reverse Cigar String
	"""
	forwardCigarList = cigarToList(forwardCigar)
	reverseCigarList = cigarToList(reverseCigar)
	if forwardCigar.endswith("D") or forwardCigar.endswith("I"):

		# If there are not enough bases in the reverse read, don't move anything
		if basesBorrowed >= len(reverseCigarList):
			return forwardCigar, reverseCigar
		for i in range(0, basesBorrowed):
			forwardCigarList.append(reverseCigarList[0])
			del reverseCigarList[0]
	else:
		# If there are not enough bases in the forward read, don't move anything
		if basesBorrowed >= len(forwardCigarList):
			return forwardCigar, reverseCigar
		for i in range(0, basesBorrowed):
			reverseCigarList.insert(0, forwardCigarList[-1])
			del reverseCigarList[-1]

	return listToCigar(forwardCigarList), listToCigar(reverseCigarList)


def processRead(read, origReads, n):
	"""
	Splits merged reads into a forward and reverse read

	If bases are shared between the two reads, they are divided as follows:
		1) If all bases in the read overlap, a single read is returned
		2) If some bases are unique to a single read, the common bases are assigned to the opposite read
		3) If the forward and reverse read each contain unique bases, "randomly" assign the shared bases to one of the reads

	Args:
		read: A pysam.AlignedSegment() object coresponding to the read of interest
		origReads: A list containing pysam.AlignedSegment() object(s) coresponding to the original read
	Returns:
		reads: An array containing pysam.AlignedSegment() objects

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
		n += 1
		return [read]
	# Ok, this read was merged by stitcher. We need to parse the elements of the cigar string, and determine which bases are unique to each read
	else:
		n += 1

		forwardCigar, forwardStart, sharedCigar, sharedStart, reverseCigar, reverseStart = overlayCigars(read.cigarstring, read.get_tag("XD"))

		# Determine if the read originates from the foward or reverse strand
		if ":-:" in read.query_name:
			forwardFlag = 163
			reverseFlag = 83
		else:
			forwardFlag = 99
			reverseFlag = 147

		# All the bases originate from a single read, and the second read was simply a subset of the first
		# In this case, simply return that read
		if not sharedCigar:
			if not reverseCigar:
				fRead = copy.deepcopy(read)
				fRead.flag = forwardFlag
				fRead.template_length = len(read.query_sequence)
				fRead.next_reference_id = read.reference_id
				fRead.next_reference_start = read.reference_start
				return [fRead]
			else:
				rRead = copy.deepcopy(read)
				rRead.flag = reverseFlag
				rRead.template_length = len(read.query_sequence)
				rRead.next_reference_id = read.reference_id
				rRead.next_reference_start = read.reference_start
				return [rRead]

		# There are no bases unique either read. They overalap completely
		# In this case, just duplicate the read, set the apropriate flags, and save to the output BAM file
		elif not forwardCigar and not reverseCigar:

			fRead = copy.deepcopy(read)
			fRead.flag = forwardFlag
			fRead.template_length = len(read.query_sequence)
			fRead.next_reference_id = read.reference_id
			fRead.next_reference_start = read.reference_start
			return [fRead]

		# Bases are unique to the reverse read
		# In this case, simpily truncate the total collapsed read to the necessary length, and use that as the forward read
		elif not forwardCigar and reverseCigar:

			# If the pedicted and actual length of the sequence unique to the forward strand do not match up, find the location of the ins/deletion
			predEnd = getEndIndex(read.reference_start, read.cigarstring)
			actualEnd = max(getEndIndex(origReads[0].reference_start, origReads[0].cigarstring), getEndIndex(origReads[1].reference_start, origReads[1].cigarstring))
			if reverseStart > sharedStart:

				if predEnd != actualEnd:
					sharedCigar, reverseCigar = findIndel(origReads[1], origReads[0], predEnd, actualEnd, sharedCigar, reverseCigar)
					# Temp measure: Remove reads which were not adjusted. Fix this later
					if sharedCigar is None:
						return []
					sharedCigar, reverseCigar = stealSomeBases(sharedCigar, reverseCigar)

				fStart = getStartIndex(read.reference_start + effectiveLengthFromCigar(sharedCigar), reverseCigar, origReads[1].cigarstring, sharedCigar)
				rStart = getStartIndex(read.reference_start + sharedStart, sharedCigar, origReads[1].cigarstring)
				fRead = createRead(read, forwardFlag, reverseStart, reverseStart + effectiveLengthFromCigar(reverseCigar), origReads[0].template_length, reverseCigar, read.reference_start, fStart)
				rRead = createRead(read, reverseFlag, sharedStart, reverseStart, origReads[1].template_length, sharedCigar, fRead.reference_start, rStart)
				fRead.next_reference_start = rRead.reference_start
			else:
				if predEnd != actualEnd:
					reverseCigar, sharedCigar = findIndel(origReads[0], origReads[1], predEnd, actualEnd, reverseCigar, sharedCigar)
					# Temp measure: Remove reads which were not adjusted. Fix this later
					if sharedCigar is None:
						return []
					reverseCigar, sharedCigar = stealSomeBases(reverseCigar, sharedCigar)

				fStart = getStartIndex(read.reference_start + effectiveLengthFromCigar(reverseCigar), sharedCigar, origReads[0].cigarstring, reverseCigar)
				rStart = getStartIndex(read.reference_start + reverseStart, reverseCigar, origReads[1].cigarstring)
				fRead = createRead(read, forwardFlag, sharedStart, sharedStart + effectiveLengthFromCigar(sharedCigar), origReads[0].template_length, sharedCigar, read.reference_start, fStart)
				rRead = createRead(read, reverseFlag, reverseStart, sharedStart, origReads[1].template_length, reverseCigar, fRead.reference_start, rStart)
				fRead.next_reference_start = rRead.reference_start

			return [fRead, rRead]

		# Some bases are unique to the foward read
		# In this case, truncate the read to apropriate length and position, and use that as the reverse read
		elif forwardCigar and not reverseCigar:

			predEnd = getEndIndex(read.reference_start, read.cigarstring)
			actualEnd = max(getEndIndex(origReads[0].reference_start, origReads[0].cigarstring), getEndIndex(origReads[1].reference_start, origReads[1].cigarstring))
			if forwardStart > sharedStart:
				# If the pedicted and actual length of the sequence unique to the forward strand do not match up, find the location of the ins/deletion
				if predEnd != actualEnd:
					forwardCigar, sharedCigar = findIndel(origReads[1], origReads[0], predEnd, actualEnd, sharedCigar, forwardCigar)
					# Temp measure: Remove reads which were not adjusted. Fix this later
					if sharedCigar is None:
						return []
					forwardCigar, sharedCigar = stealSomeBases(sharedCigar, forwardCigar)

				rStart = getStartIndex(read.reference_start + effectiveLengthFromCigar(sharedCigar), forwardCigar, origReads[0].cigarstring, sharedCigar)
				fStart = getStartIndex(read.reference_start + sharedStart, sharedCigar, origReads[1].cigarstring)
				rRead = createRead(read, forwardFlag, forwardStart, forwardStart + effectiveLengthFromCigar(forwardCigar), origReads[0].template_length, forwardCigar, read.reference_start, rStart)
				fRead = createRead(read, reverseFlag, sharedStart, forwardStart, origReads[1].template_length, sharedCigar, rRead.reference_start, fStart)
				rRead.next_reference_start = fRead.reference_start
			else:
				# If the pedicted and actual length of the sequence unique to the forward strand do not match up, find the location of the ins/deletion
				if predEnd != actualEnd:
					forwardCigar, sharedCigar = findIndel(origReads[0], origReads[1], predEnd, actualEnd, forwardCigar, sharedCigar)
					# Temp measure: Remove reads which were not adjusted. Fix this later
					if sharedCigar is None:
						return []
					forwardCigar, sharedCigar = stealSomeBases(forwardCigar, sharedCigar)

				rStart = getStartIndex(read.reference_start + effectiveLengthFromCigar(forwardCigar), sharedCigar, origReads[1].cigarstring, forwardCigar)
				fStart = getStartIndex(read.reference_start + forwardStart, forwardCigar, origReads[0].cigarstring)
				rRead = createRead(read, reverseFlag, sharedStart, sharedStart + effectiveLengthFromCigar(sharedCigar), origReads[1].template_length, sharedCigar, read.reference_start, rStart)
				fRead = createRead(read, forwardFlag, forwardStart, sharedStart, origReads[0].template_length, forwardCigar, rRead.reference_start, fStart)
				rRead.next_reference_start = fRead.reference_start

			return [fRead, rRead]

		# Bases are unique to the forward and reverse read
		elif forwardCigar and reverseCigar:

			cigarRegex = '[0-9]+[A-Z]'
			reverseCigarStart = " "
			forwardCigarStart = " "

			predEnd = getEndIndex(read.reference_start, read.cigarstring)
			actualEnd = max(getEndIndex(origReads[0].reference_start, origReads[0].cigarstring), getEndIndex(origReads[1].reference_start, origReads[1].cigarstring))

			if forwardStart > reverseStart:

				# If the pedicted and actual length of the sequence unique to the forward strand do not match up, find the location of the ins/deletion
				if predEnd != actualEnd:
					forwardCigar, reverseCigar = findIndel(origReads[1], origReads[0], predEnd, actualEnd, reverseCigar, forwardCigar)
					# Temp measure: Remove reads which were not adjusted. Fix this later
					if forwardCigar is None:
						return []
					forwardCigarStart = re.findall(cigarRegex, forwardCigar)[0]

				# TODO: Reshuffle code to be more efficient
				# Assign the shared nucleotides to either the forward or reverse read.
				# If one of the reads ends with an insertion or deletion, add the shared sequence to that read
				if "D" == forwardCigarStart[-1] or "I" == forwardCigarStart[-1]:
					rCigar = reverseCigar
					fCigar = mergeCigars(sharedCigar, forwardCigar)
				elif reverseCigar.endswith("D") or reverseCigarStart.endswith("I"):
					# Give the shared bases to the reverse read
					fCigar = forwardCigar
					rCigar = mergeCigars(reverseCigar, sharedCigar)
				# Else, Randomly assign the shared nucleotides to either the forward read or the reverse read
				elif n % 2 == 0:
					# Give the shared bases to the reverse read
					fCigar = forwardCigar
					rCigar = mergeCigars(reverseCigar, sharedCigar)
					division = forwardStart
				else:
					# Give the shared bases to the foward read
					division = sharedStart
					rCigar = reverseCigar
					fCigar = mergeCigars(sharedCigar, forwardCigar)

				rLength = effectiveLengthFromCigar(rCigar)
				division = rLength

				fStart = getStartIndex(read.reference_start + rLength, fCigar, origReads[1].cigarstring, rCigar)
				rStart = getStartIndex(read.reference_start + reverseStart, rCigar, origReads[0].cigarstring)
				rRead = createRead(read, reverseFlag, division, division + effectiveLengthFromCigar(fCigar), origReads[0].template_length, fCigar, read.reference_start, fStart)
				fRead = createRead(read, forwardFlag, reverseStart, division, origReads[1].template_length, rCigar, rRead.reference_start, rStart)
				rRead.next_reference_start = fRead.reference_start
			elif forwardStart < reverseStart:

				# If the pedicted and actual length of the sequence unique to the forward strand do not match up, find the location of the ins/deletion
				if predEnd != actualEnd:
					forwardCigar, reverseCigar = findIndel(origReads[0], origReads[1], predEnd, actualEnd, forwardCigar, reverseCigar)
					# Temp measure: Remove reads which were not adjusted. Fix this later
					if forwardCigar is None:
						return []
					reverseCigarStart = re.findall(cigarRegex, reverseCigar)[0]

				# TODO: Reshuffle code to be more efficient
				# Assign the shared nucleotides to either the forward or reverse read.
				# If one of the reads ends with an insertion or deletion, add the shared sequence to that read
				if forwardCigar.endswith("D") or forwardCigar.endswith("I"):
					rCigar = reverseCigar
					fCigarList = cigarToList(forwardCigar)
					fCigarList.extend(cigarToList(sharedCigar))
					fCigar = listToCigar(fCigarList)
				elif "D" == reverseCigarStart[-1] or "I" == reverseCigarStart[-1]:
					# Give the shared bases to the reverse read
					fCigar = forwardCigar
					rCigarList = cigarToList(sharedCigar)
					rCigarList.extend(cigarToList(reverseCigar))
					rCigar = listToCigar(rCigarList)
				# Else, Randomly assign the shared nucleotides to either the forward read or the reverse read
				elif n % 2 == 0:
					# Give the shared bases to the forward read
					rCigar = reverseCigar
					fCigarList = cigarToList(forwardCigar)
					fCigarList.extend(cigarToList(sharedCigar))
					fCigar = listToCigar(fCigarList)
				else:
					# Give the shared bases to the reverse read
					fCigar = forwardCigar
					rCigarList = cigarToList(sharedCigar)
					rCigarList.extend(cigarToList(reverseCigar))
					rCigar = listToCigar(rCigarList)

				fLength = effectiveLengthFromCigar(fCigar)
				division = fLength
				rStart = getStartIndex(read.reference_start + fLength, rCigar, origReads[1].cigarstring, fCigar)
				fStart = getStartIndex(read.reference_start + forwardStart, origReads[0].cigarstring, fCigar)
				rRead = createRead(read, reverseFlag, division, division + effectiveLengthFromCigar(rCigar), origReads[1].template_length, rCigar, read.reference_start, rStart)
				fRead = createRead(read, forwardFlag, forwardStart, division, origReads[0].template_length, fCigar, rRead.reference_start, fStart)
				rRead.next_reference_start = fRead.reference_start
			return [fRead, rRead]
	return []


def main(args=None):

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

	mergedBAM = pysam.AlignmentFile(args.input, "rb")
	originalBAM = pysam.AlignmentFile(args.unstitched_input, "rb")

	# Check to ensure the BAM files are sorted by read name
	checkSort(mergedBAM)
	checkSort(originalBAM)

	# Checks to ensure the BAM files were given in the correct order
	checkInput(mergedBAM)

	printPrefix = "PRODUSE-SPLITMERGE"
	counter = 0

	with pysam.AlignmentFile(args.output, "wb", header=mergedBAM.header) as outputBAM:

		sys.stdout.write("\t".join([printPrefix, time.strftime('%X'), "Starting...\n"]))

		originalReads = originalBAM.fetch(until_eof=True)
		for read in mergedBAM.fetch(until_eof=True):

			# Pull coresponding original reads from the original BAM file
			matchedOrig = [originalReads.__next__()]
			while matchedOrig[0].query_name != read.query_name:
				matchedOrig = [originalReads.__next__()]

			if read.has_tag("XD"):
				# This read has been stitched. Pull the next two reads from the BAM file
				matchedOrig.append(originalReads.__next__())
				# Sanity check: If this does not have the same name, then the BAMs are likely unsorted
				if matchedOrig[1].query_name != read.query_name:
					sys.stderr.write("ERROR: The reads pulled from the pre-stitched BAM do not corespond to the reads from the post-stitched BAM\n")
					sys.stderr.write("Ensure both BAMs are sorted by READ NAME (try 'samtools sort -n <BAMfile>\n\n")
					sys.stderr.write("Stitched read name: %s\n" % (read.query_name))
					sys.stderr.write("Pre-stitched read name: %s\n" % (matchedOrig[0].query_name))
					sys.exit(1)
				# Next, sort the reads so the read with the latter start is second
				if matchedOrig[0].reference_start > matchedOrig[1].reference_start:
					matchedOrig.reverse()

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
