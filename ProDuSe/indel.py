# USAGE
# Designed for use with the ProDuSe variant caller.
#
# DESCRIPTION
# This script contains the classes and methods required to relaign INDELs
# The implementation depends upon which subclass is called.
#
# For familyIndel, a mapping and haplotype-based aproach is used
#
# AUTHORs
# Christopher Rushton (ckrushto@sfu.ca)

from skbio.alignment import StripedSmithWaterman
import pysam
import copy

class familyIndel:
	"""
	Identify and realign INDELs within the same read family
	i.e. Assume that ALL sequenced originated from the same parent sequence,
	and any deviation between sequences is due to sequencing/PCR errors

	Note that using this class if the above function is NOT true will prove catastrophic
	"""

	def __init__(self, refSeq=None, reads=None, offset=0):
		"""
		Args:
			refSeq: A string containing the reference sequence, including at least 150 flanking bases
			reads: One or more pysam.AlignedSegment() objects, listing the reads to be realigned
			offset: The offset The start positions of the reads relative to the reference sequence provided
		"""

		self.reference = refSeq
		self.reads = reads
		self.refOffset = offset

	def realign(self):
		"""
		Realign the reads, and identify any INDELs (if present)

		Further description of the implemtation is ongoing
		"""

		if self.reads is None:
			raise TypeError("No sequences were provided")

		# What if all the reads in the family do not have the same sequence?
		# In this case, we will have to generate some form of consensus
		# Identify the longest read in the set. This will be used as a "seed" read, and all other reads
		# will be mapped to it
		seedSeq = None
		seedQual = None
		seedLength = 0
		remainingReads = []
		for read in self.reads:
			if read.query_length > seedLength:
				seedLength = read.query_length
				seedQual = read.query_qualities
				seedSeq = read.query_sequence
			else:
				remainingReads.append(read)

		# Create a De Brujn graph (sort of)
		graph = DeBrujnGraph(seedSeq, seedQual, seedSeq)
		# Create an alignment object for this seed read.
		seedAlign = StripedSmithWaterman(query_sequence=seedSeq, gap_open_penalty=3, gap_extend_penalty=1, score_size=2, mask_length=15, mask_auto=False)
		# Map the remaining reads to this read
		for read in remainingReads:
			# Easiest case: The reads are identical
			if read.query_sequence == seedSeq:
				graph.map(seedSeq, read.query_sequence, read.query_qualities)
			# Otherwise, map this read onto the seed read
			else:
				alignment = seedAlign(read.query_sequence)

				# Was either portion of the alignment clipped? If so, we need to restore the clipped bases
				# Was the seed clipped at the start?
				queryStartClip = alignment.query_begin
				queryEndClip = len(seedSeq) - 1 - alignment.query_end - queryStartClip
				targetStartClip = alignment.target_begin
				targetEndClip = len(read.query_sequence) - 1 - alignment.target_end_optimal - targetStartClip

				# Sanity Check: The query and target reads should not both be start or end clipped
				if (queryStartClip != 0 and targetStartClip != 0) or (queryEndClip != 0 and targetEndClip != 0):
					raise TypeError("Both the target and query reads were clipped at at the same position during SW alignment")

				# Add the trimmed edges back onto the alignments
				totalQueryAlgin = alignment.aligned_query_sequence
				totalTargetAlign = alignment.aligned_target_sequence
				if queryStartClip != 0:
					queryClip = seedSeq[0:queryStartClip]
					targetClip = "-" * queryStartClip
					totalQueryAlgin = queryClip + totalQueryAlgin
					totalTargetAlign = targetClip + totalTargetAlign

				elif targetStartClip != 0:
					targetClip = read.query_sequence[0:targetStartClip]
					queryClip = "-" * targetStartClip
					totalQueryAlgin = queryClip + totalQueryAlgin
					totalTargetAlign = targetClip + totalTargetAlign
				if queryEndClip != 0:
					queryClip = seedSeq[-1 * queryEndClip:]
					targetClip = "-" * queryEndClip
					totalQueryAlgin = totalQueryAlgin + queryClip
					totalTargetAlign = totalTargetAlign + targetClip
				elif targetEndClip != 0:
					targetClip = read.query_sequence[-1 * targetEndClip:]
					queryClip = "-" * targetEndClip
					totalQueryAlgin = totalQueryAlgin + queryClip
					totalTargetAlign = totalTargetAlign + targetClip

				graph.map(totalQueryAlgin, totalTargetAlign, read.query_qualities)

		# Next, map this consensus read to the reference
		familyConsensus, familyQual = graph.consensus()
		refAlign = StripedSmithWaterman(query_sequence=self.reference)
		alignment = refAlign(familyConsensus)

		# Has the start position of the read changed relative to the reference?
		strippedQuerySeq = alignment.aligned_query_sequence.replace("-", "")
		changedStartIndex = self.reference.find(strippedQuerySeq) - self.refOffset
		self.toAlignedSequence(alignment.aligned_target_sequence, alignment.aligned_query_sequence, familyQual, changedStartIndex)

	def toAlignedSequence(self, sequence, reference, quality, offset):
		"""
		Converts the consensus sequence into a Pysam.AlignedSegment() object

		Args:
			sequence: A String representing the consensus sequence to be converted
			reference: A String representing the reference base for each base in the sequence
			qualities: A list containing integers representing the quality of each base in sequence
			offset: An int representing the change in start index
		"""

		cigarStr = ""
		pysamSeq = ""
		pysamQual = []
		bufferLength = 0
		leadingSoftClipping = 0
		leadingBases = True  # For soft clipping
		i = 0
		for base in sequence:
			refBase = reference[i]
			qual = quality[i]

			# Matches
			if refBase == base:
				# Dump the buffer (if one exists)
				if bufferLength != 0:
					cigarStr += "M" * bufferLength
					bufferLength = 0
				cigarStr += "M"
				pysamSeq += base
				pysamQual.append(qual)
				leadingBases = False
			# Insertions
			elif refBase == "-":
				if bufferLength != 0:
					cigarStr += "M" * bufferLength
					bufferLength = 0
				cigarStr += "I"
				pysamSeq += base
				pysamQual.append(qual)
				leadingBases = False
			# Deletions
			elif base == "-":
				if bufferLength != 0:
					cigarStr += "M" * bufferLength
					bufferLength = 0
				cigarStr += "D"
				leadingBases = False
			# Mismatches
			else:
				# Do the mismatches lead the sequence? In which case, soft clip them
				if leadingBases:
					cigarStr += "S"
					leadingSoftClipping += 1
				else:
					# These could either be mismatches (M), or soft-clipped from the end of the read
					bufferLength += 1
				pysamQual.append(qual)
				pysamSeq += base
			i += 1

		if bufferLength != 0:
			cigarStr += "S" * bufferLength
			bufferLength = 0

		# Convert the string into a real cigar sequence
		cigar = ""
		previousAction = cigarStr[0]
		actionLength = 1
		for action in cigarStr[1:]:
			if action == previousAction:
				actionLength += 1
			else:
				cigar += str(actionLength) + previousAction
				previousAction = action
				actionLength = 1

		cigar += str(actionLength) + previousAction

		# Time to make the pysam object
		read = copy.deepcopy(self.reads[0])
		read.query_sequence = pysamSeq
		read.query_qualities = pysamQual
		read.cigarstring = cigar
		read.reference_start = self.reads[0].reference_start + offset + leadingSoftClipping
		print(read)


class DeBrujnGraph:
	"""
	A very basic (and probably awful) version of a debrujn graph
	"""

	def __init__(self, sequence, qualities, refSeq):
		"""
		Args:
			sequence: A string representing a nucleotide sequence
			qualities: A list of integers, representing base qualities for each base in the sequence
			refSeq: A string representing the reference genome sequence for the positions represented by "sequence"
		"""
		if len(sequence) != len(refSeq):
			raise TypeError("The reference sequence must be the same length as the template sequence")

		self.graph = []
		self.originalSeed = sequence
		for i in range(0, len(sequence)):
			currentBase = sequence[i]
			currentQual = qualities[i]
			refBase = refSeq[i]
			if currentBase == "A":
				self.graph.append({"A": [1, currentQual], "C": [0, 0], "G": [0, 0], "T": [0, 0], "-": [0, 0], "Ref": refBase})  # [Count, MaxQuality]
			elif currentBase == "C":
				self.graph.append({"A": [0, 0], "C": [1, currentQual], "G": [0, 0], "T": [0, 0], "-": [0, 0], "Ref": refBase})
			elif currentBase == "G":
				self.graph.append({"A": [0, 0], "C": [0, 0], "G": [1, currentQual], "T": [0, 0], "-": [0, 0], "Ref": refBase})
			elif currentBase == "T":
				self.graph.append({"A": [0, 0], "C": [0, 0], "G": [0, 0], "T": [1, currentQual], "-": [0, 0], "Ref": refBase})
			else:
				raise TypeError("Only nucleotide sequences [A,C,G,T] are supported, but %s was provided" % currentBase)
		self.refSkips = []  # Represents nodes which are not present in the original sequence (indicates an insertion)

	def map(self, seedSeq, newSeq, newQual):
		"""
		Add the specified sequence "newSeq" to the graph
		Args:
			seedSeq: A string representing the seed sequence following alignment with newSeq
			newSeq: A string representing the new sequence following alignment with seedSeq
			newQual: A list containing integers corespond to the base quality for each base in newSeq
		"""

		# Have gaps been introduced into the seedSeq? (i.e. is there an insertion in newSeq)
		if len(seedSeq) != len(self.originalSeed):
			# Find the location of the gaps (i.e. insertions in the new sequence)

			gaps = list(i for i in range(0, len(seedSeq)) if seedSeq[i] == "-")
			# Are these gaps already represented in the graph? If so, no problem!
			# Otherwise, they need to be added to the graph
			newGaps = list(x for x in gaps if x not in self.refSkips)
			if len(newGaps) != 0:
				for gap in newGaps:
					self.graph.insert(gap, {"A": 0, "C": 0, "G": 0, "T": 0, "-": 1, "Ref": "-"})
					for i in range(0, len(self.refSkips)):
						if self.refSkips[i] < gap:
							self.refSkips.insert(i, gap)
							break
						else:
							self.refSkips[i] = self.refSkips[i] + 1

		# Map the new sequence onto the graph
		i = 0
		j = 0
		for base in newSeq:
			# If there is a gap in the graph at this position, but no gap in the seed sequence, move onto the next node
			while i in self.refSkips and seedSeq[i] != "-":
				self.graph[i]["-"][0] += 1
				i += 1

			self.graph[i][base][0] += 1
			# Obtain the maximum quality score for this base at this position
			if base != "-":
				if newQual[j] > self.graph[i][base][1]:
					self.graph[i][base][1] = newQual[j]
					j += 1
			i += 1

	def consensus(self):
		"""
		Returns a consensus sequence using the most common nodes in the De Brujn graph

		If multiple paths have the same weight, the path which coresponds to the reference will be used. If neither path is the reference, a random path will be chosen
		"""

		consensusSeq = ""
		consensusQual = []

		for node in self.graph:
			maxBase = None
			maxWeight = 0
			maxQual = 0
			for base, attributes in node.items():
				if base == "Ref":
					continue
				weight, qual = attributes
				if weight > maxWeight:
					maxBase = base
					maxWeight = weight
					maxQual = qual
				elif weight == maxWeight and base == node["Ref"]:
					maxBase = base
					maxQual = qual

			consensusSeq += maxBase
			consensusQual.append(maxQual)

		return consensusSeq, consensusQual
