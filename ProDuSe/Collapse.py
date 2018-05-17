#!/usr/bin/env python

import argparse
import pysam
import os
import sys
import sortedcontainers
import collections
import time
import bisect
from packaging import version
from skbio import alignment
from pyfaidx import Fasta
from configobj import ConfigObj


class Family:
    """
    Stores various statistics relating to a given read pair
    """

    def __init__(self, R1, R2, barcodeLength):

        # The first time this is initiated, this "family" will have a size of 1, since only a single
        # read pair is stored here.
        # The basis for this family is this initial read pair
        self.size = 1
        self.members = [R1.query_name]

        # Set counters which will store the number of sequences which contain elements (indels,
        # soft-clipping) that could mask INDELs
        self.R1abnormal = 0
        self.R2abnormal = 0

        #  Identify which read is actually read 1 and which is read2
        if not R1.is_read1:
            tmp = R1
            R1 = R2
            R2 = tmp

        # Identify the parental strand
        if R1.is_reverse:
            self.posParent = False
        else:
            self.posParent = True

        # Obtain the family name for this read pair (from the barcode)
        self.invalidBarcode = False
        try:
            self.familyName = R1.get_tag("OX")
            if len(self.familyName) != barcodeLength:
                raise TypeError("The read pair barcode \'%s\' for pair \'%s\' is not compatible with the specified mask" % (self.familyName, R1.query_name))
        except KeyError:
            self.familyName = None
            self.invalidBarcode = True

        # Do these families map to different chromosomes?
        self.isSplit = R1.reference_name != R2.reference_name

        try:
            # Convert pysam's cigar tuples into lists for easier iteration later
            self.R1cigar = self.cigarToList(R1.cigartuples, True)
            self.R2cigar = self.cigarToList(R2.cigartuples, False)
            self.R1start = R1.reference_start
            self.R2start = R2.reference_start
            self.softClipped = False  # Does this read pair display large swaths of soft-clipping which may support an alternative alignment?

            # We need to reverse the sequence of the reads that are mapped to the reverse strand, due to the read orientation
            if R1.is_reverse:
                self.R1sequence = [R1.query_sequence[::-1]]
                self.R1qual = [R1.query_qualities[::-1]]
                self.R1pos = R1.reference_end
                self.R1cigar = [self.R1cigar[::-1]]
                # Is this read soft clipped at the end? If so, we need to add the soft-clipping back onto the read position
                # to determine the "real" position of this read
                if R1.cigartuples[-1][0] == 4:
                    self.R1pos += R1.cigartuples[-1][1]
                    self.softClipped = True
                    # Skip realignment for now
                if R1.cigartuples[0][0] == 4:
                    self.R1start -= R1.cigartuples[0][1]
            else:
                self.R1sequence = [R1.query_sequence]
                self.R1qual = [R1.query_qualities]
                self.R1pos = R1.reference_start
                self.R1cigar = [self.R1cigar]
                if R1.cigartuples[0][0] == 4:
                    self.R1pos -= R1.cigartuples[0][1]
                    self.R1start = self.R1pos
                    self.softClipped = True

            if R2.is_reverse:
                self.R2sequence = [R2.query_sequence[::-1]]
                self.R2qual = [R2.query_qualities[::-1]]
                self.R2pos = R2.reference_end
                self.R2cigar = [self.R2cigar[::-1]]
                if R2.cigartuples[-1][0] == 4:
                    self.R2pos += R2.cigartuples[-1][1]
                    self.softClipped = True
                if R2.cigartuples[0][0] == 4:
                    self.R2start -= R2.cigartuples[0][1]
            else:
                self.R2sequence = [R2.query_sequence]
                self.R2qual = [R2.query_qualities]
                self.R2pos = R2.reference_start
                self.R2cigar = [self.R2cigar]
                if R2.cigartuples[0][0] == 4:
                    self.R2pos -= R2.cigartuples[0][1]
                    self.R2start = self.R2pos
                    self.softClipped = True

            # When we mark duplxes, the coordinates of read1 and read2 will be the other way around, since read1 will map to the
            # opposite strand
            # Thus, normalize the coordinates
            if self.R1pos > self.R2pos:
                tmp = self.R1pos
                self.R1pos = self.R2pos
                self.R2pos = tmp

            self.malformed = False
        except (IndexError, TypeError):
            self.malformed = True
            self.R1cigar = [None]
            self.R2cigar = [None]

        # These are statistics which will be used after all reads are collapsed into a consensus
        self.collapseSize = []

        # These are used to create output reads
        self.R1 = R1
        self.R2 = R2
        self.name = R1.query_name

    def _cigarstringToList(self, cigarString):

        cigToOp ={
            "S":4,
            "M":0,
            "I":2,  # This is the cigar of the query sequence, so it needs to be reversed
            "D":1,  # This is the cigar of the query sequence, so it needs to be reversed
            "H":5
        }
        cigar = []
        opLen = ""
        for char in cigarString:
            if char.isdigit():
                opLen += char
            else:
                cigar.extend([cigToOp[char]] * int(opLen))
                opLen = ""

        return cigar

    def _realign(self, read, windowBuffer=200):
        """
        Aligns a fragment of a sequence to the reference using Smith-Waterman alignment
        """

        # Load the current reference window and it's position
        global refWindow
        global refStart
        global refEnd
        global refName
        global refAlignment

        # If we have switched contigs, or moved outside the reference window that is currently buffered,
        # we need to re-buffer
        readWindowStart = read.reference_start - windowBuffer
        if readWindowStart < 0:
            readWindowStart = 0
        readWindowEnd = read.reference_start + windowBuffer
        if read.reference_name != refName or readWindowStart < refStart or readWindowEnd > refEnd:
            global refGenome
            # To reduce the number of times this needs to be performed, obtain a pretty large buffer
            refStart = readWindowStart - 500
            if refStart < 0:
                refStart = 0
            refEnd = readWindowStart + 4500
            refName = read.reference_name
            refWindow = refGenome[refName][refStart:refEnd].seq.upper()
            refAlignment = alignment.StripedSmithWaterman(query_sequence=refWindow)

        # Create a reference alignment object
        mapping = refAlignment(read.query_sequence)

        # If the read is mapped identically as the previous alignment, where will it start and end?
        expectedStart = read.reference_start - refStart
        expectedEnd = read.reference_end - refStart

        # Add soft clipping
        if mapping.target_begin != 0:
            finalCigar = self._cigarstringToList(str(mapping.target_begin) + "S" + mapping.cigar)
        else:
            finalCigar = self._cigarstringToList(mapping.cigar)

        endSoftClip = len(mapping.target_sequence) - mapping.target_end_optimal - 1
        if endSoftClip != 0:
            finalCigar.extend([4] * endSoftClip)

        # If the leading and trailing soft clipping have changed, quantify the change
        startOffset = mapping.query_begin - expectedStart - mapping.target_begin
        endOffset = mapping.query_end - expectedEnd + endSoftClip + 1

        return finalCigar, startOffset, endOffset

    def cigarToList(self, cigarTuples, isRead1, counterIncremeneted=False):
        """
        Convert a pysam-style cigar sequence (A list of tuples) into a list

        In addition, if this read contains leading or trailing soft clipping, or
        INDELs, increment the abnormal read counter, so we can realign any INDELs later on
        """
        cigarList = []
        for cigarElement in cigarTuples:
            cigOp = cigarElement[0]
            cigarList.extend([cigOp] * cigarElement[1])
            if not counterIncremeneted and (cigOp == 4 or cigOp == 1 or cigOp == 2):
                counterIncremeneted = True
                if isRead1:
                    self.R1abnormal += 1
                else:
                    self.R2abnormal += 1

        # I wish I didn't have to do this, but I can't guarantee that BWA (or any other aligner) will generate
        # a record with a valid cigar sequence, so we need to check that here
        if cigarList[0] == 1 or cigarList[0] == 2 or cigarList[-1] == 1 or cigarList[-1] == 2:
            raise TypeError("Invalid Cigar Sequence %s" % (cigarTuples))
        return cigarList

    def add(self, family, mismatchThreshold=8):
        """
        Combines another family into this family
        """

        if family.R1pos != self.R1pos or family.R2pos != self.R2pos:
            thisOrigEnd = None
            familyOrigEnd = None
            if self.R1.is_reverse:
                thisOrigEnd = self.R1.reference_end
            else:
                thisOrigEnd = self.R2.reference_end

            if family.R1.is_reverse:
                familyOrigEnd = family.R1.reference_end
            else:
                familyOrigEnd = family.R2.reference_end
            if thisOrigEnd != familyOrigEnd:
                return
            else:
                raise TypeError()

        self.R1sequence.extend(family.R1sequence)
        self.R1qual.extend(family.R1qual)
        self.R1cigar.extend(family.R1cigar)

        self.R2sequence.extend(family.R2sequence)
        self.R2qual.extend(family.R2qual)
        self.R2cigar.extend(family.R2cigar)

        self.size += family.size
        self.members.extend(family.members)

    def consensus(self, abnormalityThresold=0.4):
        """
        Merge all sequences and quality scores stored in this read pair into a consensus
        """

        # Ok, I have a confession to make
        # This "Family" object can store more than one read pair
        # And I think it's about time we fixed that

        # First, because of some cigar sequence BS, we need to account for the fact that
        # some of these families may have some leading soft clipping
        # The easiest way to account for this is to simply figure out what the offset of each
        # sequence is

        R1consensus, R1qual, R1cigar, R1softClip = self._consensusByRead(self.R1sequence, self.R1qual, self.R1cigar, self.R1.is_reverse)
        R2consensus, R2qual, R2cigar, R2softClip = self._consensusByRead(self.R2sequence, self.R2qual, self.R2cigar, self.R2.is_reverse)
        self.R1sequence = [R1consensus]
        self.R2sequence = [R2consensus]
        self.R1qual = [R1qual]
        self.R2qual = [R2qual]
        self.R1cigar = [R1cigar]
        self.R2cigar = [R2cigar]
        self.R1start += R1softClip
        self.R2start += R2softClip  # Compensate for any changes in soft-clipping

    def _consensusByRead(self, sequences, qualities, cigars, reverseClip=False):

        # Generate a consensus for all reads in the family in the provided orientation

        consensusSeq = []
        consensusQual = []
        consensusCigar = []

        # Handle leading anjd trailing soft-clipping
        softClippedStart = True
        softClippedEnd = False

        # First, shift all the sequences over to account for any leading soft-clipping
        # and differences in start position
        # The end goal of this is to have the bases of all sequences line up properly
        seqNumber = len(sequences)
        baseIndex = [0] * seqNumber
        cigarIndex = [0] * seqNumber
        seqLengths = tuple(len(x) for x in cigars)
        startOffset = 0
        refLength = 0  # The change in the number of reference positions consumed by this read relative to read1

        # Lets start processing the sequence
        while True:
            qualSum = {"A": 0, "C": 0, "T": 0, "G": 0, "-": 0, "N": 0}
            qualMax = {"A": 0, "C": 0, "T": 0, "G": 0, "-": 0, "N": 0}
            insertion = []  # This data structure will store {Base: (count, maxqual, qualSum)}
            normBaseCount = {"A": 0, "C": 0, "T": 0, "G": 0, "-": 0, "N": 0}
            softClippedBaseCount = {"A": 0, "C": 0, "T": 0, "G": 0, "N": 0}

            seqCollapsed = 0.0
            softClippedNum = 0.0

            for i in range(0, seqNumber):
                if cigarIndex[i] < seqLengths[i]:
                #try:
                    # Here we have to account for the cigar operator
                    cigar = cigars[i][cigarIndex[i]]
                    if cigar == 4:
                        softClippedNum += 1

                    # If this sequence has an insertion at this position,
                    # We need to account for it, and store it seperatly for comparison
                    j = 0
                    while cigar == 1:
                        try:
                            if i == 0:
                                refLength += 1
                            base = sequences[i][baseIndex[i]]
                            qual = qualities[i][baseIndex[i]]

                            if j >= len(insertion):
                                insertion.append({base: [1, qual, qual]})
                            elif base not in insertion[j]:
                                insertion[j][base] = [1, qual, qual]
                            else:
                                insertion[j][base][0] += 1  # Base count
                                insertion[j][base][2] += qual  # Qual sum
                                if insertion[j][base][1] < qual:  # Max qual
                                    insertion[j][base][1] = qual
                            cigarIndex[i] += 1
                            baseIndex[i] += 1
                            j += 1
                            cigar = cigars[i][cigarIndex[i]]

                        except IndexError as e: # This record ends with an insertion, which is not valid
                            badCigar = self.listToCigar(cigars[i])
                            raise TypeError("Invalid Cigar Sequence %s" % badCigar) from e

                    # In the case of a deletion, all we need to do is indicate that there is
                    # a deletion at this position
                    if cigar == 2:
                        normBaseCount["-"] += 1
                        cigarIndex[i] += 1
                    # Otherwise, save the base and associated quality score
                    else:
                        base = sequences[i][baseIndex[i]]
                        qual = qualities[i][baseIndex[i]]

                        if cigar == 4:
                            softClippedBaseCount[base] += 1
                        else:
                            normBaseCount[base] += 1
                        qualSum[base] += qual
                        if qualMax[base] < qual:
                            qualMax[base] = qual
                        cigarIndex[i] += 1
                        baseIndex[i] += 1
                    seqCollapsed += 1

                #except IndexError:
                #    pass

            if seqCollapsed / seqNumber < 0.5:
                break

            # Identify the consensus base/sequence at this position
            # The base which is most frequent is used as the consensus
            maxBase = None
            maxBaseCount = 0
            maxQual = 0

            if softClippedStart and softClippedNum / seqCollapsed <= 0.5:
                softClippedStart = False
            elif not softClippedStart and not softClippedEnd and softClippedNum / seqCollapsed > 0.5:
                softClippedEnd = True

            if softClippedStart or softClippedEnd:
                for base, count in softClippedBaseCount.items():

                    if count > maxBaseCount or (count == maxBaseCount and maxQual < qualSum[base]):  # If there is a tie, use the base with the highest aggregate quality score

                        maxBase = base
                        maxBaseCount = count
                        maxQual = qualSum[base]
                maxCigOp = 4
                if softClippedStart and not reverseClip:
                    startOffset += 1
                elif softClippedEnd and reverseClip:
                    startOffset += 1

            else:
                for base, count in normBaseCount.items():

                    if count > maxBaseCount or (count == maxBaseCount and maxQual < qualSum[base]):  # If there is a tie, use the base with the highest aggregate quality score

                        maxBase = base
                        maxBaseCount = count
                        maxQual = qualSum[base]

                if maxBase == "-":
                    maxCigOp = 2
                else:
                    maxCigOp = 0

            # Next, determine if the insertion is actually more common than the current base
            # TODO: Check if this breaks anything. Because it looks like it might
            insertionWeight = sum(insertion[0][x][0] for x in insertion[0]) if len(insertion) != 0 else -1
            if insertionWeight > seqCollapsed/ 2:

                # If so, we will simply use the consensus of the insertion
                for position in insertion:
                    posWeight = sum(position[x][0] for x in position)
                    if posWeight < insertionWeight / 2:
                        break

                    maxPosBase = None
                    maxPosCount = 0
                    maxPosQual = 0
                    for base, elements in position.items():
                        if elements[0] > maxPosCount:
                            maxPosBase = base
                            maxPosCount = elements[0]
                            maxPosQual = elements[2]
                        # If there is a tie, use the base that is not soft clipped, or that with the highest aggregate quality score
                        elif elements[0] == maxPosCount and maxPosQual < qualSum[base]:
                            maxPosBase = base
                            maxPosCount = elements[0]
                            maxPosQual = elements[2]

                    # Add these bases, quality scores, and the appropriate cigar operator
                    # into the consensus sequence
                    consensusSeq.append(maxPosBase)
                    consensusQual.append(position[base][1])
                    consensusCigar.append(1)

            # Add this base, quality score, and cigar operator into the consensus
            if maxCigOp != 2:
                consensusSeq.append(maxBase)
                consensusQual.append(qualMax[maxBase])
                consensusCigar.append(maxCigOp)
            else:
                consensusCigar.append(2)
            refLength += 1

        # To handle the extremely rare cases which may cause a consensus read to start or end with a INDEL
        # Remove such events from the start or end of the read
        # From the start of the read
        i = 0
        try:
            while True:
                op = consensusCigar[i]
                if op == 0:
                    break
                elif op == 2:  # Remove this deletion from the start of the read
                    consensusCigar.pop(i)
                    startOffset += 1
                elif op == 1:  # Remove this insertion from the sequence, quality, and cigar
                    consensusCigar.pop(i)
                    consensusSeq.pop(i)
                    consensusQual.pop(i)
                else:
                    i += 1
        except IndexError:  # The entire read is soft clipped
            pass
        # From the back of the read
        i = -1
        try:
            while True:
                op = consensusCigar[i]
                if op == 0:
                    break
                elif op == 2:  # Remove this deletion from the start of the read
                    consensusCigar.pop(i)
                    refLength -= 1
                elif op == 1:  # Remove this insertion from the sequence, quality, and cigar
                    consensusCigar.pop(i)
                    consensusSeq.pop(i)
                    consensusQual.pop(i)
                else:
                    i -= 1
        except IndexError:  # The entire read is soft clipped
            pass

        # If this read is mapped to the reverse strand, and thus we are collapsing from the back, account for the
        # case where the start position changes due to differences between the consensus length and the "seed" read length
        if reverseClip:
            startOffset += len(cigars[0]) - refLength

        return "".join(consensusSeq), consensusQual, consensusCigar, startOffset

    def listToCigar(self, cigar):
        """
        Converts a list of cigar operators into a pysam-style cigar sequence
        """

        cigarTuples = []
        currentOp = cigar[0]
        opLength = 0
        for op in cigar:
            # If this is a new cigar operator, we need to add the previous operator to the output
            # tuple
            if op != currentOp:
                cigarTuples.append((currentOp, opLength))
                currentOp = op
                opLength = 1
            else:
                opLength += 1

        cigarTuples.append((currentOp, opLength))

        # Bug-hunting here. Find invalid cigar records
        #for i in range(1, len(cigarTuples) - 2):
        #    if cigarTuples[i][0] == 4:
        #        raise TypeError
        return cigarTuples

    def toPysam(self, tagOrig):
        """
        Convers this read pair into two pysam.AlignedSegments
        """

        self.R1.query_name = self.name
        self.R1.reference_start = self.R1start
        self.R2.query_name = self.name
        self.R2.reference_start = self.R2start
        if tagOrig:
            self.R1.set_tag("Zm", ",".join(self.members))
            self.R2.set_tag("Zm", ",".join(self.members))

        # Remove the barcode tag, as it is no longer required
        self.R1.set_tag("OX", None)
        self.R2.set_tag("OX", None)

        # Un-reverse the sequence of the reverse-strand-mapped, so the sequences will be in the expected format
        if self.R1.is_reverse:
            self.R1.query_sequence = self.R1sequence[0][::-1]
            self.R1.query_qualities = self.R1qual[0][::-1]
            self.R1.cigartuples = self.listToCigar(self.R1cigar[0][::-1])
        else:
            self.R1.query_sequence = self.R1sequence[0]
            self.R1.query_qualities = self.R1qual[0]
            self.R1.cigartuples = self.listToCigar(self.R1cigar[0])

        if self.R2.is_reverse:
            self.R2.query_sequence = self.R2sequence[0][::-1]
            self.R2.query_qualities = self.R2qual[0][::-1]
            self.R2.cigartuples = self.listToCigar(self.R2cigar[0][::-1])
        else:
            self.R2.query_sequence = self.R2sequence[0]
            self.R2.query_qualities = self.R2qual[0]
            self.R2.cigartuples = self.listToCigar(self.R2cigar[0])

        self.R1.next_reference_start = self.R2.reference_start
        self.R2.next_reference_start = self.R1.reference_start

        # Update the mate cigar string
        self.R1.set_tag("MC", self.R2.cigarstring)
        self.R2.set_tag("MC", self.R1.cigarstring)


class Position:
    """
    A data structure which contains all the read pairs which start at a given position
    """

    def __init__(self):
        # A dictionary for each parental strand
        self.plusFamilies = {}
        self.negFamilies = {}

    def add(self, readPair):
        """
        Store the specified read pair at this position

        This read pair will be stored in a dictionary containing barcode: readPair,
        which coresponds to the mapping strand
        """

        if readPair.posParent:  # i.e. it originated from a "+" strand family

            # If an identical barcode is already stored in the dictionary, then we can
            # be sure that the existing read pair and the new read pair originate from the same
            # family
            # In this case, simply aggregate them
            if readPair.familyName in self.plusFamilies:
                self.plusFamilies[readPair.familyName].add(readPair)
            else:
                self.plusFamilies[readPair.familyName] = readPair

        else:
            # If an identical barcode is already stored in the dictionary, then we can
            # be sure that the existing read pair and the new read pair originate from the same
            # family
            # In this case, simply aggregate them
            if readPair.familyName in self.negFamilies:
                self.negFamilies[readPair.familyName].add(readPair)
            else:
                self.negFamilies[readPair.familyName] = readPair

    def addPosition(self, otherPosition):
        """
        Merges the read pairs stored at two positions

        In addition, keeps track of which read pairs have been merged
        """

        for familyName, readPair in otherPosition.plusFamilies.items():

            # If an identical barcode is already stored in the dictionary, then we can
            # be sure that the existing read pair and the new read pair originate from the same
            # family
            # In this case, simply aggregate them
            if familyName in self.plusFamilies:
                self.plusFamilies[familyName].add(readPair)
            else:
                self.plusFamilies[familyName] = readPair

        for familyName, readPair in otherPosition.negFamilies.items():
            # If an identical barcode is already stored in the dictionary, then we can
            # be sure that the existing read pair and the new read pair originate from the same
            # family
            # In this case, simply aggregate them
            if familyName in self.negFamilies:
                self.negFamilies[familyName].add(readPair)
            else:
                self.negFamilies[familyName] = readPair

    def collapse(self, collapseIndices, collapseThreshold=3):
        """
        Identify read pairs which are duplicates, and merge them into a consensus
        """

        # Alright ladies and gentlemen, now that we have done all the boring organizing tasks, it's
        # time to do what we all came here to do

        # Lets, collapse the (+) strand families first
        # Obtain a list of all adapter sequences which occur at this position, along with how frequently that
        # barcode occurs

        barcodesInFamilies = collections.OrderedDict()  # Store the base barcode of the family, as well as all barcodes which should be collapsed into that family
        barcodesByFrequency = sortedcontainers.SortedListWithKey(list(self.plusFamilies.keys()),
                                                                 key=lambda x: self.plusFamilies[x].size)
        for barcode in reversed(barcodesByFrequency):

            # If adaptersInFamily is empty, then no barcodes have been processed as of yet
            # Since the first barcode returned is the largest, we'll consider that the name of the first family
            # Since it almost certain that it would end up as the family name anyways
            if len(barcodesInFamilies) == 0:
                barcodesInFamilies[barcode] = [barcode]
                continue

            # Consider only the specified indexes
            cBarcode = list(barcode[x] for x in collapseIndices)
            # Next, calculate the distance between the current barcode and any existing barcode family
            minBarcode = None
            minDistance = 10000
            for familyName in barcodesInFamilies.keys():
                distance = 0
                cFamilyName = list(familyName[x] for x in collapseIndices)
                for b1, b2 in zip(cFamilyName, cBarcode):
                    if b1 != b2:
                        distance += 1
                if distance < minDistance:
                    # In the case of a tie, the largest family will be taken
                    minDistance = distance
                    minBarcode = familyName

            # Next, check if the distance between this barcode and the nearest barcode is within threshold
            # If so, we can consider them as derrived from the same parental molecule, and group them
            if minDistance <= collapseThreshold:
                barcodesInFamilies[minBarcode].append(barcode)
            # Otherwise, a new family needs to be constructed for this barcode, as it is too distant from any existing
            # barcodes
            else:
                barcodesInFamilies[barcode] = [barcode]

        # Now that each barcode has been assigned to a family, we need to actually
        # do the deed, and collapse all read pairs in a given family
        tmpPlusFamilies = {}
        for familyName, members in barcodesInFamilies.items():
            # Choose the first read pair encounter as the "template", to which all
            # other read pairs that need to be collapsed will be added
            consensusPair = None
            for member in members:

                if consensusPair is None:
                    consensusPair = self.plusFamilies[member]
                    continue
                consensusPair.add(self.plusFamilies[member])

            # Finally, collapse the consensus read into an actual consensus
            consensusPair.consensus()
            tmpPlusFamilies[familyName] = consensusPair

        # Finally, to save memory, delete families which were collapsed into the consensus families
        self.plusFamilies = tmpPlusFamilies

        # Now repeat everything for (-) strand families
        # Obtain a list of all adapter sequences which occur at this position, along with how frequently that
        # barcode occurs
        barcodesInFamilies = collections.OrderedDict()  # Store the base barcode of the family, as well as all barcodes which should be collapsed into that family
        barcodesByFrequency = sortedcontainers.SortedListWithKey(list(self.negFamilies.keys()),
                                                                 key=lambda x: self.negFamilies[x].size)
        for barcode in reversed(barcodesByFrequency):

            # If adaptersInFamily is empty, then no barcodes have been processed as of yet
            # Since the first barcode returned is the largest, we'll consider thta the name of the first family
            # Since it almost certain that it would end up as the family name anyways
            if len(barcodesInFamilies) == 0:
                barcodesInFamilies[barcode] = [barcode]
                continue

            cBarcode = list(barcode[x] for x in collapseIndices)
            # Next, calculate the distance between the current barcode and any existing barcode family
            minBarcode = None
            minDistance = 10000
            for familyName in barcodesInFamilies.keys():
                distance = 0
                cFamilyName = list(familyName[x] for x in collapseIndices)
                for b1, b2 in zip(cFamilyName, cBarcode):
                    if b1 != b2:
                        distance += 1
                if distance < minDistance:
                    # In the case of a tie, the most largest family will be taken
                    minDistance = distance
                    minBarcode = familyName

            # Next, check if the distance between this barcode and the nearest barcode is within threshold
            # If so, we can consider them as derrived from the same parental molecule, and group them
            if minDistance <= collapseThreshold:
                barcodesInFamilies[minBarcode].append(barcode)
            # Otherwise, a new family needs to be constructed for this barcode, as it is too distant from any existing
            # barcodes
            else:
                barcodesInFamilies[barcode] = [barcode]

        # Now that each barcode has been assigned to a family, we need to actually
        # do the deed, and collapse all read pairs in a given family
        tmpNegFamilies = {}
        for familyName, members in barcodesInFamilies.items():
            # Choose the first read pair encounter as the "template", to which all
            # other read pairs that need to be collapsed will be added
            consensusPair = None
            for member in members:

                if consensusPair is None:
                    consensusPair = self.negFamilies[member]
                    continue
                consensusPair.add(self.negFamilies[member])

            # Finally, collapse the consensus read into an actual consensus
            consensusPair.consensus()
            tmpNegFamilies[familyName] = consensusPair
        self.negFamilies = tmpNegFamilies

    def markDuplexes(self, duplexIndices, duplexDistance=2):
        """
        Identify families which exist in a duplex

        """

        global counter
        # Here we are simply going to examine each family which originates from a (-) strand molecule
        # and try to find a coresponding (+) strand family
        processedPlusFamilies = {}
        for key, readPair in self.negFamilies.items():

            # Since these collections are negative, we need to flip the adapter sequences, as they are currently reversed.
            # This was not necessary for collapsing, since all molecules derrived from the same parental strand had the adapter
            # Sequences in the same orientation, but since we are now comparing between (+) and (-) strand families, this
            # will need to occur
            lKey = len(key)
            adapter = key[int(lKey / 2):] + key[:int(lKey / 2)]

            # Assign this family a unique name. If a (+) strand family is in duplex, it will be assigned a complementary name later
            name = adapter + ":-:" + str(readPair.size) + ":" + str(counter)
            readPair.name = name

            # If there are no families which originate from the (+) parental strand, just assign the current (-) strand families
            # a unique name
            if len(self.plusFamilies) == 0:
                counter += 1
                continue

            # Otherwise, we need to determine if any of families which originate from the parental (+) strand could originate from the
            # same parental molecule, using the adapter sequence
            # First, calculate the distance between the (-) strand family, and all possible (+) strand families, and find the (+) strand family
            # which is closest
            dKey = tuple(adapter[x] for x in duplexIndices)
            minDist = 100000
            minAdapter = None
            minSize = 0
            for familyName in self.plusFamilies:
                # Calculate the distance between these two family names
                dFamilyName = tuple(familyName[x] for x in duplexIndices)
                distance = 0
                for base1, base2 in zip(dKey, dFamilyName):
                    if base1 != base2:
                        distance += 1

                # If this family is closer than the current closest family, replace it
                # If there is a tie, use the largest family size (and if that is a tie, take the first one that is encountered)
                if distance < minDist or (distance == minDist and self.plusFamilies[familyName].size > minSize):
                    minDist = distance
                    minAdapter = familyName
                    minSize = self.plusFamilies[familyName].size

            # If the (-) strand family and closest (+) strand family are within the specified mismatch threshold, they are considered to be
            # in duplex
            if minDist <= duplexDistance:
                # In this case, ensure that the names of the families are given a complimentary names
                # so they can be identified as in duplex
                duplexPair = self.plusFamilies.pop(minAdapter)
                duplexPair.name = adapter + ":+:" + str(duplexPair.size) + ":" + str(counter)
                processedPlusFamilies[minAdapter] = duplexPair

            counter += 1

        # We have finished processing all the families which originate from the (+) parental strand
        # Any (+) families remaining must be derrived from unique molecules, as they were not paired with a (-) family
        # while processing those
        # In this case, just assign each an apropriate name
        for adapter in list(self.plusFamilies.keys()):
            readPair = self.plusFamilies.pop(adapter)
            readPair.name = adapter + ":+:" + str(readPair.size) + ":" + str(counter)
            processedPlusFamilies[adapter] = readPair
            counter += 1
        self.plusFamilies = processedPlusFamilies


class FamilyCoordinator:
    """
    Processes reads from a given BAM file, identifies duplicates, and collapses duplicates into families
    """

    def __init__(self, inputFile, reference, familyIndices, familyThreshold, duplexIndices, duplexThreshold,
                 barcodeLength, targets=None, tagOrig = False, baseBuffer=400, padding=10, printPrefix="PRODUSE-COLLAPSE"):
        self.inFile = inputFile
        self.tagOrig = tagOrig

        self.readCounter = 0
        self.pairCounter = 0
        self.familyCounter = 0
        self.malformedCigar = 0
        self.failedQC = 0
        self.missingBarcode = 0
        self.outsideCaptureSpace = 0

        # Post-collapse stats
        self.familyDistribution = {}
        self.duplexDistribution = {}
        self.depthDistribution = {}
        global refGenome
        refGenome = self._loadContigs(reference)
        global refWindow
        refWindow = None
        global refStart
        refStart = None
        global refEnd
        refEnd = None
        global refName
        refName = None

        self.targets = self._loadTargets(targets, padding)

        self.familyIndices = familyIndices
        self.familyThreshold = familyThreshold
        self.duplexIndices = duplexIndices
        self.duplexThreshold = duplexThreshold
        self.barcodeLength = barcodeLength

        self._waitingForMate = {}
        self._pairsAtPositions = sortedcontainers.SortedDict()

        self._previousPos = -1000
        self._previousChr = None
        self._baseBuffer = baseBuffer
        global counter
        counter = 0

        self.printPrefix = printPrefix

    def _loadTargets(self, targetFile, padding):
        """
        Loads BED intervals from the specified file into a dictionary

        :param targetFile: A string representing a filepath to a BED file
        """

        targets = {}

        # If no file was specified, return an empty dictionary
        if targetFile is None:
            return targets

        # Does the BAM file use "chr"-prefixed contig names? If so, we should make sure the names of these regions are consistent
        isChrPrefixed = None
        try:
            # Check the first reference sequence listed in the BAM header. This isn't perfect, but it should handle most cases
            isChrPrefixed = self.inFile.header["SQ"][0]["SN"].startswith("chr")  # TODO: Handle this better
        except AssertionError:
            pass
        try:
            with open(targetFile) as f:
                for line in f:
                    elements = line.split("\t")

                    chrom = elements[0]
                    # Make sure the "chr" prefix is consistent with the BAM file
                    if isChrPrefixed is not None:
                        if isChrPrefixed:
                            if not chrom.startswith("chr"):
                                chrom = "chr" + chrom
                        else:
                            if chrom.startswith("chr"):
                                chrom = chrom[3:]

                    start = int(elements[1]) - padding
                    end = int(elements[2]) + padding
                    if chrom not in targets:
                        targets[chrom] = []
                    targets[chrom].extend([start, end])
                    # Ensure that the BED interval is actually valid
                    if start > end:
                        sys.stderr.write(
                            "ERROR: The start position \'%s\' is greater than the end position \'%s\'\n" % (start, end))
                        exit(1)

        except IndexError:
            sys.stderr.write("ERROR: Unable to parse BED interval \'%s\'\n" % (line.strip("\n").strip("\r")))
            sys.stderr.write("Ensure the file is a valid BED file, and try again\n")
            exit(1)
        except ValueError:
            sys.stderr.write("ERROR: Unable to parse BED interval \'%s\'\n" % (line.strip("\n").strip("\r")))
            sys.stderr.write("Ensure the file is a valid BED file, and try again\n")
            exit(1)

        # Finally, to speed up access, convert each list into a set
        for chrom, locations in targets.items():
            targets[chrom] = tuple(locations)
        return targets

    def _loadContigs(self, refGenome):
        """
        Prepares contigs from the specified file

        :args refGenome: A string containing a filepath to the reference FASTA file

        :returns:
        """
        return Fasta(refGenome, one_based_attributes=False, rebuild=False)

    def _withinCaptureSpace(self, readPair):  # TODO: Fix cases where the read pair completely overlaps the capture space
        """
        Determines if a read pair overlaps the capture space

        :param readPair: A Family() object
        :return: If at least one read falls within the capture space
        """

        if readPair.R1.reference_name in self.targets:
            bisectPoint1 = bisect.bisect_left(self.targets[readPair.R1.reference_name], readPair.R1.reference_start)
            bisectPoint2 = bisect.bisect_left(self.targets[readPair.R1.reference_name], readPair.R1.reference_end)
            if bisectPoint1 != bisectPoint2 or bisectPoint1 % 2 == 1:
                return True

        if readPair.R2.reference_name in self.targets:
            bisectPoint1 = bisect.bisect_left(self.targets[readPair.R2.reference_name], readPair.R2.reference_start)
            bisectPoint2 = bisect.bisect_left(self.targets[readPair.R2.reference_name], readPair.R2.reference_end)
            if bisectPoint1 != bisectPoint2 or bisectPoint1 % 2 == 1:
                return True
        return False

    def __next__(self):
        return self.__iter__()

    def __iter__(self):

        sys.stderr.write("\t".join([self.printPrefix, time.strftime('%X'), "Starting...\n"]))

        try:
            while True:
                read = next(self.inFile)

                # Discard supplementary and secondary alignments
                if read.is_supplementary or read.is_secondary:
                    continue
                self.readCounter += 1
                # If this read and it's mate are unmapped, ignore it
                try:
                    currentPos = read.reference_start
                    currentChrom = read.reference_name
                except ValueError:
                    continue

                # If the start position of this read is not the same as the start position of the
                # previous read, and assuming the input BAM file is sorted, then we can assume that
                # no additional reads exists which start at the previous position
                # Thus, let's process the previous position
                #
                # To account for soft-clipped reads (where the start position may not be accurate),
                # include a reasonable offset before we attempt to collapse
                if currentChrom != self._previousChr or currentPos != self._previousPos + self._baseBuffer:
                    # As a sanity check, ensure that the new position is greater than the previous position
                    # (unless we have switched chromosomes)
                    # If it's not, then the input BAM file is unsorted
                    assert self._previousPos + self._baseBuffer < currentPos or currentChrom != self._previousChr

                    # Identify all previous positions that are to be processed
                    # If we have switched chromosomes, purge everything, as no more reads are coming which map to the
                    # previous chromosome. Thus, all families are complete
                    if currentChrom != self._previousChr:
                        i = len(self._pairsAtPositions) - 1
                    else:
                        # Otherwise, we need to determine which positions are safe to process
                        i = self._pairsAtPositions.bisect(currentPos - self._baseBuffer) - 1
                    while i >= 0:
                        posKey = self._pairsAtPositions.iloc[i]
                        posList = self._pairsAtPositions.pop(posKey)

                        # Now, it's time to do all the magic
                        # Identify any reads which could originate from the same family, and collapse them into a consensus
                        for posToProcess in posList.values():
                            posToProcess.collapse(self.familyIndices, self.familyThreshold)
                            posToProcess.markDuplexes(self.duplexIndices, self.duplexThreshold)

                            # Return all families stored at this position
                            for readPair in posToProcess.plusFamilies.values():
                                readPair.toPysam(self.tagOrig)
                                yield readPair.R1
                                yield readPair.R2
                            for readPair in posToProcess.negFamilies.values():
                                readPair.toPysam(self.tagOrig)
                                yield readPair.R1
                                yield readPair.R2
                        i -= 1

                    self._previousChr = currentChrom
                    self._previousPos = currentPos - self._baseBuffer

                # Have we seen this read's mate?
                if read.query_name not in self._waitingForMate:
                    # If not, lets buffer this read until we encounter it's mate
                    self._waitingForMate[read.query_name] = read
                    continue

                # First, create an object representing this read pair
                pair = Family(read, self._waitingForMate[read.query_name], self.barcodeLength)
                self.pairCounter += 1

                if self.pairCounter % 100000 == 0:
                    sys.stderr.write(
                        "\t".join(
                            [self.printPrefix, time.strftime('%X'), "Pairs Processed:" + str(self.pairCounter) + "\n"]))

                # Delete the mate from the dictionary to free up space
                del self._waitingForMate[read.query_name]

                # Perform some basic QC
                # Is this read pair missing a cigar string? If so, don't process it
                if pair.malformed:
                    self.malformedCigar += 1
                    continue

                # Is this read pair missing a barcode tag? If so, we can't process it, as we
                # won't be able to find out which family it belongs to
                if pair.invalidBarcode:
                    self.missingBarcode += 1
                    continue

                # If this read pair is split(i.e. the reads map to different chromosomes), then it's very likely that
                # one of the existing read's positions was processed a long time ago. In which case, we can no longer
                # collapse it
                if pair.isSplit:
                    continue

                # if this read paif falls outside the capture space completely, don't process it
                if self.targets:
                    withinCapture = self._withinCaptureSpace(pair)
                    if not withinCapture:
                        self.outsideCaptureSpace += 1
                        continue

                # If this read pair contains leading soft clipping, then the start position of
                # this read pair may not be accurate (due to possible leading insertions or deletions)
                # We need to set it aside for now, and try to identify any families which it
                # originates from which do not contain soft clipping when we are collapsing
                # if pair.softClipped:
                # 	continue

                # Store this read pair until we obtain all read pairs that overlap this position
                # If this is the first time we are seeing a read pair at this position,
                # we need to create the data structure which will hold all the read pairs
                if pair.R1pos not in self._pairsAtPositions:
                    self._pairsAtPositions[pair.R1pos] = {}
                    self._pairsAtPositions[pair.R1pos][pair.R2pos] = Position()
                elif pair.R2pos not in self._pairsAtPositions[pair.R1pos]:
                    self._pairsAtPositions[pair.R1pos][pair.R2pos] = Position()
                self._pairsAtPositions[pair.R1pos][pair.R2pos].add(pair)

        except StopIteration:

            # We have run out of reads. Thus, finalize and collapse any remaining positions
            for posList in self._pairsAtPositions.values():

                # Now, it's time to do all the magic
                # Identify any reads which could originate from the same family, and collapse them into a consensus
                for posToProcess in posList.values():
                    posToProcess.collapse(self.familyIndices, self.familyThreshold)
                    posToProcess.markDuplexes(self.duplexIndices, self.duplexThreshold)

                    # Return all remaining families
                    for readPair in posToProcess.plusFamilies.values():
                        readPair.toPysam(self.tagOrig)
                        yield readPair.R1
                        yield readPair.R2
                    for readPair in posToProcess.negFamilies.values():
                        readPair.toPysam(self.tagOrig)
                        yield readPair.R1
                        yield readPair.R2

            # Print out a status (or error) message, briefly summarizing the overall collapse
            if self.missingBarcode > 0 and self.familyCounter == 0:
                sys.stderr.write("ERROR: Unable to find a \'OX\' tag, which contains the degenerate barcode, for any read in the input BAM file\n")
                sys.stderr.write(
                    "Check that BWA was run using the \'-C\' option\n")

            elif self.readCounter == 0:
                sys.stderr.write("ERROR: The input BAM file is empty!")
            else:
                sys.stderr.write(
                    "\t".join([self.printPrefix, time.strftime('%X'), "Pairs Processed:" + str(self.pairCounter) + "\n"]))
                sys.stderr.write("\t".join([self.printPrefix, time.strftime('%X'), "Collapse Complete\n"]))


def validateArgs(args):
    """
    Validates that the specified set of arguments are valid
    :param args: A dictionary listing {argument: parameter}
    :return: A dictionary listing {argument: parameter} that have been validatd
    """

    # Convert the dictionary into a list, to allow the arguments to be re-parsed
    listArgs = []
    for argument, parameter in args.items():

        if parameter is None or parameter is False or parameter == "None" or parameter == "False":
            continue
        # Something was provided as an argument
        listArgs.append("--" + argument)

        # If the parameter is a boolean, ignore it, as this will be reset once the arguments are re-parsed
        if parameter == "True" or parameter is True:
            continue
        # If the parameter is a list, we need to add each element seperately
        if isinstance(parameter, list):
            for p in parameter:
                listArgs.append(str(p))
        else:
            listArgs.append(str(parameter))

    parser = argparse.ArgumentParser(description="Collapsed duplicate sequences into a consensus")
    parser.add_argument("-c", "--config", type=lambda x: isValidFile(x, parser),
                        help="An optional configuration file, which can provide one or more arguments")
    parser.add_argument("-i", "--input", type=lambda x: isValidFile(x, parser, True), required=True,
                        help="Input sorted BAM file (use \"-\" to read from stdin [use control + d to stop reading])")
    parser.add_argument("-o", "--output", required=True,
                        help="Output BAM file (use \"-\" to write to stdout). Will be unsorted")
    parser.add_argument("-fm", "--family_mask", required=True, type=str,
                        help="Positions in the barcode to consider when collapsing reads into a consensus (1=Consider, 0=Ignore)")
    parser.add_argument("-dm", "--duplex_mask", required=True, type=str,
                        help="Positions in the barcode to consider when determining if two families are in duplex (1=Consider, 0=Ignore)")
    parser.add_argument("-fmm", "--family_mismatch", required=True, type=int,
                        help="Maximum number of mismatches permitted when collapsing reads into a family")
    parser.add_argument("-dmm", "--duplex_mismatch", required=True, type=int,
                        help="Maximum number of mismatches permitted when identifying of two families are in duplex")
    parser.add_argument("--tag_family_members", action="store_true",
                        help="Store the names of all reads used to generate a consensus in the tag \'Zm\'")
    parser.add_argument("-r", "--reference", required=True, type=lambda x: isValidFile(x, parser),
                        help="Reference genome, in FASTA format")
    parser.add_argument("-t", "--targets", type=lambda x: isValidFile(x, parser),
                        help="A BED file listing targets of interest")
    validatedArgs = parser.parse_args(listArgs)
    return vars(validatedArgs)


def isValidFile(file, parser, allowStream=False):
    """
    Checks to ensure the provided file is exists, and throws an error if it is not.

    A UNIX pipe ("-") is also considered valid

    :param file: A string containing a filepath to the file of interest (or a pipe symbol)
    :param parser: An argparse.ArgumentParser() object.
    :param allowStream: If a UNIX PIPE symbol ("-"), is permitted

    :returns: The "file" variable, if the file is valid
    :raises parser.error: If the file is not valid
    """

    if file == "-" and allowStream:  # i.e. they are specifying a pipe
        return file
    elif os.path.exists(file):
        return file
    else:
        raise parser.error("Unable to locate %s. Please ensure the file exists, and try again." % (file))


# Process command line arguments
parser = argparse.ArgumentParser(description="Collapsed duplicate sequences into a consensus")
parser.add_argument("-c", "--config", type=lambda x: isValidFile(x, parser),
                    help="An optional configuration file, which can provide one or more arguments")
parser.add_argument("-i", "--input", type=lambda x: isValidFile(x, parser, True),
                    help="Input sorted BAM file (use \"-\" to read from stdin [and \"Control + D\" to stop reading])")
parser.add_argument("-o", "--output", help="Output BAM file (use \"-\" to write to stdout). Will be unsorted")
parser.add_argument("-fm", "--family_mask", type=str,
                    help="Positions in the barcode to consider when collapsing reads into a consensus (1=Consider, 0=Ignore)")
parser.add_argument("-dm", "--duplex_mask", type=str,
                    help="Positions in the barcode to consider when determining if two families are in duplex (1=Consider, 0=Ignore)")
parser.add_argument("-fmm", "--family_mismatch", type=int,
                    help="Maximum number of mismatches permitted when collapsing reads into a family")
parser.add_argument("-dmm", "--duplex_mismatch", type=int,
                    help="Maximum number of mismatches permitted when identifying of two families are in duplex")
parser.add_argument("--tag_family_members", action="store_true", help="Store the names of all reads used to generate a consensus in the tag \'Zm\'")
parser.add_argument("-r", "--reference", type=lambda x: isValidFile(x, parser),
                    help="Reference genome, in FASTA format")
parser.add_argument("-t", "--targets", type=lambda x: isValidFile(x, parser),
                    help="A BED file listing targets of interest")


def main(args=None, sysStdin=None, printPrefix="PRODUSE-COLLAPSE"):

    if args is None:
        if sysStdin is None:
            args = parser.parse_args()
        else:
            args = parser.parse_args(sysStdin)
        args = vars(args)

    # If a config file was specified, parse arguments from it
    if args["config"] is not None:
        config = ConfigObj(args["config"])
        try:
            for argument, parameter in config["collapse"].items():
                if argument in args and args[argument] is None:  # Aka this argument is used by collapse, and a parameter was not provided at the command line
                    args[argument] = parameter
        except KeyError:  # i.e. there is no section named "collapse" in the config file
            sys.stderr.write(
                "ERROR: Unable to locate a section named \'collapse\' in the config file \'%s\'\n" % (args["config"]))
            exit(1)

    # Re-parse and to validate the inputs
    # This is done here, as we can't check ahead of time what is in the config file
    args = validateArgs(args)

    # If the input or output is a pipe, prepare to open that
    if args["input"] == "-":
        args["input"] = sys.stdin
    if args["output"] == "-":
        args["output"] = sys.stdout

    # Sanity check some paramters
    if args["family_mismatch"] < 0:
        raise parser.error("\'-fmm/--family_mismatch\' threshold must be greater or equal to 0")
    elif args["duplex_mismatch"] < 0:
        raise parser.error("\'-dmm/--duplex_mismatch\' threshold must be greater or equal to 0")
    elif len(args["family_mask"]) != len(args["duplex_mask"]):
        raise parser.error(
            "The lengths of \'-fm\--family_mask\' and \'-dm\--duplex_mask\' must be the same, because the barcode size is the same!")

    # Convert the user-specified barcode sequences into a list of barcode indices
    # Double these, since the barcode from both the forward and reverse read will be considered
    familyMask = args["family_mask"] * 2
    duplexMask = args["duplex_mask"] * 2
    familyIndices = list(i for i in range(0, len(familyMask)) if familyMask[i] == "1")
    duplexIndices = list(i for i in range(0, len(duplexMask)) if duplexMask[i] == "1")
    inBAM = pysam.AlignmentFile(args["input"], "rb")

    # As of pysam V0.14.0, the header is now managed using an AlignmentHeader class.
    # Thus, support both approaches
    if version.parse(pysam.__version__) >= version.parse("0.14.0"):
        raise NotImplementedError("Pysam version 0.14.0 and above are currently not supported.")

    # Add this command (collapse) to the BAM header
    header = inBAM.header
    if "PG" not in header:
        header["PG"] = []

    # Format the input command in such a way that it can be added to the header (i.e. follow BAM specifications)
    command = "produse collapse"
    for argument, parameter in args.items():
        if not parameter:
            continue
        command += " --" + argument
        if not isinstance(parameter, bool):
            command += " " + str(parameter)
    header["PG"].append({"ID": "PRODUSE-COLLAPSE", "PN": "ProDuSe", "CL": command})

    outBAM = pysam.AlignmentFile(args["output"], "wb", template=inBAM, header=header)

    readProcessor = FamilyCoordinator(inBAM, args["reference"], familyIndices, args["family_mismatch"],
                                      duplexIndices, args["duplex_mismatch"], len(args["duplex_mask"]) * 2,
                                      args["targets"], args["tag_family_members"], printPrefix=printPrefix)

    for read in readProcessor:
        outBAM.write(read)


if __name__ == "__main__":
    main()
