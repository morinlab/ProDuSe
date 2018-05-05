#! /usr/bin/env python

import argparse
import os
import sys
import pysam
from configobj import ConfigObj
from packaging import version
import time


class ReadIterator:
    """
    Identifies read pairs which overlap, and merges the overlapping portion into a consensus
    The consensus overlap is then assigned to a single read only (thus clipping overlap)
    """

    def __init__(self, inFile, tag, printPrefix="PRODUSE-CLIPOVERLAP"):
        self.inFile = inFile
        self._waitingForMate = {}
        self._trimR1 = True
        self._pairCount = 0
        self._generateTag = tag
        self._printPrefix = printPrefix

    def _possibleOverlap(self, read1, read2):
        """
        Determine if there is a reasonable chance of overlap between the two reads

        Note that, by default, if the read pair maps to unique chromosomes, or are very far apart,
        they will be considered non-overlapping. However, in very rare cases, reads map still overlap
        (for example, large deletions or translocations)
        I would advise NOT considering these as overlapping, as clipping overlapping bases of these
        reads may result in the inability to identify translocation breakpoints or larger deletions
        downstream

        :param read1: A pysam.AlignedSegment()
        :param read2: A pysam.AlignedSegment()

        :returns: A boolean indicating if this read pair could possibly overlap
        """
        # If one of the reads is unmapped, they cannot overlap
        if read1.is_unmapped or read2.is_unmapped:
            return False

        # Next, check to see if the reads are close enough to possibly be overlapping
        if read1.reference_end >= read2.reference_start:
            return True
        else:
            return False

    def cigarToList(self, cigarTuples,  posLocation=None):
        """
        Convert a pysam-style cigar sequence (A list of tuples) into a list

        In addition, outputs the size of any leading soft-clipping
        In addition, if an reference position offset is provided, calculates the location in the cigar sequence
        which coresponds to that reference position
        """

        # If there is no cigar sequence, return nothing
        if not cigarTuples:
            return (), 0
        cigarList = []
        firstElement = True
        leadingSoftClipping = 0

        cigarOffset = 0
        seqOffset = 0
        i = 0
        for cigarElement in cigarTuples:
            cigOp = cigarElement[0]
            # Flag leading soft clipping
            if firstElement and cigOp == 4:
                leadingSoftClipping = cigarElement[1]
            # If the cigar operator at a given reference position is desired, calculate that here
            if posLocation and i < posLocation:
                if cigOp == 1:  # Insertions do not consume a reference position, so account for this
                    cigarOffset += cigarElement[1]
                elif cigOp == 2:  # Deletions do not consume sequence or quality positions, so account for this
                    if i + cigarElement[1] < posLocation:  # Account for the rare case where read2 starts halfway through a deletion
                        seqOffset -= cigarElement[1]
                        i += cigarElement[1]
                    else:
                        while i < posLocation:
                            i += 1
                            seqOffset -= 1
                elif cigOp != 4:
                    i += cigarElement[1]

            cigarList.extend([cigOp] * cigarElement[1])
            firstElement = False

        seqOffset = cigarOffset + seqOffset
        return cigarList, leadingSoftClipping, cigarOffset, seqOffset

    def listToCigar(self, cigar):
        """
        Converts a list of cigar operators into a pysam-style cigar sequence
        """
        # If no cigar string was provided, just return nothing
        if len(cigar) == 0:
            return ()

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

        return cigarTuples

    def _checkIndels(self, read1Ins, read1Type, read2Ins, read2Type, consensusSeq, consensusQual, consensusCigar, assignedBases):
        # If only one read has an indel, we can safely assume it is an artifact, and ignore them, since a simpler
        # alignment exists
        # TODO: Handle assigned bases in the case of indels
        read1IndelNum = len(read1Ins)
        read2IndelNum = len(read2Ins)
        if read1IndelNum and read2IndelNum:
            # For simplicities' sake, use the read with the lower number of events, since that is more likely to
            # be reality
            if read1IndelNum < read2IndelNum:
                events = read1Ins
                type = read1Type
            elif read1IndelNum > read2IndelNum:
                events = read2Ins
                type = read2Type
            else:  # In the case of a tie, use the read with shortest number of bases affected
                read1IndelLength = sum(len(x[2]) for x in read1Ins)
                read2IndelLength = sum(len(x[2]) for x in read2Ins)
                if read1IndelLength < read2IndelLength:
                    events = read1Ins
                    type = read1Type
                else:
                    events = read2Ins
                    type = read2Type

            # Add these INDELs back into the original read
            cOffset = 0  # To account for new indels
            sOffset = 0
            for event in events:
                cigarPos = event[0] + cOffset
                seqPos = event[1] + sOffset
                seq = event[2]
                qual = event[3]

                eventLength = len(seq)
                if qual is None:  # i.e. this event is a deletion
                    for i in range(cigarPos, cigarPos + eventLength):
                        consensusCigar[i] = 2
                        del consensusSeq[seqPos]
                        del consensusQual[seqPos]
                        del assignedBases[seqPos]
                    sOffset -= eventLength
                else:  # i.e. this event is an insertion
                    # Add the cigar
                    tmp = consensusCigar[:cigarPos]
                    tmp.extend([1] * eventLength)
                    tmp.extend(consensusCigar[cigarPos:])
                    consensusCigar = tmp

                    # Add the sequence
                    tmp = consensusSeq[:seqPos]
                    tmp.extend(seq)
                    tmp.extend(consensusSeq[seqPos:])
                    consensusSeq = tmp

                    # Add the quality scores
                    tmp = consensusQual[:seqPos]
                    tmp.extend(qual)
                    tmp.extend(consensusQual[seqPos:])
                    consensusQual = tmp

                    # Add the event type
                    tmp = assignedBases[:seqPos]
                    tmp.extend([type] * eventLength)
                    tmp.extend(assignedBases[seqPos:])
                    assignedBases = tmp

                    sOffset += eventLength
                    cOffset += eventLength

        return consensusSeq, consensusQual, consensusCigar, assignedBases

    def identifyAndClipOverlap(self, read1, read2):
        """
        Identify bases which overlap, the two reads, and generate a consensus from them
        """

        # Assuming read 1 starts first, identify where the next read starts, relative to the reference
        r2StartIndex = read2.reference_start - read1.reference_start
        consensusLength = read1.reference_end - read2.reference_start  # Use this to improve runtime by pre-computing the overlap length

        # If read2 is completely ecliped by read1, the consensus length will be the length of read 2
        if consensusLength > read2.reference_length:
            consensusLength = read2.reference_length

        # Convert the cigar tuples into lists of cigar elements
        read2Cigar, r2SoftClip, read2CigarOffset, read2SeqOffset = self.cigarToList(read2.cigartuples)
        read1Cigar, r1SoftClip, read1CigarOffset, read1SeqOffset = self.cigarToList(read1.cigartuples, r2StartIndex)

        r1CigarIndex = r2StartIndex + r1SoftClip + read1CigarOffset
        r2CigarIndex = r2SoftClip + read2CigarOffset

        # Iterate over the sequences and qualities scores, and generate a consensus
        # As python stores lists as arrays, pre-construct each list to the consensus length for drastically improved performance
        iSeq = 0
        iCigar = 0
        consensusSeq = [None] * consensusLength
        consensusQual = [None] * consensusLength
        consensusCigar = [None] * consensusLength
        assignedBases = [None] * consensusLength  # Indicate which bases were used to generate a consensus  (S=Both, F=Forward read, R=Reverse read)

        read1Ins = []  # The format of these lists are going to be ([Cigar] 1=Insertion, 2=Deletion), location, sequence, quality)
        read1Index = r1SoftClip + r2StartIndex + read1SeqOffset
        read1SeqClipPoint = read1Index
        read1CigarClipPoint = r2StartIndex + read1CigarOffset + r1SoftClip
        read2Ins = []
        read2Index = r2SoftClip + read2SeqOffset
        read1Type = "R" if read1.is_reverse else "F"
        read2Type = "R" if read2.is_reverse else "F"
        read2StartOffset = 0

        # Store these as tuples to improve indexing performance
        read1_query_sequence = tuple(read1.query_sequence)
        read1_query_qualities = tuple(read1.query_qualities)
        read2_query_sequence = tuple(read2.query_sequence)
        read2_query_qualities = tuple(read2.query_qualities)

        try:
            while True:
                # Obtain the base, quality, and cigar operator at this position
                b1 = read1_query_sequence[read1Index]
                q1 = read1_query_qualities[read1Index]
                c1 = read1Cigar[r1CigarIndex]
                b2 = read2_query_sequence[read2Index]
                q2 = read2_query_qualities[read2Index]
                c2 = read2Cigar[r2CigarIndex]

                # Easiest case, and most common. Both bases are mapped to the reference
                if c1 == 0 and c2 == 0:
                    # If both reads agree, this is easy. Use this base as the consensus
                    if b1 == b2:
                        consensusCigar[iCigar] = c1
                        consensusSeq[iSeq] = b1
                        # Avoid using python's max function. This should be faster
                        if q1 > q2:
                            consensusQual[iSeq] = q1
                        else:
                            consensusQual[iSeq] = q2
                        assignedBases[iSeq] = "S"

                    # If these bases disagree, use the base with the highest quality
                    elif q1 > q2:
                        consensusCigar[iCigar] = c1
                        consensusSeq[iSeq] = b1
                        consensusQual[iSeq] = q1
                        assignedBases[iSeq] = read1Type
                    else:  # If the bases disagree and the quality score is the same, just use read 2's base
                        consensusCigar[iCigar] = c2
                        consensusSeq[iSeq] = b2
                        consensusQual[iSeq] = q2
                        assignedBases[iSeq] = read2Type
                    r1CigarIndex += 1
                    r2CigarIndex += 1
                    read1Index += 1
                    read2Index += 1
                    iSeq += 1
                    iCigar += 1
                    if c2 != 1:
                        read2StartOffset += 1

                # If we have reached the soft-clipped segment (if applicable), stop processing
                elif c1 == 4 or c2 == 4:
                    break

                # Insertion
                elif c1 == 1 and c2 != 1:
                    # If the mate is soft-clipped, use this insertion as the consensus
                    read1Ins.append((iCigar, iSeq, [], []))

                    while c1 == 1:
                        read1Ins[-1][2].append(b1)
                        read1Ins[-1][3].append(q1)
                        read1Index += 1
                        r1CigarIndex += 1
                        b1 = read1_query_sequence[read1Index]
                        q1 = read1_query_qualities[read1Index]
                        c1 = read1Cigar[r1CigarIndex]

                elif c2 == 1 and c1 != 1:
                    read2Ins.append((iCigar, iSeq, [], []))

                    while c2 == 1:
                        read2Ins[-1][2].append(b2)
                        read2Ins[-1][3].append(q2)
                        read2Index += 1
                        r2CigarIndex += 1
                        b2 = read2_query_sequence[read2Index]
                        q2 = read2_query_qualities[read2Index]
                        c2 = read2Cigar[r2CigarIndex]

                # Handle deletions
                elif c1 == 2 and c2 == 2:
                    consensusCigar[iCigar] = 2
                    del consensusSeq[iSeq]
                    del consensusQual[iSeq]
                    del assignedBases[iSeq]
                    r1CigarIndex += 1
                    r2CigarIndex += 1
                    read2StartOffset += 1
                    iCigar += 1
                elif c1 == 2:
                    read1Ins.append((iCigar, iSeq, [], None))
                    while c1 == 2 and c2 != 2 and c2 != 4:
                        consensusSeq[iSeq] = b2
                        consensusQual[iSeq] = q2
                        consensusCigar[iCigar] = c2
                        assignedBases[iSeq] = read2Type
                        read1Ins[-1][2].append("N")
                        read2Index += 1
                        r1CigarIndex += 1
                        r2CigarIndex += 1
                        read2StartOffset += 1
                        iSeq += 1
                        iCigar += 1
                        c2 = read2Cigar[r2CigarIndex]
                        b2 = read2_query_sequence[read2Index]
                        q2 = read2_query_qualities[read2Index]
                        c1 = read1Cigar[r1CigarIndex]

                elif c2 == 2:

                    read2Ins.append((iCigar, iSeq, [], None))

                    while c2 == 2 and c1 != 2 and c1 != 4:
                        consensusSeq[iSeq] = b1
                        consensusQual[iSeq] = q1
                        consensusCigar[iCigar] = c1
                        assignedBases[iSeq] = read1Type
                        read2Ins[-1][2].append("N")
                        read1Index += 1
                        r1CigarIndex += 1
                        r2CigarIndex += 1
                        read2StartOffset += 1
                        iSeq += 1
                        iCigar += 1
                        c2 = read2Cigar[r2CigarIndex]
                        b1 = read1_query_sequence[read1Index]
                        q1 = read1_query_qualities[read1Index]
                        c1 = read1Cigar[r1CigarIndex]

                else:  # Both reads either contain an insertion
                    # If both reads agree, this is easy
                    if b1 == b2:
                        consensusCigar.insert(iCigar, c1)
                        consensusSeq.insert(iSeq, b1)
                        # Avoid using python's max function. This should be faster
                        if q1 > q2:
                            consensusQual.insert(iSeq, q1)
                        else:
                            consensusQual.insert(iSeq, q2)
                        assignedBases.insert(iSeq, "S")

                    # If these bases disagree, use the base with the highest quality
                    elif q1 > q2:
                        consensusCigar.insert(iCigar, c1)
                        consensusSeq.insert(iSeq, b1)
                        consensusQual.insert(iSeq, q1)
                        assignedBases.insert(iSeq, read1Type)
                    else:  # If the bases disagree and the quality score is the same, just use read 2's base, I guess
                        consensusCigar.insert(iCigar, c2)
                        consensusSeq.insert(iSeq, b2)
                        consensusQual.insert(iSeq, q2)
                        assignedBases.insert(iSeq, read2Type)
                    r1CigarIndex += 1
                    r2CigarIndex += 1
                    read1Index += 1
                    read2Index += 1
                    iSeq += 1
                    iCigar += 1
                    if c2 != 1:
                        read2StartOffset += 1

        except IndexError:
            pass


        # Now that an overlap consensus has been generated, lets look at any indels which we encountered
        consensusSeq, consensusQual, consensusCigar, assignedBases = self._checkIndels(read1Ins, read2Type,
                                                                        read2Ins, read2Type, consensusSeq,
                                                                        consensusQual, consensusCigar, assignedBases)


        # Assign the consensus generated to one of two reads
        if consensusSeq:

            try:
                consensusSeq = "".join(consensusSeq)
            except TypeError as e:
                raise TypeError("An error was encountered while processing read pair \'%s\'. This usually indicates that this read pair has an invalid CIGAR sequence" % read1.query_name) from e

            # Trim read 1
            read1Qual = read1.query_qualities[:read1SeqClipPoint]
            read1.query_sequence = read1.query_sequence[:read1SeqClipPoint]
            read1Cigar = read1Cigar[:read1CigarClipPoint]

            # Trim read 2
            read2Qual = read2.query_qualities[read2Index:]
            read2.query_sequence = read2.query_sequence[read2Index:]
            read2Cigar = read2Cigar[r2CigarIndex:]

            # Add the consensus to the specified read
            # We alternate the sequence that is clipped to minimize the likelihood of strand bias filtering
            # will remove a variant if it was actually supported by both reads (the overlap)

            # To improve compatibility with tools that require a non-soft clipped sequence in both reads, add
            # one base from the consensus to the read that will be clipped (here, is is read 1)
            # But make sure this doesn't cause the clipping point to impact an INDEL
            try:
                if self._trimR1:
                    conCigarClipPoint = 0
                    conSeqClipPoint = 0
                    clipOffset = 0
                    while True:
                        clippedCigar = consensusCigar[conCigarClipPoint]
                        if clippedCigar == 0:
                            conSeqClipPoint += 1
                            conCigarClipPoint += 1
                            clipOffset += 1
                            if consensusCigar[conCigarClipPoint] == 0 or consensusCigar[conCigarClipPoint] == 4:  # The next position is also a match. Split between these two
                                break
                            else:
                                continue
                        elif clippedCigar == 2:  # Deletion, no sequence consumed
                            conCigarClipPoint += 1
                            clipOffset += 1
                        elif clippedCigar == 1:  # Insertion, start position not affected
                            conCigarClipPoint += 1
                            conSeqClipPoint += 1
                        else:
                            raise TypeError(
                                "An error occured while splitting the consensus sequence. Cigar operator encountered: %s" % clippedCigar)
                else:
                    conCigarClipPoint = -1
                    conSeqClipPoint = -1
                    clipOffset = -1
                    while True:
                        clippedCigar = consensusCigar[conCigarClipPoint]
                        if clippedCigar == 0:  # This position is a normal mapping
                            conCigarClipPoint -= 1
                            conSeqClipPoint -= 1
                            clipOffset -= 1
                            if consensusCigar[conCigarClipPoint] == 0 or consensusCigar[conCigarClipPoint] == 4:  # i.e. the next position is also a match. We should split between these two
                                break
                            else:
                                continue
                        elif clippedCigar == 2:  # Deletion, no sequence consumed
                            conCigarClipPoint -= 1
                            clipOffset -= 1
                        elif clippedCigar == 1:  # Insertion, start position not affected
                            conCigarClipPoint -= 1
                            conSeqClipPoint -= 1
                        else:
                            raise TypeError(
                                "An error occured while splitting the consensus sequence. Cigar operator encountered: %s" % clippedCigar)
            except IndexError:  # Only a single base overlaps
                # Don't do anything special with the overlap
                # This shouldn't cause any problems. If it does, then for some reason the BAM file has a read that is only a single
                # nucleotide long, and if that's the case, there are FAR bigger problems than the possible incompatabilities
                # that this will cause
                pass

            # Add the consensus to read 1
            # Sequence
            try:
                read1.query_sequence = read1.query_sequence + consensusSeq[:conSeqClipPoint]
            except TypeError:  # Because pysam converts "" into None
                read1.query_sequence = consensusSeq[:conSeqClipPoint]
            # Quality score
            read1Qual.extend(consensusQual[:conSeqClipPoint])
            read1.query_qualities = read1Qual
            # Cigar string
            read1Cigar.extend(consensusCigar[:conCigarClipPoint])

            # Add the consensus to read 2
            # Sequence
            try:
                read2.query_sequence = consensusSeq[conSeqClipPoint:] + read2.query_sequence
            except TypeError:  # Because pysam converts "" into None
                read2.query_sequence = consensusSeq[conSeqClipPoint:]
            # Quality
            consensusQual = consensusQual[conSeqClipPoint:]
            consensusQual.extend(read2Qual)
            read2.query_qualities = consensusQual
            # Cigar
            consensusCigar = consensusCigar[conCigarClipPoint:]
            consensusCigar.extend(read2Cigar)
            read2Cigar = consensusCigar

            # Start position of read2
            read2.reference_start = read2.reference_start + clipOffset
            if not self._trimR1:
                read2.reference_start += read2StartOffset
            read1.next_reference_start = read2.reference_start

            if self._generateTag:  # Add a tag indicating from which read a base originated from
                try:
                    tag = "".join(assignedBases[conSeqClipPoint:]) + read2Type * len(read2.query_sequence)
                except TypeError:
                    tag = "".join(assignedBases)
                read2.set_tag("co", tag)

                try:
                    tag = read1Type * len(read1.query_sequence) + "".join(assignedBases[:conSeqClipPoint])
                except TypeError:
                    tag = "".join(assignedBases[:conSeqClipPoint])
                read1.set_tag("co", tag)

            read1.cigartuples = self.listToCigar(read1Cigar)
            read2.cigartuples = self.listToCigar(read2Cigar)

            # Finally, update the mate cigar tags with the new cigar sequences
            read1.set_tag("MC", read2.cigarstring)
            read2.set_tag("MC", read1.cigarstring)

            self._trimR1 = not self._trimR1
            """
            # Add a tag indicating from which read a base was generated from
            if self._generateTag:
                try:
                    read1Tag = read1Type * len(read1.query_sequence) + "".join(assignedBases)
                except TypeError:
                    read1Tag = "".join(assignedBases)
                try:
                    read2Tag = "".join(assignedBases) + read2Type * len(read2.query_sequence)
                except TypeError:
                    read2Tag = "".join(assignedBases)
                read1.set_tag("co", read1Tag)
                read2.set_tag("co", read2Tag)

            # Add the consensus sequence to each read
            try:
                read1.query_sequence = read1.query_sequence + consensusSeq
            except TypeError:  # Because pysam converts "" into None
                read1.query_sequence = consensusSeq

            try:
                read2.query_sequence = consensusSeq + read2.query_sequence
            except TypeError:  # Because pysam converts "" into None
                read2.query_sequence = consensusSeq

            # Add the qualities
            read1Qual.extend(consensusQual)
            read1.query_qualities = read1Qual
            consensusQual.extend(read2Qual)

            read1.query_qualities = read1Qual
            read2.query_qualities = consensusQual  # Actually read2Qual

            # Modify the cigar string to account for the soft-clipping
            if self._trimR1:
                read1Cigar.extend([4]*len(consensusSeq))

                # To avoid problems with downstream tools that require a non-soft clipped sequence, add one base of the
                # consensus to the other read
                
                consensusCigar.extend(read2Cigar)
                read2Cigar = consensusCigar
            else:
                read1Cigar.extend(consensusCigar)
                tmp = [4] * len(consensusSeq)
                tmp.extend(read2Cigar)
                read2Cigar = tmp
                read2.reference_start = read2.reference_start + read2StartOffset
                read1.next_reference_start = read2.reference_start

            read1.cigartuples = self.listToCigar(read1Cigar)
            read2.cigartuples = self.listToCigar(read2Cigar)

            self._trimR1 = not self._trimR1
            """
    def next(self):
        return self.__iter__()

    def __iter__(self):

        sys.stderr.write("\t".join([self._printPrefix, time.strftime("%X"), "Starting...\n"]))

        try:
            while True:
                read = next(self.inFile)

                # Don't try to clip supplementary or secondary alignments
                if read.is_supplementary or read.is_secondary:
                    yield read
                    continue

                # Ignore reads that map to different chromosomes, as they will never overlap
                try:
                    if read.reference_name != read.next_reference_name:
                        yield read
                        continue
                except ValueError:  # One of the reads is unmapped. They cannot overlap
                    yield read
                    continue

                # Have we encountered this read's mate before?
                if read.query_name not in self._waitingForMate:
                    # If not, buffer this read until we encounter its' mate
                    self._waitingForMate[read.query_name] = read
                    continue

                self._pairCount += 1
                if self._pairCount % 100000 == 0:
                    sys.stderr.write("\t".join([self._printPrefix, time.strftime("%X"), "Pairs Processed:%s\n" % self._pairCount]))
                # Once a read pair has been obtained, we need to generate an alignment between the read and it's mate,
                # and identify any overlapping bases

                # Assign read 1 and read 2 based upon start position
                if read.reference_start < self._waitingForMate[read.query_name].reference_start:
                    read1 = read
                    read2 = self._waitingForMate[read.query_name]
                else:  # We don't actually care if they start at the same position, since they will be considered overlapping anyways
                    read1 = self._waitingForMate[read.query_name]
                    read2 = read

                # Delete the buffered read to reduce the memory footprint
                del self._waitingForMate[read.query_name]

                # Perform clipoverlap
                if self._possibleOverlap(read1, read2):
                    self.identifyAndClipOverlap(read1, read2)
                yield read1
                yield read2


        except StopIteration:

            # Output all remaining unpaired reads
            unpairedReadNum = len(self._waitingForMate)
            for read in self._waitingForMate.values():
                yield read

            if self._pairCount % 100000 != 0:
                sys.stderr.write(
                    "\t".join([self._printPrefix, time.strftime("%X"), "Pairs Processed:%s\n" % self._pairCount]))
            if unpairedReadNum > 0:
                sys.stderr.write(
                    "\t".join([self._printPrefix, time.strftime("%X"), "Unable to find a mate for %s reads\n" % unpairedReadNum]))


def reparseArgs(args):
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
        # If the parameter is a list, we need to add each element seperately
        if isinstance(parameter, list):
            for p in parameter:
                listArgs.append(str(p))
        else:
            listArgs.append(str(parameter))

    parser = argparse.ArgumentParser(
        description="Generates a consensus for overlapping bases in a read pair, and assigns the consensus to a single read")
    parser.add_argument("-c", "--config", metavar="INI", type=lambda x: isValidFile(x, parser),
                        help="An optional configuration file, which can specify any parameters")
    parser.add_argument("-i", "--input", required=True, metavar="BAM/SAM", type=lambda x: isValidFile(x, parser, True),
                        help="Input BAM file. Does not need to be sorted (use \"-\" to read from stdin [and \"Control + D\" to stop reading])")
    parser.add_argument("-o", "--output", required=True, metavar="BAM/SAM",
                        help="A path to an output BAM file. Will be UNSORTED (use \"-\" for stdout)")
    parser.add_argument("--tag_origin", action="store_true",
                        help="Add a read tag indicating from which read a consensus base originated")
    args = parser.parse_args(listArgs)
    return vars(args)


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


parser = argparse.ArgumentParser(
    description="Generates a consensus for overlapping bases between a read pair, and assigns the consensus to a single read")
parser.add_argument("-c", "--config", metavar="INI", type=lambda x: isValidFile(x, parser),
                    help="An optional configuration file, which can specify any parameters")
parser.add_argument("-i", "--input", metavar="BAM/SAM", type=lambda x: isValidFile(x, parser, True),
                    help="Input BAM file. Does not need to be sorted (use \"-\" to read from stdin [and \"Control + D\" to stop reading])")
parser.add_argument("-o", "--output", metavar="BAM/SAM", help="A path to an output BAM file. Will be UNSORTED (use \"-\" for stdout)")
parser.add_argument("--tag_origin", action="store_true", help="Add a read tag indicating from which read a consensus base originated")


def main(args=None, sysStdin=None, printPrefix="PRODUSE-CLIPOVERLAP"):
    if args is None:
        if sysStdin is None:
            args = parser.parse_args()
        else:
            args = parser.parse_args(sysStdin)
        args = vars(args)

    # If a config file was specified, parse the arguments out of that
    if args["config"]:
        config = ConfigObj(args["config"])
        try:
            for argument, parameter in config["clipoverlap"].items():
                if argument in args and not args[argument]:  # Aka this argument is used by clipovlerap, and a parameter was not provided at the command line
                    args[argument] = parameter
        except KeyError:  # i.e. there is no section named "clipoverlap" in the config file
            sys.stderr.write("ERROR: Unable to locate a section named \'clipoverlap\' in the config file \'%s\'\n" % (
            args["config"]))
            exit(1)

    args = reparseArgs(args)

    # Open the input and output BAM files for reading
    inBAM = pysam.AlignmentFile(args["input"])  # Pysam claims to auto-detect the input format

    # Add this command (clipoverlap) to the BAM header
    # As of pysam V0.14.0, the header is now managed using an AlignmentHeader class.
    # Thus, support both approaches (when the API for AlignmentHeaders is listed!!!)
    if version.parse(pysam.__version__) >= version.parse("0.14.0"):
        raise NotImplementedError("Pysam version 0.14.0 and above are currently not supported.")

    header = inBAM.header
    if "PG" not in header:
        header["PG"] = []

    # Format the input command in such a way that it can be added to the header (i.e. follow BAM specifications)
    command = "produse clipoverlap"
    for argument, parameter in args.items():
        if not parameter:
            continue
        command += " --" + str(argument)
        if not isinstance(parameter, bool):
            command += " " + str(parameter)
    header["PG"].append({"ID": "PRODUSE-CLIPOVERLAP", "PN": "ProDuSe", "CL": command})

    # If streams were specified as input or output (represented by a pipe symbol), set those
    if args["input"] == "-":
        args["input"] = sys.stdin
        sys.stderr.write("Reading from standard input. Use Control + D to cancel\n")

    if args["output"] == "-":
        args["output"] = sys.stdout
        # Determine the output file type by the file type of the input File
        if args["input"] is sys.stdin:  # Default to BAM output
            sys.stderr.write("Unable to determine output file type from input stream. Defaulting to BAM")
            outType ="bam"
        else:
            outType = args["input"].split(".")[-1]
    else:
        # Determine the output file type via the file extension
        outType = args["output"].split(".")[-1]

    if outType == "bam":
        outBAM = pysam.AlignmentFile(args["output"], mode="wb", header=header)
    elif outType == "cram":
        outBAM = pysam.AlignmentFile(args["output"], mode="wc", header=header)
    else:
        outBAM = pysam.AlignmentFile(args["output"], mode="w", header=header)

    processor = ReadIterator(inBAM, args["tag_origin"], printPrefix)

    for read in processor:
        outBAM.write(read)


if __name__ == "__main__":
    main()
