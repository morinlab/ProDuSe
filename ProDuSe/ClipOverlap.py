#! /usr/bin/env python

import argparse
import os
import sys
import pysam
from configobj import ConfigObj
# from skbio import alignment
import time


class ReadIterator:
    """
    Identifies read pairs which overlap, and merges the overlapping portion into a consensus
    The consensus overlap is then assigned to a single read only (thus clipping overlap)
    """

    def __init__(self, inFile, noCleanup, noTag):
        self.inFile = inFile
        self._waitingForMate = {}
        self._trimR1 = True
        self._pairCount = 0
        self._cleanup = not noCleanup
        self._generateTag = not noTag

    def _possibleOverlap(self, read1, read2, maxGap=200, allowSeperateReference=False):
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
        :param maxGap: An int representing the maximum number of bases allowed between the end of the first read and the start
        of the next read before the read pair will be considered as not overlapping
        :param allowSeperateReference: Consider reads which map to different chromosomes as possibly overlapping

        :returns: A boolean indicating if this read pair could possibly overlap
        """
        # If one of the reads is unmapped, they cannot overlap
        if read1.is_unmapped or read2.is_unmapped:
            return False

        # Check to see if the reads map to the same reference
        if read1.reference_name != read2.reference_name:
            return allowSeperateReference

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
        # If only one read has an indel, we can safely assume it is an artifact, and ignore them, since a simplier
        # alignment exists
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
                        assignedBases[i] = type
                        del consensusSeq[seqPos]
                        del consensusQual[seqPos]
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
                    tmp = assignedBases[:cigarPos]
                    tmp.extend([type] * eventLength)
                    tmp.extend(assignedBases[cigarPos:])
                    assignedBases = tmp

                    sOffset += eventLength
                    cOffset += eventLength

        return consensusSeq, consensusQual, consensusCigar, assignedBases

    def identifyOverlap(self, read1, read2):
        """
        Identify bases which overlap, the two reads, and generate a consensus from them
        """

        # In the simplist case, these overlapping bases should line up properly in terms
        # of the reference
        # Lets obtain those bases

        # Assuming read 1 starts first, identify where the next read starts, relative to the reference
        r2StartIndex = read2.reference_start - read1.reference_start
        consensusLength = read1.reference_end - read2.reference_start

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
        assignedBases = [None] * consensusLength # Indicate which bases were used to generate a consensus  (S=Both, F=Forward read, R=Reverse read)

        read1Ins = []  # The format of these lists are going to be ([Cigar] 1=Insertion, 2=Deletion), location, sequence, quality)
        read1Index = r1SoftClip + r2StartIndex + read1SeqOffset
        read1SeqClipPoint = read1Index
        read1CigarClipPoint = r2StartIndex + read1CigarOffset + r1SoftClip
        read2Ins = []
        read2Index = r2SoftClip + read2SeqOffset
        read1Type = "1" if read1.is_read1 else "2"
        read2Type = "1" if read2.is_read1 else "2"
        read2StartOffset = 0

        # Store these as tuples to reduce the number of  calls, and thus improve efficiency
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
                    # If both reads agree, this is easy
                    if b1 == b2:
                        consensusCigar[iCigar] = c1
                        consensusSeq[iSeq] = b1
                        # Avoid using python's max function. This should be faster
                        if q1 > q2:
                            consensusQual[iSeq] = q1
                        else:
                            consensusQual[iSeq] = q2
                        assignedBases[iCigar] = "S"

                    # If these bases disagree, use the base with the highest quality
                    elif q1 > q2:
                        consensusCigar[iCigar] = c1
                        consensusSeq[iSeq] = b1
                        consensusQual[iSeq] = q1
                        assignedBases[iCigar] = read1Type
                    else:  # If the bases disagree and the quality score is the same, just use read 2's base, I guess
                        consensusCigar[iCigar] = c2
                        consensusSeq[iSeq] = b2
                        consensusQual[iSeq] = q2
                        assignedBases[iCigar] = read2Type
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
                    assignedBases[iCigar] = "S"
                    del consensusSeq[iSeq]
                    del consensusQual[iSeq]
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
                        assignedBases[iCigar] = read2Type
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
                        assignedBases[iCigar] = read1Type
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
                        assignedBases.insert(iCigar, "S")

                    # If these bases disagree, use the base with the highest quality
                    elif q1 > q2:
                        consensusCigar.insert(iCigar, c1)
                        consensusSeq.insert(iSeq, b1)
                        consensusQual.insert(iSeq, q1)
                        assignedBases.insert(iCigar, read1Type)
                    else:  # If the bases disagree and the quality score is the same, just use read 2's base, I guess
                        consensusCigar.insert(iCigar, c2)
                        consensusSeq.insert(iSeq, b2)
                        consensusQual.insert(iSeq, q2)
                        assignedBases.insert(iCigar, read2Type)
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

        if consensusSeq:  # Add the overlap consensus to one of the two reads

            # To ensure a valid record is produced for downstream tools, check to see if the clip point of either
            # read coresponds to an INDEL. if it does, we need to add the consensus to that sequence, to ensure a valid
            # record is generated
            if read1CigarClipPoint == 0:
                read1ClipPoint = None
            else:
                read1ClipPoint = read1Cigar[read1CigarClipPoint]
            try:
                read2ClipPoint = read2Cigar[r2CigarIndex]
            except IndexError:
                read2ClipPoint = None
            if read1ClipPoint == 1 or read1ClipPoint == 2:
                # If the clip point of both reads is an INDEL, don;t clip them....
                if read2ClipPoint == 1 or read2ClipPoint == 2:
                    return  # This could be handled better, but there are only so many edge cases I can realistically deal with
                self._trimR1 = False
            elif read2ClipPoint == 1 or read2ClipPoint == 2:
                self._trimR1 = True

            consensusSeq = "".join(consensusSeq)

            # Trim read 1
            read1Qual = read1.query_qualities[:read1SeqClipPoint]
            read1.query_sequence = read1.query_sequence[:read1SeqClipPoint]
            read1Cigar = read1Cigar[:read1CigarClipPoint]

            # Trim read 2
            read2Qual = read2.query_qualities[read2Index:]
            read2.query_sequence = read2.query_sequence[read2Index:]
            read2Cigar = read2Cigar[r2CigarIndex:]

            if self._trimR1:

                if read2.query_sequence:  # Because pysam converts "" into None
                    read2.query_sequence = consensusSeq + read2.query_sequence
                else:
                    read2.query_sequence = consensusSeq

                consensusQual.extend(read2Qual)
                read2Qual = consensusQual

                consensusCigar.extend(read2Cigar)
                read2Cigar = consensusCigar

                # If we are cleaning up the output (i.e. discarding reads which are purely soft-clipped bases), do so now
                if self._cleanup and set(read1Cigar) == {4}:
                    read1.query_sequence = None
                    read1Qual = []
                    read1Cigar = []

            else:
                if read1.query_sequence:  # Because pysam converts "" into None
                    read1.query_sequence= read1.query_sequence + consensusSeq
                else:
                    read1.query_sequence = consensusSeq
                read1Qual.extend(consensusQual)
                read1Cigar.extend(consensusCigar)
                read2.reference_start = read2.reference_start + read2StartOffset
                read1.next_reference_start = read2.reference_start

                # If we are cleaning up the output (i.e. discarding reads which are purely soft-clipped bases), do so now
                if self._cleanup and set(read2Cigar) == {4}:
                    read2.query_sequence = None
                    read2Qual = []
                    read2Cigar = []

            read1.query_qualities = read1Qual
            read2.query_qualities = read2Qual
            read1.cigartuples = self.listToCigar(read1Cigar)
            read2.cigartuples = self.listToCigar(read2Cigar)

            # Add the assigned bases flag
            if self._generateTag:
                read1.set_tag("co", "".join(assignedBases))
                read2.set_tag("co", "".join(assignedBases))

            self._trimR1 = not self._trimR1

    def next(self):
        return self.__iter__()

    def __iter__(self):

        printPrefix = "PRODUSE-CLIPOVERLAP"
        sys.stderr.write("\t".join([printPrefix, time.strftime("%X"), "Starting...\n"]))

        try:
            while True:
                read = next(self.inFile)

                # Don't try to clip supplementary or secondary alignments
                if read.is_supplementary or read.is_secondary:
                    yield read
                    continue

                # Ignore reads that map to different chromosomes, as they will never overlap
                if read.reference_name != read.next_reference_name:
                    yield read
                    continue

                # Have we encountered this read's mate before?
                if read.query_name not in self._waitingForMate:
                    # If not, buffer this read until we encounter its' mate
                    self._waitingForMate[read.query_name] = read
                    continue

                self._pairCount += 1
                if self._pairCount % 100000 == 0:
                    sys.stderr.write("\t".join([printPrefix, time.strftime("%X"), "Pairs Processed: %s\n" % self._pairCount]))
                # Once a read pair has been obtained, we need to generate an alignment between the read and it's mate,
                # and identify any overlapping bases

                # Assign read 1 and read 2 based upon start position
                if read.reference_start < self._waitingForMate[read.query_name].reference_start:
                    read1 = read
                    read2 = self._waitingForMate[read.query_name]
                else:  # We don't actually care if they start at the same position, since they will be considered overlapping anyways
                    read1 = self._waitingForMate[read.query_name]
                    read2 = read

                # Delete this to reduce the memory footprint
                del self._waitingForMate[read.query_name]

                # First, identify if this read pair could possibly overlap
                if self._possibleOverlap(read1, read2):
                    self.identifyOverlap(read1, read2)
                yield read1
                yield read2


        except StopIteration:

            # Output all remaining unpaired reads
            unpairedReadNum = len(self._waitingForMate)
            for read in self._waitingForMate.values():
                yield read

            if self._pairCount % 100000 != 0:
                sys.stderr.write(
                    "\t".join([printPrefix, time.strftime("%X"), "Pairs Processed: %s\n" % self._pairCount]))
            if unpairedReadNum > 0:
                sys.stderr.write(
                    "\t".join([printPrefix, time.strftime("%X"), "Unable to find a mate for: %s reads\n" % unpairedReadNum]))


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
        if parameter == "True:
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
    parser.add_argument("-i", "--input", required=True, metavar="BAM/SAM/CRAM", type=lambda x: isValidFile(x, parser, True),
                        help="Input BAM file. Does not need to be sorted (use \"-\" to read from stdin [and \"Control + D\" to stop reading])")
    parser.add_argument("-o", "--output", required=True, metavar="BAM/SAM/CRAM",
                        help="A path to an output BAM file. Will be UNSORTED (use \"-\" for stdout)")
    parser.add_argument("--no_cleanup", action="store_true",
                        help="Output ALL non-clipped bases, even if the post-clip read is only comprised of soft-clipped bases")
    parser.add_argument("--no_tag", action="store_true",
                        help="Do not add a read tag indicating from which read a consensus base originated")
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
parser.add_argument("-i", "--input", metavar="BAM/SAM/CRAM", type=lambda x: isValidFile(x, parser, True),
                    help="Input BAM file. Does not need to be sorted (use \"-\" to read from stdin [and \"Control + D\" to stop reading])")
parser.add_argument("-o", "--output", metavar="BAM/SAM/CRAM", help="A path to an output BAM file. Will be UNSORTED (use \"-\" for stdout)")
parser.add_argument("--no_cleanup", action="store_true", help="Output ALL non-clipped bases, even if a post-clip read is only comprised of soft-clipped bases")
parser.add_argument("--no_tag", action="store_true", help="Do not add a read tag indicating from which read a consensus base originated")


def main(args=None, sysStdin=None):
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

    # If streams were specified as input or output (represented by a pipe symbol), set those
    if args["input"] == "-":
        args["input"] = sys.stdin
        sys.stderr.write("Reading from standard input. Use Control + D to cancel\n")
    if args["output"] == "-":
        args["output"] = sys.stdout

    # Open the input and output BAM files for reading
    inBAM = pysam.AlignmentFile(args["input"])  # Pysam claims to auto-detect the input format

    # Add this command (clipoverlap) to the BAM header
    header = inBAM.header
    if "PG" not in header:
        header["PG"] = []

    # Format the input command in such a way that it can be added to the header (i.e. follow BAM specifications)
    command = "produse clipoverlap"
    for argument, parameter in args.items():
        if not parameter:
            continue
        command += " --" + argument
        if not isinstance(parameter, bool):
            command += " " + parameter
    header["PG"].append({"ID": "PRODUSE-CLIPOVERLAP", "PN": "ProDuSe", "CL": command})

    # Determine the output file type via the file extension
    outType = args["output"].split(".")[-1]
    if outType == "bam":
        outBAM = pysam.AlignmentFile(args["output"], mode="wb", header=header)
    elif outType == "cram":
        outBAM = pysam.AlignmentFile(args["output"], mode="wc", header=header)
    else:
        outBAM = pysam.AlignmentFile(args["output"], mode="w", header=header)

    processor = ReadIterator(inBAM, args["no_cleanup"], args["no_tag"])

    for read in processor:
        outBAM.write(read)


if __name__ == "__main__":
    main()
