#! /usr/bin/env python

import argparse
import os
import pysam
import sys
import sortedcontainers
import pickle
from sklearn.ensemble import RandomForestClassifier
from skbio import alignment
from scipy.stats import ttest_ind
from pyfaidx import Fasta
import time
from fisher import pvalue as fisher_exact
from statistics import mean
from configobj import ConfigObj
import bisect

# Multiprocessing
import inspect
import multiprocessing

# Import version number
try:
    import __version as pVer
    import ProDuSeExceptions as pe
except ImportError:
    from ProDuSe import __version as pVer
    from ProDuSe import ProDuSeExceptions as pe


class Haplotype(object):
    """
    TODO: An object that stores an aligned sequence ("Haplotype"), it's differences relative to the reference, and
    the NGS reads which support this haplotype

    A haplotype is generated from reads which contain indels, and represents that indel
    """

    def __init__(self, seq=None, eventFromRef=None, scoreFromRef=None):
        """

        :param seq: A string containing a nucleotide sequence which coresponds to the assembled haplotype
        :param support: A list containing strings coresponding to the reads which support this event
        :param eventFromRef: A cigar string containing the distance from this alignment to the reference
        :param
        """
        self.seq = seq
        self.support = []
        self.eventFromRef=eventFromRef
        self.alignmentStructure=alignment.StripedSmithWaterman(seq, mismatch_score=-3, match_score=1)
        self.scoreFromRef = scoreFromRef


class Position(object):
    """
    Stores the allele counts and associated characteristics at this position
    """

    # Reduce memory footprint
    __slots__ = ("alleles", "qualities", "famSizes", "posParent", "mapStrand", "distToEnd", "mismatchNum",
                 "mappingQual", "readNum", "readName", "alt", "altAlleles", "depth", "altDepth",
                 "leftWindow", "ref", "rightWindow", "nearbyVar", "nWindow", "baseCounts",
                 "alleleStrongCounts", "strandCounts", "posMolec", "negMolec", "pMapStrand", "nMapStrand",
                 "alleleMapQual", "alleleBaseQual", "alleleStrandBias", "alleleVafs", "molecDepth", "alleleMismatch",
                 "alleleAvFamSize", "alleleEndDist", "duplexCounts")
    def __init__(self, refBase):
        self.ref = refBase

        # Pre-initialize these lists to a sensible size to improve performance in low-depth regions
        self.alleles = [None] * 5
        self.qualities = [None] * 5
        self.famSizes = [None] * 5
        self.posParent = [None] * 5
        self.mapStrand = [None] * 5
        self.distToEnd = [None] * 5
        self.mismatchNum = [None] * 5
        self.mappingQual = [None] * 5
        self.readNum = [None] * 5
        self.readName = [None] * 5
        self.alt = False
        self.altAlleles = {}
        self.depth = 0
        self.altDepth = 0

    def add(self, base, qual, size, posParent, mapStrand, dist, mapQual, readNum, readName):

        # Add a new base, as well as the corresponding characteristics, to this position
        try:
            self.alleles[self.depth] = base
        except IndexError:  # We have exceeded the current max buffered positions. Increase the buffer
            self.alleles.extend([None] * 150)
            self.qualities.extend([None] * 150)
            self.famSizes.extend([None] * 150)
            self.posParent.extend([None] * 150)
            self.mapStrand.extend([None] * 150)
            self.distToEnd.extend([None] * 150)
            self.mismatchNum.extend([None] * 150)
            self.mappingQual.extend([None] * 150)
            self.readNum.extend([None] * 150)
            self.readName.extend([None] * 150)
            self.alleles[self.depth] = base
        self.qualities[self.depth] = qual
        self.famSizes[self.depth] = size
        self.posParent[self.depth] = posParent
        self.mapStrand[self.depth] = mapStrand
        self.distToEnd[self.depth] = dist
        self.mappingQual[self.depth] = mapQual
        self.readNum[self.depth] = readNum
        self.readName[self.depth] = readName
        self.depth += 1

        # Does this base support a variant?
        if self.ref != "N" and self.ref != base:
            self.alt = True
            try:
                self.altAlleles[base] += 1
            except KeyError:
                self.altAlleles[base] = 1
            self.altDepth += 1
            return True
        else:
            return False

    def addMismatchNum(self, misMatchNum):
        self.mismatchNum[self.depth - 1] = misMatchNum

    def processVariant(self, leftWindow, rightWindow, nearbyVariants, nWindow):
        self.leftWindow = leftWindow
        self.rightWindow = rightWindow
        self.nearbyVar = nearbyVariants
        self.nWindow = nWindow  # How many positions were examined to find nearby variants?

    def summarizeVariant(self, minAltDepth=4, strongMoleculeThreshold=3):
        """
        Aggregate all statistics at this position
        :param strongMoleculeThreshold:
        :return:
        """

        # Molecule counts and type
        self.baseCounts =  {"DPN": {"A": 0, "C": 0, "G": 0, "T": 0},
                            "DpN": {"A": 0, "C": 0, "G": 0, "T": 0},
                            "DPn": {"A": 0, "C": 0, "G": 0, "T": 0},
                            "Dpn": {"A": 0, "C": 0, "G": 0, "T": 0},
                            "SN" : {"A": 0, "C": 0, "G": 0, "T": 0},
                            "SP" : {"A": 0, "C": 0, "G": 0, "T": 0},
                            "Sp" : {"A": 0, "C": 0, "G": 0, "T": 0},
                            "Sn" : {"A": 0, "C": 0, "G": 0, "T": 0}
                            }
        self.strandCounts = {"A": 0, "C": 0, "G":0, "T": 0}
        self.alleleStrongCounts = {"A": 0, "C": 0, "G":0, "T": 0}
        self.duplexCounts = {"A": 0, "C": 0, "G":0, "T": 0}

        # Parental strand counts
        self.posMolec = {"A": 0, "C": 0, "G": 0, "T": 0}
        self.negMolec = {"A": 0, "C": 0, "G": 0, "T": 0}

        # Map strand counts, for strand bias calculations
        self.pMapStrand = {"A": 0, "C": 0, "G": 0, "T": 0}
        self.nMapStrand = {"A": 0, "C": 0, "G": 0, "T": 0}

        # Handle overlapping read pairs
        readNameIndex = {}
        readNameBaseIndex = {"A": {}, "C": {}, "G": {}, "T": {}}

        # Identify and flag reads in duplex
        # If a read ID occurs once, that read is a singleton
        # If that read ID occurs twice, that read is in duplex
        singletonIndex = {}
        duplexIDs = {}
        depth = 0

        # Since the features of this position are "paired" (i.e. index 0 of each list corresponds to the statistics
        # of the same read), cycle through all stored attributes
        i = -1
        for base, qual, famSize, posParent, map, distToEnd, mismatchNum, mapQual, readID, readName \
                in zip(self.alleles, self.qualities, self.famSizes, self.posParent, self.mapStrand, self.distToEnd,
                       self.mismatchNum, self.mappingQual, self.readNum, self.readName):

            i += 1
            # If the base in "None", then we have started reading into the buffer. There are no more
            # positions to process
            if base is None:
                break

            # What strand do these reads map to? Used to calculate strand bias info later
            if map:  # Negative map strand
                self.nMapStrand[base] += 1
            else:
                self.pMapStrand[base] += 1

            # Have we examined another read with the same name?
            if readName in readNameIndex:
                # If so, this read pair overlaps
                # In this case, only use the base with the highest quality score
                otherIndex = readNameIndex[readName]
                otherQual = self.qualities[otherIndex]  # Obtain the other quality score
                otherBase = self.alleles[otherIndex]  # The other base
                if qual > otherQual:  # This read has a higher quality score
                    # Use this read instead of the other read
                    readNameIndex[readName] = i
                    del readNameBaseIndex[otherBase][readName]
                    readNameBaseIndex[base][readName] = i
                    depth -= 1
                elif qual < otherQual:  # The other read has a higher quality score. Ignore this one
                    continue
                else:
                    # In the case of a tie, use the read with the largest family size
                    otherFSize = self.famSizes[otherIndex]
                    if otherFSize > famSize:
                        continue
                    elif otherFSize < famSize:
                        readNameIndex[readName] = i
                        del readNameBaseIndex[otherBase][readName]
                        readNameBaseIndex[base][readName] = i
                        depth -= 1
                    else:
                        # In the case of a tie, use the reference base to be conservative
                        if otherBase == self.ref:
                            continue
                        elif base == self.ref:
                            readNameIndex[readName] = i
                            del readNameBaseIndex[otherBase][readName]
                            readNameBaseIndex[base][readName] = i
                            depth -= 1
                        else:  # Both bases have the same quality score, and support different alternate alleles
                            # Choose the alternate allele with the most support
                            try:
                                otherWeight = self.altAlleles[otherBase]
                            except KeyError:
                                otherWeight = 0
                            try:
                                weight = self.altAlleles[base]
                            except KeyError:
                                weight = 0
                            if otherWeight > weight:
                                continue
                            elif weight > otherWeight:
                                readNameIndex[readName] = i
                                del readNameBaseIndex[otherBase][readName]
                                readNameBaseIndex[base][readName] = i
                                depth -= 1
                            else:
                                # These bases have the same quality score, and support different alternate alleles which
                                # are equally supported
                                # I give up. Just take the first read encountered
                                continue

            else:
                # This is the first time we have encountered this read name. Store it, in case we encounter the mate
                readNameIndex[readName] = i
                readNameBaseIndex[base][readName] = i

            # Handle duplexes
            # Check if there is another read with this ID
            if readID in singletonIndex:  # This read is in duplex
                # Store this duplex
                duplexIDs[readID] = (readName, self.readName[singletonIndex[readID]])  # Store the readNames
                del singletonIndex[readID]
            else:  # We have not encountered a read with this ID as of yet
                singletonIndex[readID] = i
                depth += 1  # Total molecule count


        # Process duplexes
        for duplexName1, duplexName2 in duplexIDs.values():
            index1 = readNameIndex[duplexName1]
            index2 = readNameIndex[duplexName2]

            base1 = self.alleles[index1]
            base2 = self.alleles[index2]

            # If the bases disagree, then we should discard this duplex
            if base1 != base2:
                del readNameIndex[duplexName1]
                del readNameIndex[duplexName2]
                del readNameBaseIndex[base1][duplexName1]
                del readNameBaseIndex[base2][duplexName2]

                # Update strand and molecule counts
                depth -= 1
                continue

            # What type of duplex is this?
            posStrong = False
            negStrong = False
            if self.posParent[index1]:  # Read 1 of the duplex is from the positive parental strand
                if self.famSizes[index1] > strongMoleculeThreshold:
                    posStrong = True
            else:  # Read 1 of the duplex is from the negative parental strand
                if self.famSizes[index1] > strongMoleculeThreshold:
                    negStrong = True
            # Repeat for the second read
            if self.posParent[index2]:  # Read 2 of the duplex is from the positive parental strand
                if self.famSizes[index2] > strongMoleculeThreshold:
                    posStrong = True
            else:  # Read 2 of the duplex is from the negative parental strand
                if self.famSizes[index2] > strongMoleculeThreshold:
                    negStrong = True

            # Add the duplex to the count
            if posStrong and negStrong:
                self.baseCounts["DPN"][base1] += 1
            elif posStrong:
                self.baseCounts["DPn"][base1] += 1
            elif negStrong:
                self.baseCounts["DpN"][base1] += 1
            else:
                self.baseCounts["Dpn"][base1] += 1

            self.posMolec[base1] += 1
            self.negMolec[base1] += 1

            self.duplexCounts[base1] += 1
            # To prevent this duplex from being double-counted downstream, remove one of the reads from processing
            # and keep the read with the largest family size
            if self.famSizes[index1] > self.famSizes[index2]:
                del readNameIndex[duplexName2]
                del readNameBaseIndex[base2][duplexName2]
            else:
                del readNameIndex[duplexName1]
                del readNameBaseIndex[base1][duplexName1]

        # Process singletons
        for readName, index in singletonIndex.items():

            base = self.alleles[index]
            # What type of molecule is this?
            # What type of molecule is this?
            if self.posParent[index]:
                self.posMolec[base] += 1
                if self.famSizes[index] > strongMoleculeThreshold:
                    self.baseCounts["SP"][base] += 1
                else:
                    self.baseCounts["Sp"][base] += 1
            else:
                self.negMolec[base] += 1
                if self.famSizes[index] > strongMoleculeThreshold:
                    self.baseCounts["SN"][base] += 1
                else:
                    self.baseCounts["Sn"][base] += 1

        # Generate summary stats
        self.alleleMapQual = {}
        self.alleleBaseQual = {}
        self.alleleStrandBias = {}
        self.alleleVafs = {}
        self.molecDepth = sum(len(readNameBaseIndex[x]) for x in readNameBaseIndex)
        self.alleleMismatch = {}
        self.alleleAvFamSize = {}
        self.alleleEndDist = {}

        # Process each variant allele individually
        for base, reads in readNameBaseIndex.items():
            if not reads:  # i.e. the dictionary is empty. There are no reads which support this allele. Use placeholder values
                self.alleleMapQual[base] = (0)
                self.alleleBaseQual[base] = (0)
                self.alleleStrandBias[base] = 1
                self.alleleVafs[base] = 0
                self.alleleMismatch[base] = (0)
                self.alleleAvFamSize[base] = (0)
                self.alleleEndDist[base] = (0)
            else:
                self.alleleMapQual[base] = tuple(self.mappingQual[x] for x in reads.values())
                self.alleleBaseQual[base] = tuple(self.qualities[x] for x in reads.values())
                self.alleleMismatch[base] = tuple(self.mismatchNum[x] for x in reads.values())
                self.alleleAvFamSize[base] = tuple(self.famSizes[x] for x in reads.values())
                self.alleleEndDist[base] = tuple(self.distToEnd[x] for x in reads.values())
                # Calculate strand bias
                self.alleleStrandBias[base] = fisher_exact(self.pMapStrand[base], self.nMapStrand[base],
                                                      sum(self.pMapStrand[x] for x in ("A", "C", "G", "T")),
                                                       sum(self.nMapStrand[x] for x in ("A", "C", "G", "T"))).two_tail

                counts = len(reads.values())
                self.strandCounts[base] = counts
                self.alleleVafs[base] = counts / self.molecDepth


        # If no alternate alleles pass the minimum depth filter, than this position is not a real variant
        failedDepth = []
        for allele in self.altAlleles.keys():
            if self.strandCounts[allele] < minAltDepth:
                failedDepth.append(allele)
        for allele in failedDepth:
            del self.altAlleles[allele]
        if not self.altAlleles:
            return False
        return True

    def leftFlankProp(self, base):
        """
        Calculates the proportion of bases upstream of this variant that which are "base"

        :param base: The nucleotide to use "A,C,G,T"
        :return: A float indication the proportion of flanking bases that are "base"

        Note: For performance reasons, this will not check that base is a nucleotide
        """

        i = 0.0
        for b in self.leftWindow:
            if b == base:
                i += 1

        return i / len(self.leftWindow)

    def rightFlankProp(self, base):
        """
        Calculates the proportion of bases downstream of this variant that which are "base"

        :param base: The nucleotide to use "A,C,G,T"
        :return: A float indication the proportion of right-flanking bases that are "base"

        Note: For performance reasons, this will not check that base is a nucleotide
        """

        i = 0.0
        for b in self.rightWindow:
            if b == base:
                i += 1

        return i / len(self.rightWindow)


class IndelPos(object):
    """
    Similar to a position object, but stores indel events, rather than SNV/SNPs

    Only stores variant positions, not reference positions
    """

    __slots__ = ("alleles", "qualities", "famSizes", "posParent", "mapStrand", "distToEnd", "mismatchNum",
                 "mappingQual", "readNum", "readName", "alt", "altAlleles", "depth", "altDepth",
                 "leftWindow", "ref", "rightWindow", "nearbyVar", "nWindow", "baseCounts",
                 "alleleStrongCounts", "strandCounts", "posMolec", "negMolec", "pMapStrand", "nMapStrand",
                 "alleleMapQual", "alleleBaseQual", "alleleStrandBias", "alleleVafs", "molecDepth", "alleleMismatch",
                 "alleleAvFamSize", "alleleEndDist", "type", "length", "alt", "disagreeingDuplexes", "duplexCounts")

    def __init__(self, type, length, alt=None):

        self.type = type  # Is this an insertion or deletion?
        self.length = length  # How long is this event?

        self.ref = None
        self.alt = alt

        # Statistics coresponding to the reads which support this event
        self.alleles = [None] * 5
        self.qualities = [None] * 5
        self.famSizes = [None] * 5
        self.posParent = [None] * 5
        self.mapStrand = [None] * 5
        self.distToEnd = [None] * 5
        self.mismatchNum = [None] * 5
        self.mappingQual = [None] * 5
        self.readNum = [None] * 5
        self.readName = [None] * 5
        self.depth = 0
        self.altDepth = 0

    def add(self, allele, qual, posParent, famSize, mapStrand, mapQual, dist, readNum, readName, isAlt=False):

        # Add a new base, as well as the corresponding characteristics, to this position
        try:
            self.alleles[self.depth] = allele
        except IndexError:  # We have exceeded the current max buffered positions. Increase the buffer
            self.alleles.extend([None] * 100)
            self.qualities.extend([None] * 100)
            self.famSizes.extend([None] * 100)
            self.posParent.extend([None] * 100)
            self.mapStrand.extend([None] * 100)
            self.distToEnd.extend([None] * 100)
            self.mismatchNum.extend([None] * 100)
            self.mappingQual.extend([None] * 100)
            self.readNum.extend([None] * 100)
            self.readName.extend([None] * 100)
            self.alleles[self.depth] = allele
        self.qualities[self.depth] = qual
        self.famSizes[self.depth] = famSize
        self.posParent[self.depth] = posParent
        self.mapStrand[self.depth] = mapStrand
        self.distToEnd[self.depth] = dist
        self.mappingQual[self.depth] = mapQual
        self.readNum[self.depth] = readNum
        self.readName[self.depth] = readName
        self.depth += 1

        if isAlt:
            self.altDepth += 1

    def addMismatchNum(self, misMatchNum):
        self.mismatchNum[self.depth - 1] = misMatchNum

    def processVariant(self, leftFlank, rightFlank, nearbyVar, noiseWindow):
        """
        Store some basic statistics about this position

        :param leftFlank: A list containing the nucleotides upstream of this position
        :param rightFlank: A list containing the nucleotides downstream of this position
        :param nearbyVar: A list storing all adjacent indels that fall within noiseWindow
        :param noiseWindow: How many bases upstream and downstream were searched for nearby indels
        :return:
        """

        self.leftWindow = leftFlank
        self.rightWindow = rightFlank
        self.nearbyVar = nearbyVar
        self.nWindow = noiseWindow

    def summarizeVariant(self, minAltDepth=4, strongMoleculeThreshold=3, delBaseQual=38):
        """
        Summarize the statistics for this variant
        Generate a by-allele summary of the number of reads/molecules which support each
        variant, base quality, mapping quality, strand bias, and so forth

        This function handles overlapping read pairs. If only one read in the pair supports an indel,
        the other read (i.e. the one that does not support the indel) is used instead.

        :param minAltDepth: The minimum number of molecules which must support an alternate allele before it is filtered out
        :param strongMoleculeThreshold: How many members of a read family before it is considered a strong family?
        :param delBaseQual: What quality score should be assigned to deletions, since there are no bases in deletions?
        """

        # Data structures to hold the various summary statistics
        self.baseCounts =  {"DPN": {"ALT": 0, "REF": 0},
                            "DpN": {"ALT": 0, "REF": 0},
                            "DPn": {"ALT": 0, "REF": 0},
                            "Dpn": {"ALT": 0, "REF": 0},
                            "SN" : {"ALT": 0, "REF": 0},
                            "SP" : {"ALT": 0, "REF": 0},
                            "Sp" : {"ALT": 0, "REF": 0},
                            "Sn" : {"ALT": 0, "REF": 0}
                            }
        self.strandCounts = {"ALT": 0, "REF": 0}
        self.posMolec = {"ALT": 0, "REF": 0}
        self.negMolec = {"ALT": 0, "REF": 0}
        self.disagreeingDuplexes = {"ALT": 0, "REF": 0}
        self.duplexCounts = {"ALT": 0, "REF": 0}

        self.pMapStrand = {"ALT": 0, "REF": 0}
        self.nMapStrand = {"ALT": 0, "REF": 0}

        depth = 0.0

        # If clipOverlap was not run on the BAM file, than we need to account for cases where a single read pair overlaps
        # Store the index of each read name here
        readNameIndex = {}
        readNameBaseIndex = {"ALT": {}, "REF": {}}

        # We need to identify and flag reads that are in duplex
        # Count the number of times each read number "A counter that is unique to each read, generated by collapse" occurs
        # If it occurs once, the read is a singleton. If it occurs twice, it is a duplex
        singletonIndex = {}
        duplexIDs = {}

        # Since the characteristics of each read is in order (i.e. index 0 corresponds to the first read that overlaps
        # this position for every list, cycle through all stored attributes
        i = -1
        for base, qual, fSize, posParent, map, distToEnd, mismatchNum, mapQual, readID, readName \
                in zip(self.alleles, self.qualities, self.famSizes, self.posParent, self.mapStrand, self.distToEnd, self.mismatchNum, self.mappingQual, self.readNum, self.readName):

            # If this base is None, then we have started reading into the buffer, and should stop processing bases
            if base is None:
                break

            i += 1

            # Does this base support the indel?
            # Bases which support indels will have a quality score that is a list, instead of an int
            supportsIndel = isinstance(qual, list)
            if supportsIndel:
                alt = "ALT"
            else:
                alt = "REF"

            # Calculate Strand bias info
            if map:  # i.e. this base originated from a read that mapped to the negative strand
                self.nMapStrand[alt] += 1
            else:
                self.pMapStrand[alt] += 1


            # Check to see if there is already a read with the same name at this position
            # If this is the case, then there is an overlapping read pair.
            # We should take the consensus of this overapping read pair
            if readName in readNameIndex:
                # If one of these reads does not support the indel, discard the read supporting the indel, as
                # it is likely the result of sequencing error
                if not supportsIndel:
                    readNameIndex[readName] = i
                    readNameBaseIndex["REF"][readName] = i
                    try:
                        depth -= 1
                        del readNameBaseIndex["ALT"][readName]
                    except KeyError:
                        pass
                else:
                    continue

            else:
                readNameIndex[readName] = i
                readNameBaseIndex[alt][readName] = i

            # Check to see if there is another read that is the duplex mate of this read
            if readID in singletonIndex:  # There is another read with this ID. It is in duplex
                # Store the read names. We will determine which index these corespond to once all overlaps are resolved
                duplexIDs[readID] = (readName, self.readName[singletonIndex[readID]])
                del singletonIndex[readID]
            else:
                singletonIndex[readID] = i
                # Total molecule count (for VAF)
                depth += 1

        # Now that all duplexes and overlapping read pairs have been identified, sum molecule counts
        for duplexRead1, duplexRead2 in duplexIDs.values():
            index1 = readNameIndex[duplexRead1]
            index2 = readNameIndex[duplexRead2]

            base1 = self.alleles[index1]
            base2 = self.alleles[index2]

            # Does this support an insertion, or deletion?
            if isinstance(self.qualities[index1], int):
                alt1 = "REF"
            else:
                alt1 = "ALT"

            if isinstance(self.qualities[index2], int):
                alt2 = "REF"
            else:
                alt2 = "ALT"

            # If the bases between the duplexes disagree, then discard the duplex, as it is unreliable
            if alt1 != alt2:
                del readNameIndex[duplexRead1]
                del readNameIndex[duplexRead2]
                del readNameBaseIndex[alt1][duplexRead1]
                del readNameBaseIndex[alt2][duplexRead2]
                depth -= 1
                continue
            posStrong = False
            negStrong = False
            # What type of duplex is this?
            if self.posParent[index1]:
                if self.famSizes[index1] > strongMoleculeThreshold:
                    posStrong = True
            else:
                if self.famSizes[index1] > strongMoleculeThreshold:
                    negStrong = True
            if self.posParent[index2]:
                if self.famSizes[index2] > strongMoleculeThreshold:
                    posStrong = True
            else:
                if self.famSizes[index2] > strongMoleculeThreshold:
                    negStrong = True

            base = self.alleles[index1]
            if posStrong and negStrong:
                self.baseCounts["DPN"][alt1] += 1
            elif posStrong:
                self.baseCounts["DPn"][alt1] += 1
            elif negStrong:
                self.baseCounts["DpN"][alt1] += 1
            else:
                self.baseCounts["Dpn"][alt1] += 1

            self.posMolec[alt1] += 1
            self.negMolec[alt1] += 1


            # Count this duplex
            self.duplexCounts[alt1] += 1

            # To avoid double-counting this molecule during downstream processing, remove one of the read names from
            # the read ID index
            # Use the family with the largest size
            if self.famSizes[index1] > self.famSizes[index2]:
                del readNameIndex[duplexRead2]
                del readNameBaseIndex[alt2][duplexRead2]

            else:
                del readNameIndex[duplexRead1]
                del readNameBaseIndex[alt1][duplexRead1]

        # Add singleton counts
        for readName, index in singletonIndex.items():
            alt = "REF" if isinstance(self.qualities[index], int) else "ALT"

            # What type of molecule is this?
            if self.posParent[index]:
                self.posMolec[alt] += 1
                if self.famSizes[index] > strongMoleculeThreshold:
                    self.baseCounts["SP"][alt] += 1
                else:
                    self.baseCounts["Sp"][alt] += 1
            else:
                self.negMolec[alt] += 1
                if self.famSizes[index] > strongMoleculeThreshold:
                    self.baseCounts["SN"][alt] += 1
                else:
                    self.baseCounts["Sn"][alt] += 1

        # Generate summary stats
        self.alleleMapQual = {}
        self.alleleBaseQual = {}
        self.alleleStrandBias = {}
        self.alleleVafs = {}
        self.alleleMismatch = {}
        self.alleleAvFamSize = {}
        self.alleleEndDist = {}
        self.molecDepth = sum(len(readNameBaseIndex[x]) for x in readNameBaseIndex)
        # Since this is an indel, we need to generate a consensus for the indel sequence
        # The consensus will be the most frequent size, and the consensus sequence (insertions only)
        refSeq = ""
        altSeq = ""
        altQual = []

        for base, reads in readNameBaseIndex.items():

            if not reads:  # i.e. the dictionary is empty. There are no values to process
                self.alleleMapQual[base] = [0]
                self.alleleBaseQual[base] = [0]
                self.alleleStrandBias[base] = 1
                self.alleleVafs[base] = 0
                self.alleleMismatch[base] = [0]
                self.alleleAvFamSize[base] = [0]
                self.alleleEndDist[base] = [0]
            else:
                if base == "ALT":
                    # Is this an insertion? Or deletion
                    if self.type == "D":  # Deletion
                        # The reference sequence is stored in "alleles"
                        # Find an instance where the reference seq is long enough
                        for index in reads.values():
                            try:
                                refSeq = self.alleles[index][0:self.length]
                                break
                            except IndexError:
                                continue
                        self.alleleBaseQual[base] = delBaseQual
                    else:  # Insertion
                        # The inserted seq is stored in "alleles"
                        # Generate a consensus
                        consensus = []
                        consensusQual = []
                        for index in reads.values():  # Loop through each read and summarize all the stats
                            i = 0
                            for b, q in zip(self.alleles[index], self.qualities[index]):
                                if i == len(consensus):
                                    consensus.append(sortedcontainers.SortedDict())
                                    consensusQual.append(sortedcontainers.SortedDict())
                                if b not in consensus[i]:
                                    consensus[i][b] = 1
                                    consensusQual[i][b] = q
                                else:
                                    consensus[i][b] += 1
                                    if q > consensusQual[i][b]:
                                        consensusQual[i][b] = q
                                i += 1
                        # Now use the most frequent base at each index as the consensus for the insertion
                        try:
                            for i in range(0, self.length):
                                maxBase = ""
                                maxCount = 0
                                for iBase, count in consensus[i].items():
                                    if count > maxCount:
                                        maxBase = iBase
                                        maxCount = count
                                altSeq += maxBase
                                altQual.append(consensusQual[i][maxBase])
                        except IndexError:  # i.e. the reads which supported the indel were not counted, likely due to duplex handling
                            return False
                        # Take an average of the quality scores
                        self.alleleBaseQual[base] = int(mean(altQual))
                else:
                    self.alleleBaseQual[base] = tuple(self.qualities[x] for x in reads.values())
                self.alleleMapQual[base] = tuple(self.mappingQual[x] for x in reads.values())

                # Calculate strand bias
                self.alleleStrandBias[base] = fisher_exact(self.pMapStrand[base], self.nMapStrand[base], sum(self.pMapStrand[x] for x in ("ALT", "REF")), sum(self.nMapStrand[x] for x in ("ALT", "REF"))).two_tail

                counts = len(reads.values())
                self.strandCounts[base] = counts
                self.alleleVafs[base] = counts / self.molecDepth

                # Aggregate remaining stats
                self.alleleMismatch[base] = tuple(self.mismatchNum[x] for x in reads.values())
                self.alleleAvFamSize[base] = tuple(self.famSizes[x] for x in reads.values())
                self.alleleEndDist[base] = tuple(self.distToEnd[x] for x in reads.values())

        # Finally, finalize the reference and alternate alleles for this event
        if self.type == "D":
            # A deletion
            self.alt = self.ref
            self.ref = self.ref + refSeq
        else:
            # An insertion
            self.alt = self.ref + altSeq

        # Basic depth filtering
        if self.strandCounts["ALT"] < minAltDepth:
            return False
        return True

    def leftFlankProp(self, base):
        """
        Calculates the proportion of bases upstream of this variant that which are "base"

        :param base: The nucleotide to use "A,C,G,T"
        :return: A float indication the proportion of flanking bases that are "base"

        Note: For performance reasons, this will not check that base is a nucleotide
        """

        i = 0.0
        for b in self.leftWindow:
            if b == base:
                i += 1

        return i / len(self.leftWindow)

    def rightFlankProp(self, base):
        """
        Calculates the proportion of bases downstream of this variant that which are "base"

        :param base: The nucleotide to use "A,C,G,T"
        :return: A float indication the proportion of right-flanking bases that are "base"

        Note: For performance reasons, this will not check that base is a nucleotide
        """

        i = 0.0
        for b in self.rightWindow:
            if b == base:
                i += 1

        return i / len(self.rightWindow)


class PileupEngine(object):
    """
    Generates a custom pileup using family characteristics
    """

    def __init__(self, inBAM, refGenome, targetRegions, minAltDepth=1, homopolymerWindow=7, noiseWindow=150,
                 pileupWindow=1000, oBAM=None, normalBAM=None, printPrefix="PRODUSE-CALL\t", softClipUntilIndel=25):
        try:
            self._inFile = pysam.AlignmentFile(inBAM, require_index=True)
            if normalBAM:
                self._normalFile = pysam.AlignmentFile(normalBAM, require_index=True)
            else:
                self._normalFile = None
        except FileNotFoundError as e:
            raise pe.IndexNotFound("Unable to locate BAM index \'%s.bai\': No such file or directory" % inBAM) from e
        self._refGenome = Fasta(refGenome, read_ahead=20000)
        self._captureSpace = self._loadCaptureSpace(targetRegions)
        self.candidateIndels = {}
        self.rawIndels = {}
        self.pileup = {}
        self.candidateVar = {}
        self.filteredVar = {}
        self._minAltDepth = minAltDepth
        self._homopolymerWindow = homopolymerWindow
        self._noiseWindow = noiseWindow
        self._softClipUntilIndel = softClipUntilIndel  # How many bases must be consecutively soft clipped before realignment is performed?
        # The number of positions to store before collapsing previous positions
        # Reduce this number to lower memory footprint, but be warned that, if this falls below the length of a read,
        # Some positions may not be tallied correctly
        self._pileupWindow = pileupWindow

        self._refStart = 0
        self._chrom = None
        self._refWindow = ""

        self._bufferPos = 0
        self._realignBuffer = 600

        # For local realignment
        self._indelReads = []
        self._bufferedReads = []

        # Debugging
        if oBAM is not None:
            self._oBAM = pysam.AlignmentFile(oBAM, mode="wb", template=self._inFile)
        else:
            self._oBAM = None

        # Status updates
        self._printPrefix = printPrefix
        self.readCount = 0
        self.varCount = 0

        # Output vcf genotype fields
        self._genotypeFormat = "DP:AD:ADF:ADR"

    def reset(self):
        """
        Resets all positions and buffers back to 0
        :return:
        """
        self._refStart = 0
        self._chrom = None
        self._bufferPos = 0
        self._refWindow = ""

    def _loadCaptureSpace(self, file):
        """
        Parse genomic regions from a specified BED file, and load them into a dictionary

        :param file: A string containing a filepath to the BED file which specifies the capture space
        :return: A dictionary listing "chrom" : (start, end, start, end, start, end)
        """

        if file is None:
            return None

        # Does the BAM file use "chr"-prefixed contig names? If so, we should make sure the names of these regions are consistent
        isChrPrefixed = None
        try:
            # Check the first reference sequence listed in the BAM header. This isn't perfect, but it should handle most cases
            isChrPrefixed = self._inFile.header["SQ"][0]["SN"].startswith("chr")
        except AssertionError:
            pass

        try:
            targets = {}
            with open(file) as f:
                for line in f:
                    elements = line.split("\t")
                    chrom = elements[0]
                    # Make sure the "chr" prefix used here is consistent with the BAM file
                    if isChrPrefixed is not None:
                        if isChrPrefixed:
                            if not chrom.startswith("chr"):
                                chrom = "chr" + chrom
                        else:
                            if chrom.startswith("chr"):
                                chrom = chrom[3:]

                    # Since BED files use 1-based numbering, while pysam uses 0, include an offset
                    start = int(elements[1]) - 1
                    end = int(elements[2]) - 1
                    if chrom not in targets:
                        targets[chrom] = []
                    targets[chrom].extend([start, end])
                    # Ensure that the BED interval is actually valid
                    if start > end:
                        sys.stderr.write(
                            "ERROR: The start position \'%s\' is greater than the end position \'%s\'\n" % (
                                start, end))
                        exit(1)
            # Finally, to speed up access, convert each list into a tuple
            for chrom, locations in targets.items():
                targets[chrom] = tuple(locations)
            return targets
        except (IndexError, ValueError):
            sys.stderr.write("ERROR: Unable to parse BED interval \'%s\'\n" % (line.strip("\n").strip("\r")))
            sys.stderr.write("Ensure the file is a valid BED file, and try again\n")
            exit(1)

    def _insideCaptureSpace(self, read):
        """
        Does this read fall within our capture space? Note that this excludes soft-clipping
        :param read: A pysam.AlignedSegment() object, coresponding to the read of interest
        :return: A boolean indicating if this read overlaps the capture space
        """
        # If no capture space was specified, then considering this read "within" the capture space
        if self._captureSpace is None:
            return True
        # Check Chromosome
        if read.reference_name not in self._captureSpace:
            return False

        # Check read coordinates
        bisectPoint1 = bisect.bisect_left(self._captureSpace[read.reference_name], read.reference_start)
        bisectPoint2 = bisect.bisect_left(self._captureSpace[read.reference_name], read.reference_end)
        if bisectPoint1 == bisectPoint2 and bisectPoint1 % 2 != 1:
            return False
        return True

    def _writeVcfHeader(self, file, filtThreshold):
        """
        Prints a VCF header to the specified file

        :param file: A file object, open to writing
        :param filtThreshold: A float indicating the variant call confidence threshold for the random forest filter
        """

        header = ["##fileformat=VCFv4.3",  # Mandatory
                  "##source=ProDuSe",
                  "##source_version=" + pVer.__version__,
                  "##analysis_date=" + time.strftime('%Y%m%d'),
                  "##reference=" + self._refGenome.filename]

        # Add the contig information
        for contig, features in self._refGenome.faidx.index.items():
            line = "##contig=<ID=" + contig + ",length=" + str(len(features)) + ">"
            header.append(line)

        # Add the info field info
        header.extend([
            '##INFO=<ID=DPN,Number=R,Type=Integer,Description="Duplex Support with Strong Positive and Strong Negative Consensus">',
            '##INFO=<ID=DPn,Number=R,Type=Integer,Description="Duplex Support with Strong Positive and Weak Negative Consensus">',
            '##INFO=<ID=DpN,Number=R,Type=Integer,Description="Duplex Support with Weak Positive and Strong Negative Consensus">',
            '##INFO=<ID=Dpn,Number=R,Type=Integer,Description="Duplex Support with Weak Positive and Weak Negative Consensus">',
            '##INFO=<ID=SP,Number=R,Type=Integer,Description="Singleton Support with Strong Positive Consensus">',
            '##INFO=<ID=Sp,Number=R,Type=Integer,Description="Singleton Support with Weak Positive Consensus">',
            '##INFO=<ID=SN,Number=R,Type=Integer,Description="Singleton Support with Strong Negative Consensus">',
            '##INFO=<ID=Sn,Number=R,Type=Integer,Description="Singleton Support with Weak Negative Consensus">',
            '##INFO=<ID=MC,Number=R,Type=Integer,Description="Total Molecule counts for Each Allele">',
            '##INFO=<ID=STP,Number=R,Type=Float,Description="Probability of Strand Bias during sequencing for each Allele">',
            '##INFO=<ID=PC,Number=R,Type=Integer,Description="Positive Strand Molecule Counts">',
            '##INFO=<ID=NC,Number=R,Type=Integer,Description="Negative Strand Molecule Counts">',
            '##INFO=<ID=VAF,Number=R,Type=Float,Description="Variant allele fraction of alternate allele(s) at this locus">',
            '##INFO=<ID=FILT,Number=R,Type=Float,Description="Variant call confidence for alternate allele(s)">',
            '##FILTER=<ID=PASS,Description="Variant call confidence is above random forest filter confidence threshold %s">' % (filtThreshold),
            '##FILTER=<ID=LOW_CONF,Description="Variant call confidence is below the random forest filter confidence threshold %s">' % (filtThreshold),
            '##FILTER=<ID=NO_DUPLEX,Description="Variant does not have duplex support. Only used when the duplex filter is enabled">',
            '##FILTER=<ID=REPEAT,Description="Variant indel is an expansion or contraction of an adjacent repeat">',
            '##FILTER=<ID=GERMLINE,Description="Variant allele frequency is not significantly higher in the tumour relative to the normal sample">',
            '##FILTER=<ID=NO_DEPTH_NORMAL,Description="Insufficient depth in the normal to assign a variant as germline or somatic">',
            ])

        # Genotype fields
        header.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total MOLECULE count">')
        header.append('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Total number of MOLECULES supporting each allele">')
        header.append('##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Number of READS mapped to the forward strand for each allele">')
        header.append('##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Number of READS mapped to the reverse strand for each allele">')

        header.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "TUMOR" + os.linesep]))
        file.write(os.linesep.join(header))


    def varToFilteringStats(self, pos, allele):
        """
        Summarizes various statistics that a variant is to be filtered based upon

        :param pos: A Position or IndelPos object to be filtered
        :param allele: A string containing a nucleotide, or (for IndelPos) "REF" or "ALT", which is the allele to be filtered
        :return: A tuple containing the various variant stats
        """

        # If this is an indel, we need to set the reference base as "REF"
        if isinstance(pos, IndelPos):
            ref = "REF"
            baseQualBias = 1  # Isn't really applicable for indels
            meanBaseQual = 38  # Since a deletion does not have a base quality, and base quality doesn't apply for indels either
        else:
            ref = pos.ref

            # Mean base quality
            meanBaseQual = mean(pos.alleleBaseQual[allele])
            # Calculate differences in base quality btwn reference and alternate alleles
            try:
                baseQualBias = ttest_ind(pos.alleleBaseQual[ref], pos.alleleBaseQual[allele], equal_var=False).pvalue
                if baseQualBias != baseQualBias:  #i.e the base quality bias is NaN. Set it to 1
                    baseQualBias = 1
            except FloatingPointError:
                baseQualBias = 1

        # Calculate differences in read mapping quality between reference and alternate alleles
        try:
            mapQualBias = ttest_ind(pos.alleleMapQual[ref], pos.alleleMapQual[allele], equal_var=False).pvalue
            if mapQualBias != mapQualBias:
                mapQualBias = 1
        except FloatingPointError:
            mapQualBias = 1

        # Quantify differences in the number of mismatches supported by each read
        try:
            mismatchBias = ttest_ind(pos.alleleMismatch[ref], pos.alleleMismatch[allele], equal_var=False).pvalue
            if mismatchBias != mismatchBias:
                mismatchBias = 1
        except FloatingPointError:
            mismatchBias = 1

        # Quantify differences between the family size supporting the reference and alternate alleles
        try:
            fSizeBias = ttest_ind(pos.alleleAvFamSize[ref], pos.alleleAvFamSize[allele], equal_var=False).pvalue
            if fSizeBias != fSizeBias:
                fSizeBias = 1
        except FloatingPointError:
            fSizeBias = 1

        # Does this mutation indicate DNA damage?
        # C -> A and G -> T mutations are commonly affiliated with DNA damage
        try:
            if (pos.ref == "C" and "A" in pos.altAlleles) or (pos.ref == "G" and "T" in pos.altAlleles):
                dnaDamageMut = True
            else:
                dnaDamageMut = False
        except AttributeError:  # i.e. this is an indel. This filter does not apply
            dnaDamageMut = False


        posStats = (
            pos.strandCounts[allele],
            pos.alleleStrandBias[allele],  # Strand bias pVal
            meanBaseQual,  # Average base qual
            baseQualBias,
            mean(pos.alleleMapQual[allele]),  # Average read mapping quality
            mapQualBias,
            pos.molecDepth,
            max(pos.leftFlankProp(allele), pos.rightFlankProp(allele)),
            len(pos.nearbyVar) / pos.nWindow,
            mean(pos.alleleMismatch[allele]),
            mismatchBias,
            mean(pos.alleleAvFamSize[allele]),
            fSizeBias,
            mean(pos.alleleEndDist[allele]),
            pos.duplexCounts[allele],
            dnaDamageMut
        )

        return posStats

    def filterAndWriteVariants(self, outFile, filter, unfilteredOut=None, filtThreshold=0.6, onlyDuplex=False, indelRepeatThresh=4, writeHeader=False):
        """

        Filters candidate variants based upon specified characteristics, and writes the output to the output file
        This is a placeholder for the real filter, which will be developed at a later time

        :param outFile: A string containing a filepath to the output VCF file, to which variants will be written
        :param filter: A sklearn.ensemble.RandomForestClassifier() which is to be used to filter variants
        :param unfilteredOut: A sting containing a filepath to an output VCF file which will contain ALL variants, even those that do not pass filters
        :param filtThreshold: A float representing the classification threshold in which to filter variants
        :param onlyDuplex: A boolean indicating if only variants with duplex support should pass filters. Not recommended
        :param indelRepeatThresh: An int indicating the number of times a repeat has to occur before an indel which is an expansion/contraction of this repeat is filtered out
        :param writeHeader: A boolean. If true, any existing file outFile will be overwritten, and a VCF header will be added to the file
        :return:
        """

        def _generateVCFEntry(pos, chrom, start, filterCon, filterField):
            """
            Generates a VCF entry of the specified variant

            See https://samtools.github.io/hts-specs/VCFv4.3.pdf for more details

            :param pos: A Position or IndelPos object
            :param chrom: A string containing the reference name (i.e. chromosome) of this variant
            :param start: An int corresponding to the start position of this variant
            :param filterCon: A string (or iterable of strings) indicating how confident the filter is that this variant is real (btwn 0 and 1)
            :param filterField: A string to be used for the info field
            :return: A string coresponding to the VCF entry of the provided position
            """

            # Is this an indel?
            if isinstance(pos, IndelPos):
                # If so, we only need to examine one variant
                qual = pos.alleleBaseQual["ALT"]
                alleles = ("REF", "ALT")
                alt = pos.alt
            else:
                # Otherwise, we need to arrange the alternate alleles by their frequency
                altAlleles = []
                altWeight = []
                for allele in pos.altAlleles.keys():
                    weight = pos.strandCounts[allele]
                    if not altWeight:  # This is the first allele we are examining
                        altAlleles.append(allele)
                        altWeight.append(weight)
                    else:
                        # Otherwise, sort the alternate alleles
                        i = 0
                        for oWeight in altWeight:
                            if weight > oWeight:
                                break
                            i += 1
                        altWeight.insert(i, weight)
                        altAlleles.insert(i, allele)
                alt = ",".join(altAlleles)
                alleles = [pos.ref]
                alleles.extend(altAlleles)
                qual = int(mean(pos.alleleBaseQual[altAlleles[0]]))

            # Generate an info column for this variant
            info = []
            for molecule, counts in pos.baseCounts.items():
                molecFields = molecule + "=" + ",".join(str(counts[x]) for x in alleles)
                info.append(molecFields)

            info.append("STP=" + ",".join(str(pos.alleleStrandBias[x]) for x in alleles))  # Strand bias
            info.append("VAF=" + ",".join(str(pos.alleleVafs[x])[:6] for x in alleles))  # VAF, truncate to 4 decimal places
            info.append("FILT=" + ",".join(str(x) for x in filterCon))

            # Generate genotype fields
            genotype = [
                str(int(pos.molecDepth)),
                ",".join(str(pos.strandCounts[x]) for x in alleles),
                ",".join(str(pos.pMapStrand[x]) for x in alleles),
                ",".join(str(pos.nMapStrand[x]) for x in alleles)
            ]
            variant =   [chrom,  # CHROM
                        str(start + 1),   # POS, add 1 to account for differences in indexing btwn pysam and vcf spec
                        ".",     # ID  (No ID)
                        pos.ref, # REF
                        alt,  # ALT
                        str(qual),     # QUAL
                        filterField,   # FILTER
                        ";".join(info)  ,   # INFO
                         self._genotypeFormat,  # Genotype format field
                         ":".join(genotype) # GENOTYPE
                        ]

            return "\t".join(variant) + os.linesep

        if writeHeader:
            mode = "w"  # i.e. overwrite existing files
        else:
            mode = "a"

        with open(outFile, mode) as o, open(unfilteredOut if unfilteredOut is not None else os.devnull, mode) as u:

            if writeHeader:
                # Write the standard VCF file header to the output file(s)
                self._writeVcfHeader(o, filtThreshold)
                self._writeVcfHeader(u, filtThreshold)

            # Start processing variants
            for chrom, positions in self.candidateVar.items():
                loci = tuple(positions.keys())
                for position in loci:
                    candidateSNV = positions[position]
                    # Is there an indel that occurs before this variant? If so, we should process that
                    # to maintain sorted order of the output VCF files
                    try:
                        posToDelete = []
                        for iPos in self.candidateIndels[chrom].irange(maximum = position, inclusive=(True, False)):
                            indel = self.candidateIndels[chrom][iPos]
                            altAllele = indel.summarizeVariant()

                            if altAllele:  # Passed basic alt depth filters

                                # Filter indel
                                stats = [self.varToFilteringStats(indel, "ALT")]
                                filterResults = filter.predict_proba(stats)[0][0]

                                passesConfFilter = True
                                duplexSupportFilt = True
                                # Set the appropriate FILTER attribute for this variant
                                # PASS = Passes both filters
                                # LOW_CONF = Low random forest confidence
                                # NO_DUPLEX = No duplex support
                                filterField = "PASS"
                                if filterResults < filtThreshold:
                                    passesConfFilter = False
                                    filterField = "LOW_CONF"

                                if (onlyDuplex and indel.duplexCounts["ALT"] == 0):  # Are we filtering based on duplex support?
                                    duplexSupportFilt = False
                                    if filterField == "LOW_CONF":
                                        filterField = filterField + ";NO_DUPLEX"
                                    else:
                                        filterField = "NO_DUPLEX"

                                # Check to see if this indel is an expansion or contraction of a repeat
                                if len(indel.ref) > 1:  # i.e. a deletion
                                    isRepeat = self.checkForRepeat(chrom, iPos, indel.ref[1:], indelRepeatThresh, True)
                                else:
                                    isRepeat = self.checkForRepeat(chrom, iPos, indel.alt[1:], indelRepeatThresh, False)
                                if isRepeat:  # This variant is likely a sequencing artifact
                                    passesConfFilter = False
                                    if filterField == "PASS":
                                        filterField = "REPEAT"
                                    else:
                                        filterField = filterField + ";REPEAT"

                                # Is this event germline?
                                isGermline = False
                                if self._normalFile is not None:
                                    isGermline = self._checkIfGermlineIndel(chrom, iPos, indel)
                                    if isGermline is True:
                                        if filterField == "PASS":
                                            filterField = "GERMLINE"
                                        else:
                                            filterField = filterField + ";GERMLINE"
                                    elif isGermline is None:  # Insufficient depth in normal.
                                        if filterField == "PASS":
                                            filterField = "NO_DEPTH_NORMAL"
                                        else:
                                            filterField = filterField + ";NO_DEPTH_NORMAL"
                                # Write out the unfiltered variant
                                vcfEntry = _generateVCFEntry(indel, chrom, iPos, [str(filterResults)], filterField)
                                u.write(vcfEntry)

                                if passesConfFilter and duplexSupportFilt and isGermline is False:  # This variant passes filters
                                    o.write(vcfEntry)

                            posToDelete.append(iPos)
                            self.varCount += 1
                            if self.varCount % 10000 == 0:
                                sys.stderr.write(
                                    "\t".join([self._printPrefix, time.strftime('%X'),
                                               chrom + ":Positions Filtered:" + str(self.varCount) + "\n"]))

                        # Delete to reduce memory usage
                        for iPos in posToDelete:
                            del self.candidateIndels[chrom][iPos]
                    except KeyError:  # i.e. There are no indels on this chromosome
                        pass

                    altAllele = candidateSNV.summarizeVariant()
                    self.varCount += 1
                    if altAllele:

                        # Summarize the filter stats of each allele
                        alleleStats = tuple(self.varToFilteringStats(candidateSNV, allele) for allele in candidateSNV.altAlleles.keys())

                        filterResults = filter.predict_proba(alleleStats)

                        # Filter alleles
                        failedAlleles = []
                        allFiltResults = []
                        allFiltConf = []
                        passFiltConf = []
                        passFiltResults = []
                        for result, allele in zip(filterResults, candidateSNV.altAlleles.keys()):
                            result = result[0]
                            passesConfFilter = True
                            duplexSupportFilt = True
                            # PASS = Passes both filters
                            # LOW_CONF = Low random forest confidence
                            # NO_DUPLEX = No duplex support
                            filterField = "PASS"
                            if result < filtThreshold:  # Fails filter confidence threshold
                                passesConfFilter = False
                                filterField = "LOW_CONF"
                            if onlyDuplex and candidateSNV.duplexCounts[allele] == 0:
                                duplexSupportFilt = False
                                if filterField == "LOW_CONF":
                                    filterField = filterField + ";NO_DUPLEX"
                                else:
                                    filterField = "NO_DUPLEX"
                            # Is this variant germline? Check if it's in the normal sample (if provided)
                            isGermline = False
                            if self._normalFile is not None:
                                isGermline = self._checkIfGermlineSNV(chrom, position, allele, candidateSNV)
                                if isGermline is True:
                                    if filterField == "PASS":
                                      filterField = "GERMLINE"
                                    else:
                                        filterField = filterField + ";GERMLINE"
                                elif isGermline is None:  # Insufficient depth in normal.
                                    if filterField == "PASS":
                                        filterField = "NO_DEPTH_NORMAL"
                                    else:
                                        filterField = filterField + ";NO_DEPTH_NORMAL"

                            if not passesConfFilter or not duplexSupportFilt or isGermline is True or isGermline is None:  # Fails one or more filters
                                failedAlleles.append(allele)
                            else:
                                passFiltResults.append(filterField)
                                passFiltConf.append(result)
                            allFiltResults.append(filterField)
                            allFiltConf.append(result)

                        vcfEntry = _generateVCFEntry(candidateSNV, chrom, position, allFiltConf, ";".join(allFiltResults))
                        u.write(vcfEntry)

                        # Remove failed alleles
                        for fAllele in failedAlleles:
                            del candidateSNV.altAlleles[fAllele]

                        # If any alt alleles passed filters, print them out
                        if candidateSNV.altAlleles:
                            vcfEntry = _generateVCFEntry(candidateSNV, chrom, position, passFiltConf, ";".join(passFiltResults))
                            o.write(vcfEntry)

                    del positions[position]  # Delete variants after processing them to reduce memory consumption

                    # Status messages
                    if self.varCount % 10000 == 0:
                        sys.stderr.write(
                            "\t".join([self._printPrefix, time.strftime('%X'),
                                       chrom + ":Positions Filtered:" + str(self.varCount) + "\n"]))
                # Process any remaining indels on this chromosome
                try:
                    for iPos, indel in self.candidateIndels[chrom].items():
                        altAllele = indel.summarizeVariant()

                        if altAllele:
                            # Filter indel
                            stats = [self.varToFilteringStats(indel, "ALT")]
                            filterResults = filter.predict_proba(stats)[0][0]

                            passesConfFilter = True
                            duplexSupportFilt = True
                            # Set the appropriate FILTER attribute for this variant
                            # PASS = Passes both filters
                            # LOW_CONF = Low random forest confidence
                            # NO_DUPLEX = No duplex support
                            filterField = "PASS"
                            if filterResults < filtThreshold:
                                passesConfFilter = False
                                filterField = "LOW_CONF"

                            if (onlyDuplex and indel.duplexCounts[
                                "ALT"] == 0):  # Are we filtering based on duplex support?
                                duplexSupportFilt = False
                                if filterField == "LOW_CONF":
                                    filterField = filterField + ";NO_DUPLEX"
                                else:
                                    filterField = "NO_DUPLEX"

                            # Check to see if this indel is an expansion or contraction of a repeat
                            if len(indel.ref) > 1:  # i.e. a deletion
                                isRepeat = self.checkForRepeat(chrom, iPos, indel.ref[1:], indelRepeatThresh, True)
                            else:
                                isRepeat = self.checkForRepeat(chrom, iPos, indel.alt[1:], indelRepeatThresh, False)
                            if isRepeat:  # This variant is likely a sequencing artifact
                                passesConfFilter = False
                                if filterField == "PASS":
                                    filterField = "REPEAT"
                                else:
                                    filterField = filterField + ";REPEAT"

                            # Is this event germline?
                            isGermline = False
                            if self._normalFile is not None:
                                isGermline = self._checkIfGermlineIndel(chrom, iPos, indel)
                                if isGermline is True:
                                    if filterField == "PASS":
                                        filterField = "GERMLINE"
                                    else:
                                        filterField = filterField + ";GERMLINE"
                                elif isGermline is None:  # Insufficient depth in normal.
                                    if filterField == "PASS":
                                        filterField = "NO_DEPTH_NORMAL"
                                    else:
                                        filterField = filterField + ";NO_DEPTH_NORMAL"

                            # Write out the unfiltered variant
                            vcfEntry = _generateVCFEntry(indel, chrom, iPos, [str(filterResults)], filterField)
                            u.write(vcfEntry)

                            if passesConfFilter and duplexSupportFilt and isGermline is False:  # This variant passes filters
                                o.write(vcfEntry)

                    del self.candidateIndels[chrom]
                except KeyError:
                    pass
            # Finally, process any indels on remaining contigs
            for chrom in self.candidateIndels:
                for iPos, indel in self.candidateIndels[chrom].items():
                    hasAlt = indel.summarizeVariant()
                    if hasAlt:
                        # Filter indel
                        stats = [self.varToFilteringStats(indel, "ALT")]
                        filterResults = filter.predict_proba(stats)[0][0]

                        passesConfFilter = True
                        duplexSupportFilt = True
                        # Set the appropriate FILTER attribute for this variant
                        # PASS = Passes both filters
                        # LOW_CONF = Low random forest confidence
                        # NO_DUPLEX = No duplex support
                        filterField = "PASS"
                        if filterResults < filtThreshold:
                            passesConfFilter = False
                            filterField = "LOW_CONF"

                        if (onlyDuplex and indel.duplexCounts[
                            "ALT"] == 0):  # Are we filtering based on duplex support?
                            duplexSupportFilt = False
                            if filterField == "LOW_CONF":
                                filterField = filterField + ",NO_DUPLEX"
                            else:
                                filterField = "NO_DUPLEX"

                        # Check to see if this indel is an expansion or contraction of a repeat
                        if len(indel.ref) > 1:  # i.e. a deletion
                            isRepeat = self.checkForRepeat(chrom, iPos, indel.ref[1:], indelRepeatThresh, True)
                        else:
                            isRepeat = self.checkForRepeat(chrom, iPos, indel.alt[1:], indelRepeatThresh, False)
                        if isRepeat:  # This variant is likely a sequencing artifact
                            passesConfFilter = False
                            if filterField == "PASS":
                                filterField = "REPEAT"
                            else:
                                filterField = filterField + ",REPEAT"
                        # Write out the unfiltered variant
                        vcfEntry = _generateVCFEntry(indel, chrom, iPos, [str(filterResults)], filterField)
                        u.write(vcfEntry)

                        if passesConfFilter and duplexSupportFilt:  # This variant passes filters
                            o.write(vcfEntry)

            self.candidateVar = {}
            self.candidateIndels = {}

    def checkForRepeat(self, chrom, position, seq, threshold, offsetSearch):
        """
        Check to see if a given sequence occurs at a given position in the reference more than x number of times
        :param chrom: A string listing the name of the contig to be examined
        :param position: An int listing the position in the contig to be examined
        :param seq: The sequence of the variant to be checked
        :param threshold: An int listing the number of times the variant should occur
        :param offsetSearch: A boolean. When searching, offset the starting position by the length of seq
        :return: A boolean indicating if a repeat occurs
        """

        seqLength = len(seq)
        if offsetSearch:
            position = position + seqLength

        endPos = position + seqLength * threshold
        # Obtain the reference sequence to be examined for repeats
        try:
            refSeq = self._refGenome[chrom][position:endPos].seq
        except IndexError:  # i.e. this indel occurs near the end of the genome. This filter does not apply
            return False

        # Now, check each slice for the repeats
        for j in range(threshold):
            cPos = j * seqLength
            cEnd = cPos + seqLength
            if refSeq[cPos:cEnd] != seq: # This sequence does not match. This sequence does not occur the specified number of times
                return False

        # The specified seq occurs more than threshold number of times
        return True

    def _checkIfGermlineSNV(self, chrom, position, altAllele, altPosition, pValThresh=0.05, minReadThreshold=6, minVafThreshold=0.33, maxVafThreshold=0.08):
        """
        Determines if a given snv is a germline mutation

        Compares the number of reads which support the alternate allele between the tumour sample and the matched normal
        sample using a fisher's exact test. If the number of reads supporting the alternate allele in the tumour is not
        significantly higher than in the normal, this variant is flagged as germline

        WARNING: This implementation is VERY basic for performance and code maintinance readsons. It will not take into
        account family sizes or any duplexes. That said, if the normal BAM was collapsed, each family will only be
        counted once

        :param chrom: A string indicating the name of the contig the read is mapped against
        :param position: An integer indicating the position of the variant in the specified contig
        :param altAllele: A string containing the alternate allele
        :param altPosition: A position() object containing the variant position and aggregated stats
        :param pValThresh: A float indicating the p-value threshold. If the fisher's pvalue is above this threshold, the variant is flagged as germline
        :param minReadThreshold: An int indicating the minimum normal depth required to confidently assign a variant as germline or not
        :param minVafThreshold: A float. If the VAF in the normal sample is higher than this threshold, the variant is flagged as germline
        :param maxVafThreshold: A flat. If the normal sample VAF is lower than this threshold, the variant is flagged as somatic
        :return: A boolean indicating if the variant is germline, or None, if there is insufficient depth
        """

        # Count the number of reads which support the reference and alternate alleles in the normal
        refAllele = altPosition.ref
        refCount = 0
        altCount = 0
        pileup = self._normalFile.pileup(contig=chrom, start=position, stop=position+1, maxDepth=20000)
        for pileupPos in pileup:
            # Since pysam's pileup generates an iterable of all positions covered by reads which overlap the specified
            # position, we need to go to the position we care about
            if pileupPos.pos == position:
                # Obtain the base by all reads which overlap this position
                for read in pileupPos.pileups:
                    # If this read has a deletion at this position, don't count it
                    if read.query_position is None:
                        continue

                    base = read.alignment.query_sequence[read.query_position]
                    if base == refAllele:
                        refCount += 1
                    elif base == altAllele:
                        altCount += 1

        # Is there sufficient normal depth?
        if refCount + altCount < minReadThreshold:
            return None

        # Calculate the VAF of the mutation in the normal
        normVaf = altCount / (refCount + altCount)
        # Based on the VAF, is this variant germline?
        if normVaf > minVafThreshold:  # Germline
            return True
        elif normVaf < maxVafThreshold:  # Somatic
            return False

        # Perform a fisher's exact test to determine the mutation occurs more frequently in the tumour
        pVal = fisher_exact(refCount, altCount, altPosition.strandCounts[refAllele], altPosition.strandCounts[altAllele]).right_tail
        if pVal < pValThresh:
            return False  # Somatic
        else:
            return True  # Germline

    def _checkIfGermlineIndel(self, chrom, position, altPosition, pValThresh=0.05, minReadThreshold=6, sequenceWindow=200, minVafThreshold=0.33, maxVafThreshold=0.08):
        """
        Determines if a given indel is a germline mutation

        Performs local realignment of the reads in the matched normal which could overlap this indel. Creates a haplotype
        for the reference and alternate alleles, then maps all reads against those alleles to determine the best mappings

        :param chrom: A string indicating the name of the contig the read is mapped against
        :param position: An integer indicating the position of the variant in the specified contig
        :param altPosition: A position() object containing the variant position and aggregated stats
        :param pValThresh: A float indicating the p-value threshold. If the fisher's pvalue is above this threshold, the variant is flagged as germline
        :param minReadThreshold: An int indicating the minimum normal depth required to confidently assign a variant as germline or not
        :param minVafThreshold: A float. If the VAF in the normal sample is higher than this threshold, the variant is flagged as germline
        :param maxVafThreshold: A flat. If the normal sample VAF is lower than this threshold, the variant is flagged as somatic
        :return: A boolean indicating if the variant is germline, or None, if there is insufficient depth
        """

        class discountPysamRead:
            # This is a placeholder for a read, so we can provide a query_sequence
            def __init__(self, seq, position, length):
                self.query_sequence = seq
                self.reference_start = position
                self.query_alignment_length = length

        # Since reads in the normal could be soft-clipped at this indel, include the soft clipping offset
        pileupStart = position - self._softClipUntilIndel
        if pileupStart < 0:
            pileupStart = 0
        indelSize = len(altPosition.alt) - len(altPosition.ref)
        pileupEnd = position - indelSize if indelSize < 0 else position + 1
        pileupEnd += self._softClipUntilIndel
        windowStart = 0
        first = True
        windowEnd = 0
        positions = self._normalFile.pileup(contig=chrom, start=pileupStart, stop=pileupEnd, maxDepth=20000)
        # Parse all reads which overlap this position
        read1Reads = {}
        read2Reads = {}
        for pileupPos in positions:
            if first:
                windowStart = pileupPos.pos
                first = False
            windowEnd = pileupPos.pos
            # Obtain all the reads which overlap this position
            # If we have already included this read, don't double count it
            for read in pileupPos.pileups:
                read = read.alignment
                if read.is_read1:
                    read1Reads[read.query_name] = read
                elif read.is_read2:
                    read2Reads[read.query_name] = read

        normalReads = tuple(read1Reads.values()) + tuple(read2Reads.values())

        # If there is too little coverage at this locus, then we can't accurately determine if this mutation is germline
        if len(normalReads) < minReadThreshold:
            return None

        windowStart -= sequenceWindow
        windowEnd += sequenceWindow
        # Set the reference buffer correctly so the correct haplotypes are generated
        self._refWindow = self._refGenome[self._chrom][windowStart - self._realignBuffer:windowEnd + self._realignBuffer].seq
        self._refStart = windowStart - self._realignBuffer

        # If this is a deletion, remove bases from the sequence
        if indelSize < 0:
            altSeq = self._refGenome[chrom][windowStart:position].seq + self._refGenome[chrom][position - indelSize: windowEnd - indelSize].seq
        else:
            # Add the insertion into the sequence
            altSeq = self._refGenome[chrom][windowStart-1:position-1].seq + altPosition.alt + \
                     self._refGenome[chrom][position: windowEnd].seq

        simRead = discountPysamRead(altSeq, windowStart, windowEnd - windowStart)
        # Create the relevant haplotypes
        refHap, altHap = self.generateHaplotypes([simRead], softclippingBuffer=0)

        # Align all reads in the normal against these haplotypes, and determine which reads support which haplotype
        refReads = []
        altReads = []
        for read in normalReads:
            refAlign = refHap.alignmentStructure(read.query_sequence)
            altAlign = altHap.alignmentStructure(read.query_sequence)
            # If this read maps equally well to both the reference and alternate haplotype, then it probably doesn't
            # overlap the indel, and we should ignore it
            if altAlign.optimal_alignment_score > refAlign.optimal_alignment_score:
                altReads.append(read.query_name)
            elif altAlign.optimal_alignment_score < refAlign.optimal_alignment_score:
                refReads.append(read.query_name)

        # Avoid double-counting any read pairs
        # In the case that one read supports the reference allele, while the other supports the alternate allele, keep
        # the read which supports the alternate allele
        altReads = set(altReads)
        refReads = set(x for x in refReads if x not in altReads)

        altCount = len(altReads)
        refCount = len(refReads)

        # Check for sufficient coverage again
        if altCount + refCount < minReadThreshold:
            return None

        # Calculate normal VAF
        normVaf = altCount / (refCount + altCount)
        # Based on the VAF, can we assign this mutation as somatic or germline?
        # We set a VAF threshold, as very low VAF variants in the tumour will be flagged as germline even if there
        # are no reads supporting the alternate allele in the tumour
        if normVaf > minVafThreshold:
            return True  # Germline
        elif normVaf < maxVafThreshold:
            return False  # Somatic

        # Perform a fisher's exact test to determine the mutation occurs more frequently in the tumour
        pVal = fisher_exact(refCount, altCount, altPosition.strandCounts["REF"], altPosition.strandCounts["ALT"]).right_tail
        if pVal < pValThresh:
            return False  # Somatic
        else:
            return True  # Germline

    def processRead(self, read, rCigar):
        """
        Add all positions covered by this read to the pileup

        :param read: A pysam.AlignedSegment() object
        :param rCigar: A list containing cigar operators
        :param mismatchMax: The lead length multiplied this is the maximum number of mismatched permitted before the read is flagged as supporting an indel
        """

        if self._oBAM:
            self._oBAM.write(read)
        # Is this read valid? Check the naming scheme of the read, and identify each feature stored in the name
        try:
            nameElements = read.query_name.split(":")
            # Barcode is element 1. This is no longer needed, so ignore it
            # Parental strand is element 2.
            assert nameElements[1] == "+" or nameElements[1] == "-"
            rPosParent = nameElements[1] == "+"
            # Family size is element 3
            rFSize = int(nameElements[2])
            # Counter is 4. This is used to identify duplexes
            counter = nameElements[3]

        except (TypeError, IndexError, AssertionError):
            sys.stderr.write(
                "ERROR: The names of the reads in this BAM file are not consistent with those generated by ProDuSe Collapse\n")
            sys.stderr.write(
                "We expected something line \'TAATGCATCTTGATTTGGTTGCGAGTTGCAAT:+:207:0\', and instead saw %s\n" % (
                    read.query_name))
            exit(1)

        # Obtain generic read characteristics
        rMappingQual = read.mapping_quality

        # All positions covered by this read
        positions = []

        # The number of variants supported by this read
        numMismatch = 0

        # Iterate through all bases and cigar indexes
        cigarIndex = 0  # Keep this seperate to account for deletions
        readIterator = zip(read.query_sequence, read.query_qualities)
        position = read.reference_start

        # Indel attributes
        indelPos = 0
        indelType = None
        indelQual = []
        lastCigar = 0
        indelSeq = ""
        indelEndDist = []
        indelRef = ""  # Ref base where the indel starts

        try:
            while True:

                cigar = rCigar[cigarIndex]
                cigarIndex += 1

                if cigar == 0:  # Matched base
                    # This is the most common case by far
                    base, qual = next(readIterator)
                    # Calculate the distance from the end of the read for this variant
                    dist2Start = position - read.reference_start
                    dist2End = read.reference_end - position
                    if dist2Start < dist2End:
                        distFromEnd = dist2Start
                    else:
                        distFromEnd = dist2End

                    if lastCigar != 1 and lastCigar != 2:
                        # If this position overlaps a candidate indel, we also need to add this position to those indels
                        # But don't do this if the last position was an indel, as we do not want to double-count this read
                        try:
                            if position in self.rawIndels[self._chrom]:
                                for indel in self.rawIndels[self._chrom][position ].values():
                                    indel.add(base, qual, rPosParent, rFSize, read.is_reverse, rMappingQual,
                                              distFromEnd,
                                              counter,
                                              read.query_name)
                                    positions.append(indel)
                        except KeyError:
                            pass

                    # If the previous event was an indel, we need to process that, as we have reached the end of the indel
                    if lastCigar != 0:
                        # Add the indel as a candidate indel
                        try:
                            self.rawIndels[self._chrom][indelPos][indelType].add(indelSeq, indelQual, rPosParent,
                                                                                 rFSize, read.is_reverse, rMappingQual,
                                                                                 indelEndDist, counter, read.query_name,
                                                                                 isAlt=True)
                            positions.append(self.rawIndels[self._chrom][indelPos][indelType])
                            if self.rawIndels[self._chrom][indelPos][indelType].ref is None:
                                self.rawIndels[self._chrom][indelPos][indelType].ref = indelRef
                        except KeyError:  # No indel was discovered here during haplotype reassembly. This is likely (?) an artifact?
                            pass
                        # Clear the indel buffer
                        # It will actually be overwritten, so we don't need to do anything except reset the quality list
                        indelQual = []
                        lastCigar = 0

                    # Obtain the pileup object corresponding to this position
                    # If this position has not been covered before, we need to generate a new pileup at this position
                    if position not in self.pileup[self._chrom]:
                        refBase = self._refWindow[position - self._refStart]
                        self.pileup[self._chrom][position] = Position(refBase)

                    pileupPos = self.pileup[self._chrom][position]
                    position += 1

                    # Ignore bases that are "N"
                    if base == "N":
                        continue
                    isAlt = pileupPos.add(base, qual, rFSize, rPosParent, read.is_reverse, distFromEnd, rMappingQual,
                                          counter,
                                          read.query_name)

                    # Check to see if this position now has evidence of a variant
                    if isAlt:
                        numMismatch += 1
                    positions.append(pileupPos)

                elif cigar == 4:  # Soft clipped. Ignore this
                    next(readIterator)
                    continue
                elif cigar == 2:  # A deletion

                    # Was the last cigar operator a normal match?
                    if lastCigar == 0:
                        # Create a new indel
                        lastCigar = 2
                        indelPos = position
                        indelSeq = self._refWindow[position - self._refStart]
                        indelType = "D"
                        indelQual.append(None)
                        indelEndDist = min(abs(position - read.reference_start), abs(read.reference_end - position))
                        indelRef = self._refWindow[position - self._refStart]
                    elif lastCigar == 2:
                        # Continue the previous deletion
                        indelSeq += self._refWindow[position - self._refStart]
                        indelQual.append(None)
                    elif lastCigar == 1:
                        # This is an alignment artifact caused by skbio's SSW algorithm. Skip it
                        indelPos = 0
                        indelType = None
                        indelQual = []
                        lastCigar = 0
                        indelSeq = ""
                        indelEndDist = 0
                    position += 1
                    continue
                elif cigar == 1:  # An insertion.
                    base, qual = next(readIterator)
                    distFromEnd = min(abs(position - read.reference_start), abs(read.reference_end - position))
                    # Was the last operator a normal match?
                    if lastCigar == 0:
                        # Create a new indel
                        lastCigar = 1
                        indelPos = position
                        indelSeq = base
                        indelType = "I"
                        indelQual.append(qual)
                        indelEndDist = distFromEnd
                        indelRef = self._refWindow[position - self._refStart - 1]
                    elif lastCigar == 1:  # We are continuing the previous insertion
                        indelSeq += base
                        indelQual.append(qual)
                    elif lastCigar == 2:
                        # This is an alignment artifact caused by local realignment. Skip it
                        indelPos = 0
                        indelType = None
                        indelQual = []
                        lastCigar = 0
                        indelSeq = ""
                        indelEndDist = 0
                    continue

        except IndexError:
            pass
        except StopIteration as e:
            raise pe.MalformedReadException("Read \'%s\' appears to be malformed" % read.query_name) from e

        # Store the number of total mismatches for this read at each variant
        for pos in positions:
            pos.addMismatchNum(numMismatch)

    def generateHaplotypes(self, indelReads, softclippingBuffer=100):
        """
        Generate a reference using all reads which contain an indel

        :return: A tuple of reference haplotypes, or None if there are no alternative haplotypes
        """

        # If no reads in this window support an indel, then there are no alternative haplotypes
        if not indelReads:
            return None

        if not self._bufferedReads or indelReads[0].reference_start < self._bufferedReads[0].reference_start:
            self._bufferPos = indelReads[0].reference_start - softclippingBuffer
            if self._bufferPos < 0:
                self._bufferPos = 0
            refStart = self._bufferPos - self._refStart
        else:
            self._bufferPos = self._bufferedReads[0].reference_start - softclippingBuffer
            if self._bufferPos < 0:
                self._bufferPos = 0
            refStart = self._bufferPos - self._refStart

        # What is the general length of each read?
        readLength = indelReads[-1].query_alignment_length
        refEnd = refStart + self._realignBuffer + softclippingBuffer + readLength

        referenceSeq = self._refWindow[refStart:refEnd]
        refHaplotype = Haplotype(referenceSeq)

        candidateHaplotypes = {referenceSeq: refHaplotype}

        # First pass: Align all reads against the reference haplotype using SW-alignment to identify any differences
        # between this read and the reference
        for read in indelReads:

            # Align this read against the reference
            altAlignment = refHaplotype.alignmentStructure(read.query_sequence)
            altHaplotype = []
            # Generate a haplotype for this alignment
            # Also determine the location of any alternative events
            i = altAlignment.query_begin - 1
            eventSize = 0
            eventType = None
            eventLoc = sortedcontainers.SortedDict()
            for b1, b2 in zip(altAlignment.aligned_query_sequence, altAlignment.aligned_target_sequence):
                i += 1
                if b2 == "-" and b1 == "-":
                    i -= 1
                elif b2 == "-": # Deletion
                    eventSize += 1
                    eventType = "D"
                elif b1 == "-": # Insertion
                    altHaplotype.append(b2)
                    eventType = "I"
                    eventSize += 1
                else:  # Match/mismatch. Use the reference base
                    altHaplotype.append(b1)
                    if eventType:
                        eventLoc[i-eventSize] = [eventSize, eventType]
                        if eventType == "D":
                            i -= eventSize
                        eventType = None
                        eventSize = 0

            altAssembledHaplotype = referenceSeq[:altAlignment.query_begin] + "".join(altHaplotype) + referenceSeq[altAlignment.query_end+1:]
            if altAssembledHaplotype not in candidateHaplotypes:  # This haplotype does not match an existing haplotype perfectly
                # See if this haplotype matches any other haplotype
                altMatchFound = False
                for name, haplotype in candidateHaplotypes.items():
                    alignBtwnHap = haplotype.alignmentStructure(altAssembledHaplotype)
                    if "I" not in alignBtwnHap.cigar and "D" not in alignBtwnHap.cigar:
                        # i.e. there are no indels between these two haplotypes. They are likely the same event
                        # Keep the event which is closer to the reference
                        if haplotype.scoreFromRef is not None and altAlignment.optimal_alignment_score > haplotype.scoreFromRef:
                            candidateHaplotypes[name] = Haplotype(altAssembledHaplotype, eventLoc, altAlignment.optimal_alignment_score)
                        altMatchFound = True
                        break

                if not altMatchFound:
                    # Store this new haplotype
                    candidateHaplotypes[altAssembledHaplotype] = Haplotype(altAssembledHaplotype, eventLoc, altAlignment.optimal_alignment_score)

        return tuple(candidateHaplotypes.values())

    def _cigarFromAlignment(self, align, haplotypeDist):
        """q
        For a given read (query sequence), generate a pysam-style cigar tuple
        Incorperate any differences between this haplotype and the reference as well
        :return:
        """

        cigarLength = align.target_begin
        cigar = [4] * align.target_begin
        length = ""

        for x in align.cigar:
            if x == "M":
                l = int(length)
                cigar.extend([0] * l)
                cigarLength += l
                length = ""
            elif x == "D":
                l = int(length)
                cigar.extend([1] * int(length))  # Switch operators because the output cigar will refer to the read, not the reference
                # Since insertions do not consume reference positions, we need to start an offset here
                cigarLength += l
                length = ""
            elif x == "I":
                cigar.extend([2] * int(length))
                length = ""
            else:
                length += x

        # Add trailing soft clipping
        if cigarLength < len(align.target_sequence):
            cigar.extend([4] * (len(align.target_sequence) - cigarLength))


        readStartOffset = 0
        # Add the event that is represented by this haplotype
        for index in reversed(haplotypeDist):
            event = haplotypeDist[index]
            try:
                index = index - align.query_begin + align.target_begin
                length, type = event

                # Ignore events at the start or end of the read (including soft-clipping
                # If these events actually existed, they would be flanked by mapped bases in the cigar
                if index <= 0:
                    if type == "D":  # This is a deletion which consumes reference positions
                        readStartOffset += length
                    elif type == "I":  # This is an insertion which does not consume reference positions
                        readStartOffset -= length
                    continue
                elif cigar[index + 1] == 4:
                    continue

                if cigar[index] != 4:  # i.e. this base is not soft clipped
                    if type == "D":  # This event is a deletion. We need to insert it into the cigar
                        # Since skbio's implementation of SSW isn't the most reliable thing, in some cases inferior
                        # alignments can be given a higher score (somehow)
                        # If this is the case, correct the cigar
                        for i in range(0, length):
                            if cigar[i + index] == 1:
                                cigar[i + index] = 0
                            else:
                                cigar.insert(i + index, 2)
                    else:  # This event is an insertion. We need to replace the existing cigar operators
                        for i in range(0, length):
                            # Once again, account for skbio being weird
                            if cigar[i+index] == 2:
                                del cigar[i + index]
                            elif cigar[i+index] == 1:  # i.e. There is already an insertion here.
                                # We need to find a mapped base. The location of this shouldn't matter
                                try:
                                    iIndex = index
                                    while True:
                                        iIndex += 1
                                        x = cigar[i+iIndex]
                                        if x == 0:  # A mapped base! Just make sure the next base is not soft-clipped
                                            if cigar[i+iIndex+1] == 4:  # It's soft-clipped. Try the other end
                                                raise IndexError
                                            cigar[i + iIndex] = 1
                                            break
                                        elif x == 4:  # We have ended up in soft-clipping
                                            raise IndexError
                                except IndexError:
                                    # We reached the end of the read when looking in the forward direction. Try the reverse
                                    try:
                                        iIndex = index
                                        while True:
                                            iIndex -= 1
                                            x = cigar[i + iIndex]
                                            if x == 0:  # A mapped base! Just make sure the next base is not soft-clipped
                                                if cigar[i + iIndex - 1] == 4:  # It's soft-clipped. Try the other end
                                                    raise IndexError
                                                cigar[i + iIndex] = 1
                                                break
                                            elif x == 4:  # We have ended up in soft-clipping
                                                raise IndexError
                                    except IndexError:
                                        pass  # This read will be malformed, but there is no way to generate a valid alignment from it
                                        # This SHOULD never occur
                            else:
                                cigar[i+index] = 1
            except IndexError:
                continue

        # Convert the cigar list to a pysam-style list of tuples
        length = 0
        oldOp = None
        pysamCigar = []
        for op in cigar:
            if op != oldOp:
                if oldOp is not None:
                    pysamCigar.append((oldOp, length))
                length = 1
                oldOp = op
            else:
                length += 1
        pysamCigar.append((op, length))

        return pysamCigar, cigar, readStartOffset

    def _realign(self, haplotypes, endPoint=None, indelWindow = 200):
        """
        Realign all reads stored in the buffer against these haplotypes, and add the realigned reads to the pileip

        :param haplotypes: A list of Haplotype objects, which will be aligned against
        :param endPoint: An integer. If the reference_start position of the read is greater than this, the read will not be processed
        :param indelWindow: An integer specifying a window offset to be applied when processing indel reads.

        To account for the situtation whereby a haplotype may not be generated because all reads that supported that
        haplotype were processed in an earlier window, don't process all indel-supporting reads in the earlier window
        :return:
        """

        def readToIndel(pos, event, length, read):

            # Add position
            if pos not in self.rawIndels[read.reference_name]:
                self.rawIndels[read.reference_name][pos] = {}
            # Add event
            if event not in self.rawIndels[read.reference_name][pos]:
                self.rawIndels[read.reference_name][pos][event] = IndelPos(event, length)

        i = 0
        for read in self._indelReads:

            if endPoint is not None and read.reference_start > endPoint - indelWindow:  # i.e. we have processed all the reads we need to so far.
                break
            maxScore = 0
            maxAlignment = None
            maxHap = None
            # Find the best haplotype mapping for this read
            for hap in haplotypes:
                align = hap.alignmentStructure(read.query_sequence)
                if align.optimal_alignment_score > maxScore:
                    maxScore = align.optimal_alignment_score  # TODO: Add score offset here
                    maxAlignment = align
                    maxHap = hap

            # If the best haplotype mapping is the reference, keep the original alignment
            if maxHap.eventFromRef:
                # Use the alt alignment
                pysamCigar, listCigar, startOffset = self._cigarFromAlignment(maxAlignment, maxHap.eventFromRef)
                maxHap.support.append(read)

                # Update read statistics
                read.reference_start = maxAlignment.query_begin + self._bufferPos + startOffset
                read.cigartuples = pysamCigar

                # Create an indel position for these events, if it does not exist already
                for pos, event in maxHap.eventFromRef.items():

                    length, indelType = event
                    # Where is this event?
                    # We need to subtract 1, since we will be referencing the indel from the proceeding mapped base
                    pos = self._bufferPos + pos

                    readToIndel(pos, indelType, length, read)

            else:
                # For reads in which we are using the original alignment, generate a cigar list for easy processing
                listCigar = self.cigarToTuple(read.cigartuples)
            # Add this realigned read to the pileup
            self.processRead(read, listCigar)

            i += 1

        # Remove all reads that have been processed from the indel buffer
        self._indelReads = self._indelReads[i:]

        i = 0
        # Now, process all normal reads
        for read in self._bufferedReads:
            # Don't process all reads, in case reads near the end of the buffer map to a haplotype that doesn't exist yet
            if endPoint is not None and read.reference_start > endPoint:
                break
            maxScore = 0
            maxAlignment = None
            maxHap = None
            for hap in haplotypes:
                align = hap.alignmentStructure(read.query_sequence)
                if align.optimal_alignment_score > maxScore:
                    maxScore = align.optimal_alignment_score
                    # Bias in favor of the reference alignment, since skbio's implementation of SSW does strange things
                    if maxHap is None:
                        maxScore += 2
                    maxAlignment = align
                    maxHap = hap

            # If the best haplotype mapping is the reference, keep the original alignment
            if maxHap.eventFromRef:
                pysamCigar, listCigar, startOffset = self._cigarFromAlignment(maxAlignment, maxHap.eventFromRef)
                maxHap.support.append(read)

                # Update read statistics
                read.reference_start = maxAlignment.query_begin + self._bufferPos + startOffset
                read.cigartuples = pysamCigar

                # Create an indel position for these events, if it does not exist already
                for pos, event in maxHap.eventFromRef.items():
                    length, indelType = event
                    # Where is this event?
                    pos = self._bufferPos + pos

                    readToIndel(pos, indelType, length, read)
            else:
                # For reads in which we are using the original alignment, generate a cigar list for easy processing
                listCigar = self.cigarToTuple(read.cigartuples)
            # Add this realigned read to the pileup
            self.processRead(read, listCigar)
            i += 1
        # Remove all realigned reads from the buffer
        self._bufferedReads = self._bufferedReads[i:]

    def cigarToTuple(self, cigarTuples):
        """
        Converts a pysam-style cigar tuples into a regular tuple, where 1 operator = 1 base

        :param cigarTuples: A list of tuples
        :return: A tuple listing expanded cigar operators
        """

        cigarList = []
        for cigarElement in cigarTuples:
            cigOp = cigarElement[0]
            # Ignore soft-clipped bases
            cigarList.extend([cigOp] * cigarElement[1])

        return tuple(cigarList)

    def generatePileup(self, chrom=None, readProcessBuffer=300):
        """

        :param chrom: A string indicating which contig to process. If None, all contigs will be processed
        :param readProcessBuffer:An int specifying how large (in terms of genomic coordinates) the read buffer should be before reads are added to the pileup
        Since this is (probably) cfDNA, fragment sizes are smaller and the buffer could be made smaller, but I am playing it safe here
        :return:
        """
        def processWindow(pos, coordinate, coordinateEnd, candidateVar):

            # Ignore non-variant positions
            if pos.altDepth >= self._minAltDepth:

                # Is this position within the capture space?
                if self._captureSpace is not None:
                    if self._chrom not in self._captureSpace or bisect.bisect_right(self._captureSpace[self._chrom],
                                                                                    coordinate) % 2 != 1:
                        return  # This position does not fall within the capture space. Ignore it

                # If there is a variant, we need to obtain some general characteristics to be used for filtering
                # First, store the sequence of the surrounding bases. This will be used to identify errors
                # due to homopolymers
                leftFlank = []
                rightFlank = []

                windowPos = coordinate - self._refStart  # Where, relative to the reference window, are we located?
                windowPosEnd = coordinateEnd - self._refStart
                if windowPos < 0 :  # i.e. we are processing bases outside the current window position
                    # This is ugly, but we need to parse the reference sequence directly from the FASTA file
                    # Thus, we need to parse the sequence here
                    for j in range(coordinate - self._homopolymerWindow, coordinate):
                        leftFlank.append(self._refGenome[self._chrom][j].seq)
                    for j in range(coordinateEnd + 1, coordinateEnd + self._homopolymerWindow + 1):
                        rightFlank.append(self._refGenome[self._chrom][j].seq)
                else:
                    # Obtain the sequence from the left flank
                    for j in range(windowPos - self._homopolymerWindow, windowPos):
                        try:
                            leftFlank.append(self._refWindow[j])
                        except IndexError:  # In the rare case where we overflow outside the window
                            leftFlank.append(self._refGenome[self._chrom][j + self._refStart].seq)
                    # Obtain the right flank sequence
                    for j in range(windowPosEnd + 1, windowPosEnd + self._homopolymerWindow + 1):
                        try:
                            rightFlank.append(self._refWindow[j])
                        except IndexError:  # We overflowed outside the window. This case is a lot more likely for the right window
                            rightFlank.append(self._refGenome[self._chrom][j + self._refStart].seq)


                # Next, examine how "noisy" the neighbouring region is (i.e. how many candidiate variants there are)
                # Flag all nearby candidate variants. To avoid flagging the same variant twice, add this
                # candidate variant position to all variants behind it, and vise-versa
                if self._chrom not in candidateVar:
                    candidateVar[self._chrom] = sortedcontainers.SortedDict()

                nearbyVar = []
                for j in candidateVar[self._chrom].irange(coordinate - self._noiseWindow):
                    candidateVar[self._chrom][j].nearbyVar.append(pos)
                    nearbyVar.append(candidateVar[self._chrom][j])

                # Finally, save all these stats we have just generated inside the variant, to be examined during
                # filtering
                pos.processVariant(leftFlank, rightFlank, nearbyVar, self._noiseWindow)
                candidateVar[self._chrom][coordinate] = pos

        def realignReadsAndPileup(windowEndPoint=None):
            """
            Add all reads stored in the buffer to the pileup

            First, generate all possible haplotypes for this region using the reads which contain an indel
            Next, realign ALL buffered reads in this region against this indel and the normal reference
            Finally, process all reads, and generate a pileup from the positions
            :return:
            """
            # First, generate the haplotypes using reads with INDELs
            haplotypes = self.generateHaplotypes(self._indelReads)
            # Map all reads in the buffer against these haplotypes
            if haplotypes is not None:
                self._realign(haplotypes, windowEndPoint)

                # Next, generate an indel position object for each haplotype
                for hap in haplotypes:
                    if len(hap.support) < self._minAltDepth:
                        continue

            # If there are no alternative haplotypes, just add all reads to the pileup
            else:
                i = 0
                for read in self._indelReads:
                    if windowEndPoint is not None and read.reference_start > windowEndPoint:  # i.e. we have processed all the reads we need to so far.
                        break
                    cigarTuple = self.cigarToTuple(read.cigartuples)
                    self.processRead(read, cigarTuple)
                    i += 1
                self._indelReads = self._indelReads[i:]
                i = 0
                for read in self._bufferedReads:
                    if windowEndPoint is not None and read.reference_start > windowEndPoint:  # i.e. we have processed all the reads we need to so far.
                        break
                    cigarTuple = self.cigarToTuple(read.cigartuples)
                    self.processRead(read, cigarTuple)
                    i += 1
                self._bufferedReads = self._bufferedReads[i:]
            # Finally, move the alignment buffer up to reflect the window that was just processed
            self._bufferPos = windowEndPoint

        lastPos = -1

        # Sanity check that the realignment buffer is larger than the actual processing buffer, to ensure that
        # we are confident of the haplotypes that have been generated
        if self._realignBuffer < readProcessBuffer:
            raise AttributeError("Based on the current buffer sizes, reads will be added to the pileup before they are properly realigned")

        for read in self._inFile.fetch(contig=chrom):

            # Prior to adding each base in this read onto the pileup, we need to analyze the current read,
            # and obtain general characteristics (family size, mapping quality etc)

            # Print out status messages
            self.readCount += 1
            if self.readCount % 100000 == 0:
                if chrom is not None:
                    sys.stderr.write(
                    "\t".join([self._printPrefix, time.strftime('%X'), chrom + ":Reads Processed:%s\n" % self.readCount]))
                else:                    sys.stderr.write(
                    "\t".join([self._printPrefix, time.strftime('%X'), "Reads Processed:%s\n" % self.readCount]))

            # Ignore unmapped reads
            if read.is_unmapped:
                continue

            # Ignore reads without a cigar string, as they are malformed or empty
            if not read.cigartuples:
                continue

            # Ignore reads that are secondary or supplemental alignments
            if read.is_secondary or read.is_supplementary:
                continue

            # Ignore reads with a mapping quality of 0
            if read.mapping_quality == 0:
                continue

            # Have we advanced positions? If so, we may need to change the loaded reference window, and process
            # positions which no more reads will be mapped to
            if read.reference_name != self._chrom or read.reference_start != lastPos:

                # Sanity check that the input BAM is sorted
                if read.reference_name == self._chrom and read.reference_start < lastPos:
                    raise pe.UnsortedInputException("Input SAM/BAM file does not appear to be sorted")

                # Process positions at which no more reads are expected to map to (assuming the input is actually sorted)
                # We do this to drastically reduce the memory footprint of the pileup, since we don't care about non-variant
                # positions very much
                try:
                    if read.reference_name != self._chrom and self._chrom is not None:

                        # If we have switched chromosomes, process all variants stored on the previous chromosome
                        # First, assemble haplotypes for the reads aligned in this region
                        realignReadsAndPileup(None)
                        assert not self._bufferedReads
                        self._bufferPos = read.reference_start

                        # Process SNVs
                        posToProcess = self.pileup[self._chrom].keys()
                        for coordinate in posToProcess:
                            processWindow(self.pileup[self._chrom][coordinate], coordinate, coordinate, self.candidateVar)
                        del self.pileup[self._chrom]

                        # Process indels
                        # We only need to do this once per chromosome, as there are a lot less indels than possible SNVs
                        posToProcess = self.rawIndels[self._chrom].keys()
                        for coordinate in posToProcess:

                            for type, indelObj in self.rawIndels[self._chrom][coordinate].items():
                                # What is the end position of this indel?
                                # If it's an insertion, the start and end points are the same
                                if type == "I":
                                    endPos = coordinate
                                else: # If it's a deletion, add the length of the deletion
                                    endPos = coordinate + indelObj.length  # Add the length of the deletion at this position

                                processWindow(indelObj, coordinate, endPos, self.candidateIndels)
                        del self.rawIndels[self._chrom]

                        # Since this is the first read we are processing from this chromosome, we need to create the dictionary which
                        # will store all bases on this chromosome
                        self.pileup[read.reference_name] = sortedcontainers.SortedDict()
                        self.rawIndels[read.reference_name] = sortedcontainers.SortedDict()

                    # Otherwise, check to see if previous positions fall outside the buffer window. If so, process them
                    elif self._chrom is not None and read.reference_start > self._bufferPos + self._realignBuffer:
                        realignReadsAndPileup(read.reference_start - readProcessBuffer)
                        posToProcess = tuple(self.pileup[self._chrom].irange(minimum=None,
                                                                             maximum=read.reference_start - self._realignBuffer))
                        for coordinate in posToProcess:
                            processWindow(self.pileup[self._chrom][coordinate], coordinate, coordinate, self.candidateVar)
                        for coordinate in posToProcess:
                            del self.pileup[self._chrom][coordinate]
                    elif self._chrom is None:
                        self.pileup[read.reference_name] = sortedcontainers.SortedDict()
                        self.rawIndels[read.reference_name] = sortedcontainers.SortedDict()
                        self._bufferPos = read.reference_start

                except KeyError as e:  # i.e. All reads which mapped to this contig or previous positions failed QC. There are no positions to process
                    pass

                # Load the current reference window and it's position
                # If we have switched contigs, or moved outside the window that is currently buffered,
                # we need to re-buffer the reference
                readWindowStart = read.reference_start - self._realignBuffer
                if readWindowStart < 0:
                    readWindowStart = 0
                readWindowEnd = read.reference_end + self._realignBuffer

                if read.reference_name != self._chrom or readWindowStart < self._refStart or readWindowEnd > refEnd:
                    # To reduce the number of times this needs to be performed, obtain a pretty large buffer
                    self._refStart = readWindowStart - self._realignBuffer - self._homopolymerWindow
                    if self._refStart < 0:
                        self._refStart = 0
                    refEnd = readWindowStart + 8000
                    self._chrom = read.reference_name
                    self._refWindow = self._refGenome[self._chrom][self._refStart:refEnd].seq

                lastPos = read.reference_start

            # Does this read overlap the capture space? If not, ignore it
            # Note that this will keep reads which partially overlap the capture space
            # We will perform a more precise restricting of the regions which overlap the capture space later
            if not self._insideCaptureSpace(read):
                continue

            # If this read supports an INDEL, we need to set it aside and generate a haplotype for this read
            if "I" in read.cigarstring or "D" in read.cigarstring:
                self._indelReads.append(read)
            else:
                self._bufferedReads.append(read)  # No indels

        # Now that we have finished processing all reads in the input file, process all remaining positions
        try:
            realignReadsAndPileup(None)
            posToProcess = self.pileup[self._chrom].keys()
            for coordinate in posToProcess:
                processWindow(self.pileup[self._chrom][coordinate], coordinate, coordinate, self.candidateVar)
            # Process indels
            posToProcess = self.rawIndels[self._chrom].keys()
            for coordinate in posToProcess:
                for type, indelObj in self.rawIndels[self._chrom][coordinate].items():
                    # What is the end position of this indel?
                    # If it's an insertion, the start and end points are the same
                    if type == "I":
                        endPos = coordinate
                    else:  # If it's a deletion, add the length of the deletion
                        endPos = coordinate + indelObj.length  # Add the length of the deletion at this position
                    processWindow(indelObj, coordinate, endPos, self.candidateIndels)
            del self.rawIndels[self._chrom]
        except KeyError:
            pass


def runCallMultithreaded(callArgs):
    """
    Parallelize the pileup engine by chromosome

    :return:
    """

    pileupArgs = callArgs[0]
    filterArgs = callArgs[1]
    contig = callArgs[2]

    pileup = PileupEngine(*pileupArgs)

    pileup.generatePileup(chrom=contig)

    pileup.filterAndWriteVariants(*filterArgs)


def isValidFile(file, parser):
    """
    Checks to ensure a specified file exists, and throws an error if it does not
    :param file: A string containing a filepath
    :param parser: An argparse.argumentparser() object
    :return: The variable "file", if the file exists
    :raises: parser.error() if the file does not exist
    """
    if os.path.exists(file):
        return file
    else:
        raise parser.error("Unable to locate \'%s\': No such file or directory" % file)


def validateArgs(args):
    """
    Checks that the specified set of arguments are valid
    :param args: A dictionary listing {argument: parameter}
    :return: A dictionary listing {argument: parameter} that have been validated
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

    # Reparse the sanitized arguments
    # List requirements here
    parser = argparse.ArgumentParser(description="Identifies and calls variants")
    parser.add_argument("-c", "--config", metavar="INI", type=lambda x: isValidFile(x, parser),
                        help="An optional configuration file which can provide one or more arguments")
    parser.add_argument("-i", "--input", metavar="BAM", required=True, type=lambda x: isValidFile(x, parser),
                        help="Input post-collapse or post-clipoverlap BAM file")
    parser.add_argument("-n", "--normal", metavar="BAM", type=lambda x: isValidFile(x, parser),
                        help="Optional matched normal BAM file, for removing germline variants")
    parser.add_argument("-o", "--output", metavar="VCF", required=True, help="Output VCF file")
    parser.add_argument("-u", "--unfiltered", metavar="VCF",
                        help="An additional output VCF file, which lists all candidate variants, including those that failed filters")
    parser.add_argument("-r", "--reference", metavar="FASTA", required=True, type=lambda x: isValidFile(x, parser),
                        help="Reference Genome, in FASTA format")
    parser.add_argument("-t", "--targets", metavar="BED", type=lambda x: isValidFile(x, parser),
                        help="A BED file containing regions in which to restrict variant calling")
    parser.add_argument("-f", "--filter", metavar="PICKLE", required=True, type=lambda x: isValidFile(x, parser),
                        help="A pickle file containing the trained filter. Can be generated using \'produse train\'")
    parser.add_argument("-j", "--jobs", metavar="INT", type=int, default=1,
                        help="How many chromosomes to process simultaneously")
    parser.add_argument("--threshold", metavar="FLOAT", type=float, default=0.65,
                        help="Filtering threshold (lower=more lenient) [Default: 0.65]")
    parser.add_argument("--repeat_count_threshold", metavar="INT", type=int, default=4,
                        help="If an indel is an expansion or contraction of a nearby repeat which occurs at least this many times, filter it [Default: 4]")
    parser.add_argument("--duplex_support_only", action="store_true", help="Only output variants with duplex support")
    parser.add_argument("--min_alt_depth", metavar="INT", type=int, default=3,
                        help="Minimum number of reads required to even consider an alternate allele as possibly real [Default: 3]")
    parser.add_argument("--realigned_BAM", metavar="BAM", help="Optional output BAM/SAM file for realigned reads")
    validatedArgs = parser.parse_args(listArgs)

    # Sanity check
    if 1 < validatedArgs.threshold or 0 > validatedArgs.threshold:
        raise parser.error("\'--threshold\' must be a float between 0 and 1")
    if validatedArgs.repeat_count_threshold < 1 or validatedArgs.repeat_count_threshold > 10:
        raise parser.error("\'--repeat_count_threshold\' must be between 1 and 10")
    return vars(validatedArgs)


parser = argparse.ArgumentParser(description="Identifies and calls variants")
parser.add_argument("-c", "--config", metavar="INI", type=lambda x: isValidFile(x, parser),
                    help="An optional configuration file which can provide one or more arguments")
parser.add_argument("-i", "--input", metavar="BAM", type=lambda x: isValidFile(x, parser),
                    help="Input post-collapse or post-clipoverlap BAM file")

parser.add_argument("-n", "--normal", metavar="BAM", type=lambda x:isValidFile(x, parser),
                    help="Optional matched normal BAM file, for removing germline variants")
parser.add_argument("-o", "--output", metavar="VCF", help="Output VCF file, listing all variants which passed filters")
parser.add_argument("-u", "--unfiltered", metavar="VCF",
                    help="An additional output VCF file, which lists all candidate variants, including those that failed filters")
parser.add_argument("-r", "--reference", metavar="FASTA", type=lambda x: isValidFile(x, parser),
                    help="Reference Genome, in FASTA format")
parser.add_argument("-t", "--targets", metavar="BED", type=lambda x: isValidFile(x, parser),
                    help="A BED file containing regions in which to restrict variant calling")
parser.add_argument("-f", "--filter", metavar="PICKLE", type=lambda x: isValidFile(x, parser),
                    help="A pickle file containing a trained filter. Can be generated using \'produse train\'")
parser.add_argument("-j", "--jobs", metavar="INT", type=int, help="How many chromosomes to process simultaneously")
parser.add_argument("--threshold", metavar="FLOAT", type=float,
                    help="Filtering threshold (lower=more lenient) [Default: 0.65]")
parser.add_argument("--repeat_count_threshold", metavar="INT", type=int, help="If an indel is an expansion or contraction of a nearby repeat which occurs at least this many times, filter it [Default: 4]")
parser.add_argument("--duplex_support_only", action="store_true", help="Only output variants with duplex support")
parser.add_argument("--min_alt_depth", metavar="INT", type=int,
                    help="Minimum number of reads required to even consider an alternate allele as possibly real [Default: 3]")
parser.add_argument("--realigned_BAM", metavar="BAM", help="Optional output BAM/SAM file for realigned reads")


def main(args=None, sysStdin=None, printPrefix="PRODUSE-CALL\t"):
    if args is None:
        if sysStdin is None:
            args = parser.parse_args()
        else:
            args = parser.parse_args(sysStdin)
        args = vars(args)

    # If a config file was specified, parse the arguments from
    if args["config"] is not None:
        config = ConfigObj(args["config"])
        try:
            for argument, parameter in config["call"].items():
                if argument in args and not args[
                    argument]:  # Aka this argument is used by call, and a parameter was not provided at the command line
                    args[argument] = parameter
        except KeyError:  # i.e. there is no section named "call" in the config file
            sys.stderr.write(
                "ERROR: Unable to locate a section named \'call\' in the config file \'%s\'\n" % (args["config"]))
            exit(1)


    args = validateArgs(args)

    # Load filter
    try:
        with open(args["filter"], "r+b") as o:
            filterModel = pickle.load(o)
    except KeyError as e:
        # A KeyError will be thrown if loading a RandomForestModel generated with sklearn <0.18.0, while
        # using >0.18.0
        # All we can do is warn the user
        from sklearn import __version__ as currentVer
        raise AttributeError("The filter model \'%s\' was generated by an earlier version of sklearn, " 
                             "and is not compatible with sklearn version %s" % (args["filter"], currentVer)) from e
    # Print status messages
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Starting...\n"]))
    contigNames = Fasta(args["reference"]).records.keys()

    # How much should we multithread?
    if args["jobs"] == 0:  # i.e. use as many threads as possible
        threads = None  # The multiprocessing module will chose the number of threads
    else:
        threads = min(len(contigNames), args["jobs"], os.cpu_count())

    if threads > 1 or threads is None:
        # Multi-threaded. Append output file names with different prefixes to maintain sorted order later on

        # Process arguments
        # This is going to get ugly, since we need to pass seperate arguments to both generatePileup and filterAndWriteVariants
        # Thus, we are going to need a list of list of tuples.....
        # Also, since map doesn't support kwargs, we are going to need to obtain a lot of defaults.

        # Obtain defaults for pileup engine
        arguments = inspect.signature(PileupEngine)
        outBAMIndex = -1
        defaults = []
        i = 0
        for name, value in arguments.parameters.items():
            if name == "inBAM":
                defaults.append(args["input"])
            elif name == "refGenome":
                defaults.append(args["reference"])
            elif name == "targetRegions":
                defaults.append(args["targets"])
            elif name == "minAltDepth":
                defaults.append(args["min_alt_depth"])
            elif name == "oBAM":
                outBAMIndex = i
                defaults.append(None)
            elif name == "printPrefix":
                defaults.append(printPrefix)
            elif name == "normalBAM":
                defaults.append(args["normal"])
            else:
                defaults.append(value.default)
            i += 1

        # Obtain defaults for filterAndWriteVariants
        arguments = inspect.signature(PileupEngine.filterAndWriteVariants)
        outVCFIndex = -1
        outUnfiltIndex = -1
        headerIndex = -1
        vcfDefaults = []
        i = 0
        for name, value in arguments.parameters.items():
            if name == "self":
                i -= 1
            elif name == "outFile":
                outVCFIndex = i
                vcfDefaults.append(None)
            elif name == "unfilteredOut":
                outUnfiltIndex = i
                vcfDefaults.append(None)
            elif name == "writeHeader":
                vcfDefaults.append(False)
                headerIndex = i
            elif name == "filter":
                vcfDefaults.append(filterModel)
            else:
                vcfDefaults.append(value.default)
            i += 1

        # Now, merge these arguments into the multiprocessing list
        # The format will be [[EngineArgs1, filterArgs1, contig1], [EngineArgs2, filterArgs2, contig2]]
        multithreadArgs = tuple((defaults[:], vcfDefaults[:], contig) for contig in contigNames)
        multithreadArgs[0][1][headerIndex] = True
        # Now, add the unique filenames
        i = 0
        bamFiles = []
        oVCFFiles = []
        uVCFFiles = []
        for contig in contigNames:
            # Add output BAM name
            oBAMName = args["realigned_BAM"] + "." + contig
            bamFiles.append(oBAMName)
            multithreadArgs[i][0][outBAMIndex] = oBAMName
            # Add output VCF name
            oVCFName = args["output"] + "." + contig
            oVCFFiles.append(oVCFName)
            if args["unfiltered"] is not None:
                uVCFName = args["unfiltered"] + "." + contig
                uVCFFiles.append(uVCFName)
            else:
                uVCFName = None
            multithreadArgs[i][1][outVCFIndex] = oVCFName
            multithreadArgs[i][1][outUnfiltIndex] = uVCFName
            i += 1
        i = 0

        # Run the jobs
        processPool = multiprocessing.Pool(processes=threads)
        try:
            processPool.map_async(runCallMultithreaded, multithreadArgs)
            processPool.close()
            processPool.join()
        except KeyboardInterrupt as e:
            processPool.terminate()
            processPool.join()
            raise e

        # Finally, merge the output files
        pysam.merge("-f", args["realigned_BAM"], *bamFiles)
        for bFile in bamFiles:
            os.remove(bFile)
        with open(args["output"], "w") as o:
            for oVCFFile in oVCFFiles:
                with open(oVCFFile) as v:
                    for line in v:
                        o.write(line)
                os.remove(oVCFFile)
        if args["unfiltered"] is not None:
            with open(args["unfiltered"], "w") as o:
                for uVCFFile in uVCFFiles:
                    with open(uVCFFile) as v:
                        for line in v:
                            o.write(line)
                    os.remove(uVCFFile)

    else:  # Singe-threaded
        pileup = PileupEngine(args["input"], args["reference"], args["targets"], minAltDepth=args["min_alt_depth"],
                              oBAM=args["realigned_BAM"], normalBAM=args["normal"], printPrefix=printPrefix)
        first = True
        for contig in contigNames:
            # Find candidate variants
            pileup.generatePileup(chrom=contig)

            # Filter variants
            if first:
                pileup.filterAndWriteVariants(args["output"], filterModel, args["unfiltered"], args["threshold"],
                                              args["duplex_support_only"], args["repeat_count_threshold"], writeHeader=True)
                first = False
            else:
                pileup.filterAndWriteVariants(args["output"], filterModel, args["unfiltered"], args["threshold"],
                                              args["duplex_support_only"], args["repeat_count_threshold"])
            pileup.reset()

    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Call Complete\n"]))

if __name__ == "__main__":
    main()
