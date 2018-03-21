#! /usr/bin/env python

import argparse
import os
import pysam
import sys
import sortedcontainers
import pickle
from sklearn.ensemble import RandomForestClassifier
from skbio import alignment
from pyfaidx import Fasta
import time
from fisher import pvalue
from configobj import ConfigObj
import bisect


class Position:
    """
    Stores the allele counts and associated characteristics at this position
    """

    def __init__(self, refBase):
        self.ref = refBase

        # Pre-initialize these lists to a sensible size to improve performance in high-depth regions
        self.alleles = [None] * 5
        self.qualities = [None] * 5
        self.famSizes = [None] * 5
        self.posParent = [None] * 5
        self.mapStrand = [None] * 5
        self.distToEnd = [None] * 5
        self.mismatchNum = [None] * 5
        self.mappingQual = [None] * 5
        self.readNum = [None] * 5
        self.alt = False
        self.altAlleles = set()
        self._altIndex = {}  # Used to remove altAlleles if necessary
        self.depth = 0

        self.skipIndexes = set()

        # Statistics when summarizing variants

    def add(self, base, qual, size, posParent, mapStrand, dist, mapQual, readNum):

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
            self.alleles[self.depth] = base
        self.qualities[self.depth] = qual
        self.famSizes[self.depth] = size
        self.posParent[self.depth] = posParent
        self.mapStrand[self.depth] = mapStrand
        self.distToEnd[self.depth] = dist
        self.mappingQual[self.depth] = mapQual
        self.readNum[self.depth] = readNum

        self.depth += 1

        # Does this base support a variant?
        if self.ref != "N" and self.ref != base:
            self.alt = True
            if base not in self.altAlleles:
                self.altAlleles.add(base)
                self._altIndex[base] = self.depth - 1
            return True
        else:
            return False

    def addPos(self, position):
        """
        Merges the statistics of another position into this one
        :param position: Another Position object
        :return:
        """

        # Remove the Nones from the lists
        self.qualities = self.qualities[:self.depth]
        self.famSizes = self.famSizes[:self.depth]
        self.posParent = self.posParent[:self.depth]
        self.mapStrand = self.mapStrand[:self.depth]
        self.distToEnd = self.distToEnd[:self.depth]
        self.mappingQual = self.mappingQual[:self.depth]
        self.readNum = self.readNum[:self.depth]

        # Add the other position
        self.depth += position.depth
        self.qualities.extend(position.qualities)
        self.famSizes.extend(position.famSizes)
        self.posParent.extend(position.posParent)
        self.mapStrand.extend(position.mapStrand)
        self.distToEnd.extend(position.distToEnd)
        self.mappingQual.extend(position.mappingQual)
        self.readNum.extend(position.readNum)


    def removeLast(self):

        # Remove the last base stored at this position
        i = self.depth - 1
        del self.qualities[i]
        del self.famSizes[i]
        del self.posParent[i]
        del self.mapStrand[i]
        del self.distToEnd[i]
        del self.mappingQual[i]
        del self.readNum[i]
        del self.mismatchNum[i]
        allele = self.alleles[i]

        # Remove this base from the alternate allele list, if necessary
        try:
            if allele != self.ref and self.ref != "N" and self._altIndex[allele] == i:
                self.altAlleles.remove(allele)
                if not len(self.altAlleles):
                    self.alt = False
            del self.alleles[i]
        except KeyError as e:
            raise TypeError(self.readNum) from e
        self.depth -= 1

    def addMismatchNum(self, misMatchNum):
        self.mismatchNum[self.depth -1] = misMatchNum

    def processVariant(self, refWindow, nearbyVariants, nWindow):
        self.refWindow = refWindow
        self.nearbyVar = nearbyVariants
        self.nWindow = nWindow  # How many positions were examined to find nearby variants?

    def collapsePosition(self, minDepth=0, minAltDepth=0, strandBiasThreshold=0, strongMoleculeThreshold=3):
        """
        Merges variant stats at a given position

        :return:
        """

        self.baseCounts = { "DPN":  [0, 0, 0, 0, 0],
                            "DpN":  [0, 0, 0, 0, 0],
                            "DPn":  [0, 0, 0, 0, 0],
                            "Dpn":  [0, 0, 0, 0, 0],
                            "SN" :  [0, 0, 0, 0, 0],
                            "SP" :  [0, 0, 0, 0, 0],
                            "Sp" :  [0, 0, 0, 0, 0],
                            "Sn" :  [0, 0, 0, 0, 0]
                            }
        self.strandCounts = [0.0, 0.0, 0.0, 0.0, 0.0]
        self.posMolec = [0, 0, 0, 0, 0]
        self.negMolec = [0, 0, 0, 0, 0]
        self.disagreeingDuplexes = [0, 0, 0, 0, 0]

        baseToIndex = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4}
        indexToBase = {0: "A", 1: "C", 2: "G", 3: "T", 4: "-"}

        pMapStrand = [0, 0, 0, 0, 0]
        nMapStrand = [0, 0, 0, 0, 0]

        depth = 0.0
        # What is the maximum quality of the alternate allele?
        self.maxQual = [0, 0, 0, 0, 0]

        # We need to identify and flag reads that are in duplex
        # Count the number of times each read number "A counter that is unique to each read, generated by collapse" occurs
        # If it occurs once, the read is a singleton. If it occurs twice, it is a duplex
        # That said, we also need to handle the case whereby clipoverlap was NOT run on the input BAM file. A read pair
        # will have the same number
        # To deal with this, create a dictionary which will store both the parental strand, and index, as well as the index and base for that strand
        self.readParentalStrand = {}
        readBaseIndex = {}

        # Aggregate the stats for each read for filtering
        readStats = {}

        # Since the characteristics of each read is in order (i.e. index 0 corresponds to the first read that overlaps
        # this position for every list, cycle through all stored attributes
        for base, qual, fSize, posParent, map, distToEnd, mismatchNum, mapQual, readID \
                in zip(self.alleles, self.qualities, self.famSizes, self.posParent, self.mapStrand, self.distToEnd, self.mismatchNum, self.mappingQual, self.readNum):

            # If this base is None, then we have started reading into the buffer, and should stop processing bases
            if base is None:
                break

            try:
                baseIndex = baseToIndex[base]
            except KeyError:

                # Handle insertions
                if len(base) > 1:
                    baseIndex = 4
                    baseToIndex[base] = 4
                elif base == "-":
                    baseIndex = 4  # We can use the same index, since an insertion and deletion will never be listed by the same variant
                    baseToIndex["-"] = 4
                else:
                    # i.e. There is an odd base here. Probably an "N". In this case, ignore it
                    continue

            if qual > self.maxQual[baseIndex]:
                self.maxQual[baseIndex] = qual
            # What type of read is this, in terms of parental strand
            if posParent:
                if fSize >= strongMoleculeThreshold:
                    pType = "P"
                else:
                    pType = "p"
                self.posMolec[baseIndex] += 1
            else:
                if fSize >= strongMoleculeThreshold:
                    pType = "N"
                else:
                    pType = "n"
                self.negMolec[baseIndex] += 1

            # Calculate Strand bias info
            if map == "S":  # This base is a consensus from both the strand mapped to the (+) and (-) strand
                pMapStrand[baseIndex] += 1
                nMapStrand[baseIndex] += 1
            elif map== "F":
                pMapStrand[baseIndex] += 1
            else:
                nMapStrand[baseIndex] += 1

            # Check to see if there is another read that is the duplex mate of this read
            if readID in self.readParentalStrand:  # There is another read in duplex
                if readBaseIndex[readID] != baseToIndex[base]:  # i.e. the base at this position disagrees between the duplex. Discard this duplex
                    otherBaseIndex = readBaseIndex[readID]
                    self.disagreeingDuplexes[otherBaseIndex] += 1
                    self.disagreeingDuplexes[baseToIndex[base]] += 1
                    depth -= 1
                    self.strandCounts[otherBaseIndex] -= 1
                    del self.readParentalStrand[readID]
                    del readBaseIndex[readID]
                    if readID in readStats:
                        del readStats[readID]
                else:  # The bases agree between the two parental molecules. They are a valid duplex
                    self.readParentalStrand[readID] = self.readParentalStrand[readID] + pType
                    if base != self.ref:
                        readStats[readID].append([qual, fSize, distToEnd, mismatchNum, mapQual])  # TODO: Make this implementation less hacky
            else:
                self.readParentalStrand[readID] = pType
                readBaseIndex[readID] = baseToIndex[base]
                self.strandCounts[baseIndex] += 1
                # Total molecule count (for VAF)
                depth += 1
                if base != self.ref:
                    readStats[readID] = [[qual, fSize, distToEnd, mismatchNum, mapQual]]

        # Generate a count of each parental molecule at this position
        for readID, pStrand in self.readParentalStrand.items():
            baseIndex = readBaseIndex[readID]
            pStrand = set(pStrand)  # To deal with cases where ClipOverlap was NOT run, and there is double-counting
            if pStrand == {"N"}:
                self.baseCounts["SN"][baseIndex] += 1
            elif pStrand == {"P"}:
                self.baseCounts["SP"][baseIndex] += 1
            elif pStrand == {"n"}:
                self.baseCounts["Sn"][baseIndex] += 1
            elif pStrand == {"p"}:
                self.baseCounts["Sp"][baseIndex] += 1
            elif pStrand == {"P", "N"} or pStrand == {"N", "P"}:
                self.baseCounts["DPN"][baseIndex] += 1
            elif pStrand == {"P", "n"} or pStrand == {"n", "P"}:
                self.baseCounts["DPn"][baseIndex] += 1
            elif pStrand == {"p", "N"} or pStrand == {"N", "p"}:
                self.baseCounts["DpN"][baseIndex] += 1
            elif pStrand == {"p", "n"} or pStrand == {"n", "p"}:
                self.baseCounts["Dpn"][baseIndex] += 1
            else:
                raise TypeError("It seems the developer messed up. You should send him an angry e-mail! (Mentioning this message please)")

        # Calculate maximum alt quality
        self.maxAltQual = max(self.maxQual[baseToIndex[x]] for x in self.altAlleles)

        # Calculate strand bias
        self.strandBias = []
        if len(self.ref) == 1:
            refIndex = baseToIndex[self.ref]
            for index in range(0, 5):

                if pMapStrand[index] + nMapStrand[index] == 0:
                    self.strandBias.append(1.0)
                else:
                    self.strandBias.append(pvalue(pMapStrand[index], nMapStrand[index], pMapStrand[refIndex], nMapStrand[refIndex]).two_tail)
        else:  # This event is a deletion
            self.strandBias = [1,1,1,1]
            self.strandBias.append(
                pvalue(pMapStrand[4], nMapStrand[4], sum(pMapStrand[x] for x in range(0,4)), sum(nMapStrand[x] for x in range(0,4))).two_tail)
        self.vafs = []
        for index in range(0,5):
            try:
                self.vafs.append(self.strandCounts[index]/depth)
            except ZeroDivisionError:  # The only reads which overlapped this position were in duplex, and disagreed at this position
                self.vafs = [0,0,0,0,0]

        # Finalize the filtering stats
        completeStats = []
        altBases = []
        for readID, stats in readStats.items():
            baseIndex = readBaseIndex[readID]
            altCount = self.strandCounts[baseIndex]
            if len(stats) == 2:
                inDuplex = 0
            else:
                inDuplex = 1
            strandBias = self.strandBias[baseIndex]
            for stat in stats:
                stat.extend((altCount, depth, altCount/depth, strandBias, inDuplex))
                completeStats.append(stat)
                altBases.append(indexToBase[baseIndex])

        return completeStats, altBases

class PileupEngine:
    """
    Generates a custom pileup using family characteristics
    """
    def __init__(self, inBAM, refGenome, targetRegions, homopolymerWindow=5, noiseWindow=150, pileupWindow = 1000, printPrefix="PRODUSE-CALL\t"):
        self._inFile=pysam.AlignmentFile(inBAM)
        self._refGenome = Fasta(refGenome, read_ahead=20000)
        self._captureSpace = self._loadCaptureSpace(targetRegions)
        self._candidateIndels = {}
        self.pileup = {}
        self.candidateVar = {}
        self.filteredVar = {}
        self._homopolymerWindow = homopolymerWindow
        self._noiseWindow = noiseWindow
        # The number of positions to store before collapsing previous positions
        # Reduce this number to lower memory footprint, but be warned that, if this falls below the length of a read,
        # Some positions may not be tallied correctly
        self._pileupWindow = pileupWindow

        self._refStart = 0
        self._chrom = None
        self._refWindow = ""

        self._indelReads = []

        self._printPrefix = printPrefix

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

    def _insideCaptureSpace(self, read):  # TODO: Once INDELS are handled, reconsider this to include reads which may contain an INDEL at the capture space boundry
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

    def __writeVcfHeader(self, file):
        """
        Prints a VCF header to the specified file

        :param file: A file object, open to writing
        """

        header = ["##fileformat=VCFv4.3",
                  "##reference=" + self._refGenome.filename]

        # Add the contig information
        for contig, features in self._refGenome.faidx.index.items():
            line = "##contig=<ID=" + contig + ",length=" + str(len(features)) + ">"
            header.append(line)

        # Add the info field info
        header.append('##INFO=<ID=DPN,Number=R,Type=Integer,Description="Duplex Support with Strong Positive and Strong Negative Consensus">')
        header.append('##INFO=<ID=DPn,Number=R,Type=Integer,Description="Duplex Support with Strong Positive and Weak Negative Consensus">')
        header.append('##INFO=<ID=DpN,Number=R,Type=Integer,Description="Duplex Support with Weak Positive and Strong Negative Consensus">')
        header.append('##INFO=<ID=Dpn,Number=R,Type=Integer,Description="Duplex Support with Weak Positive and Weak Negative Consensus">')
        header.append('##INFO=<ID=SP,Number=R,Type=Integer,Description="Singleton Support with Strong Positive Consensus">')
        header.append('##INFO=<ID=Sp,Number=R,Type=Integer,Description="Singleton Support with Weak Positive Consensus">')
        header.append('##INFO=<ID=SN,Number=R,Type=Integer,Description="Singleton Support with Strong Negative Consensus">')
        header.append('##INFO=<ID=MC,Number=R,Type=Integer,Description="Total Molecule counts for Each Allele">')
        header.append('##INFO=<ID=STP,Number=R,Type=Float,Description="Probability of Strand Bias during sequencing for each Allele">')
        header.append('##INFO=<ID=PC,Number=R,Type=Integer,Description="Positive Strand Molecule Counts">')
        header.append('##INFO=<ID=NC,Number=R,Type=Integer,Description="Negative Strand Molecule Counts">')
        header.append('##INFO=<ID=VAF,Number=R,Type=Float,Description="Variant allele fraction of alternate allele(s) at this locus">')

        header.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO\n"]))
        file.write(os.linesep.join(header))

    def _processDeletion(self, delPos, moleculeTypes, chrom, position):
        """
        Merge positions that support a deletion into a consensus
        :param delPos: An iterable containing one or more position objects which overlap the deletion
        :return: outLine: A string representing a VCF entry for this deletion
        """
        def _generateInfoCol(startPos, endPos):
            infoCol = []
            for molecule in moleculeTypes:
                x = molecule + "=" + ",".join(
                    str(int((x + y) / 2)) for x, y in zip(startPos.baseCounts[molecule], endPos.baseCounts[molecule]))
                infoCol.append(x)

            infoCol.append("MC=" + ",".join(str(int((x + y) / 2)) for x, y in
                                            zip(startPos.strandCounts, endPos.strandCounts)))  # Parental Strand counts
            infoCol.append("STP=" + ",".join(
                str((x + y) / 2) for x, y in zip(startPos.strandBias, endPos.strandBias)))  # Strand bias
            infoCol.append("PC=" + ",".join(str((x + y) / 2) for x, y in zip(startPos.posMolec, endPos.posMolec)))
            infoCol.append("NC=" + ",".join(str((x + y) / 2) for x, y in zip(startPos.negMolec, endPos.negMolec)))
            infoCol.append("VAF=" + ",".join(str((x + y) / 2) for x, y in zip(startPos.vafs, endPos.vafs)))  # VAF
            return infoCol


        # My general strategy here is to use an average of the boundries of the deletion to obtain molecule counts

        # First, generate the reference base
        refSeq = "".join(x.ref for x in delPos)


        # Next, average the molecule counts on either side of the deletion, and use those for the output
        # Generate the info column
        infoCol = _generateInfoCol(delPos[0], delPos[-1])

        # Next, process the portion of this deletion which passed filters
        startIndex = 0
        endIndex = len(delPos) - 1
        delLength = len(delPos)
        while not delPos[startIndex].delPassedFilters and startIndex < delLength - 1:
            startIndex += 1
        while not delPos[endIndex].delPassedFilters and endIndex >= 0:
            endIndex -= 1
        endIndex += 1

        # If this deletion failed filters, we don't need to do anything else except state the reason why it failed filters
        if startIndex >= endIndex:
            # Use the reason listed at the first position
            fail = delPos[0].delFailReason
        else:
            fail = "PASS"
        rawOutLine = "\t".join([chrom, str(position + 1), ".", refSeq, "-", ".", fail, ";".join(infoCol)])

        # If this variant passed filters, we need to generate a new record for the passed deletion if it is not the same
        # as the old one
        if fail != "PASS":
            passedOutLine = None
        elif startIndex == 0 and endIndex == len(delPos):
            passedOutLine = rawOutLine
        else:
            refSeq = "".join(x.ref for x in delPos[startIndex:endIndex])
            infoCol = _generateInfoCol(delPos[startIndex], delPos[endIndex - 1])
            passedOutLine = "\t".join([chrom, str(position + startIndex + 1), ".", refSeq, "-", ".", "PASS", ";".join(infoCol)])

        return rawOutLine, passedOutLine

    def filterAndWriteVariants(self, outfile, filter, unfilteredOut=None, strandBiasThreshold=0.05):
        """

        Filters candidate variants based upon specified characteristics, and writes the output to the output file
        This is a placeholder for the real filter, which will be developed at a later time

        :param outfile: A string containing a filepath to the output VCF file
        :param filter: A sklearn.ensemble.RandomForestClassifier() which is to be used to filter variants
        :param unfilteredOut: A sting containing a filepath to an output VCF file which will contain ALL variants, even those that do not pass filters
        :param strandBiasThreshold: A float representing the strand bias p-value threshold
        :return:
        """

        sys.stderr.write(
            "\t".join([self._printPrefix, time.strftime('%X'), "Filtering Candidate Variants\n"]))
        moleculeTypes = ("DPN", "DPn", "DpN", "Dpn", "SP", "Sp", "SN", "Sn")
        with open(outfile, "w") as o, open(unfilteredOut, "w") if unfilteredOut is not None else open(os.devnull, "w") as u:
            self.__writeVcfHeader(o)
            self.__writeVcfHeader(u)

            for chrom in self.candidateVar:
                for position, stats in self.candidateVar[chrom].items():

                    try:
                        posStats, altAlleles = stats.collapsePosition(strandBiasThreshold)

                    # A keyerror will occur if the Reference base is not an "A,C,G,T" (i.e. an "N" or "M")
                    except KeyError:
                        continue # For now

                    if len(altAlleles) == 0:
                        continue

                    # Generate the info column
                    infoCol = []
                    for molecule in moleculeTypes:
                        x = molecule + "=" + ",".join(str(x) for x in stats.baseCounts[molecule])
                        infoCol.append(x)

                    infoCol.append("MC=" + ",".join(str(int(x)) for x in stats.strandCounts))  # Parental Strand counts
                    infoCol.append("STP=" + ",".join(str(x) for x in stats.strandBias))  # Strand bias
                    infoCol.append("PC=" + ",".join(str(x) for x in stats.posMolec))
                    infoCol.append("NC=" + ",".join(str(x) for x in stats.negMolec))
                    infoCol.append("VAF=" + ",".join(str(x) for x in stats.vafs))  # VAF

                    # See if this variant passes filters or not
                    # Aggregate the stats required for the filter
                    filterResults = filter.predict(posStats)  # 0=Passed filters, 1=Failed Filters
                    # Check to see which reads passed and failed filters, and assign variants using reads which passed filters
                    passedAltAlleles = set()
                    for altAllele, result in zip(altAlleles, filterResults):
                        if result == 0:
                            passedAltAlleles.add(altAllele)

                    if len(passedAltAlleles) > 0:
                        filterResult="PASS"
                    else:
                        filterResult="FAIL"
                    # Add 1 to the position, to compensate for the fact that pysam uses 0-based indexing, while VCF files use 1-based
                    outLine = "\t".join([chrom, str(position + 1), ".", stats.ref, ",".join(set(altAlleles)), str(stats.maxAltQual), filterResult, ";".join(infoCol)])
                    u.write(outLine + os.linesep)
                    if filter == "PASS":
                        outLine = "\t".join([chrom, str(position + 1), ".", stats.ref, ",".join(passedAltAlleles),
                                             str(stats.maxAltQual), filterResult, ";".join(infoCol)])
                        o.write(outLine + os.linesep)

            for chrom in self._candidateIndels:

                for position, stats in self._candidateIndels[chrom].items():

                    filtersFailed = stats.summarizeVariant(minMolecules, minAltMolecules, strandBiasThreshold)

                    # Generate the info column
                    infoCol = []
                    for molecule in moleculeTypes:
                        x = molecule + "=" + ",".join(str(x) for x in stats.baseCounts[molecule])
                        infoCol.append(x)

                    infoCol.append("MC=" + ",".join(str(int(x)) for x in stats.strandCounts))  # Parental Strand counts
                    infoCol.append("STP=" + ",".join(str(x) for x in stats.strandBias))  # Strand bias
                    infoCol.append("PC=" + ",".join(str(x) for x in stats.posMolec))
                    infoCol.append("NC=" + ",".join(str(x) for x in stats.negMolec))
                    infoCol.append("VAF=" + ",".join(str(x) for x in stats.vafs))  # VAF

                    failedFilters = False
                    if len(filtersFailed) == 0:
                        filter = "PASS"
                    else:
                        filter = ",".join(filtersFailed)
                        failedFilters = True

                    # Add 1 to the position, to compensate for the fact that pysam uses 0-based indexing, while VCF files use 1-based
                    outLine = "\t".join([chrom, str(position + 1), ".", stats.ref, ",".join(stats.altAlleles), str(stats.maxAltQual), filter, ";".join(infoCol)])
                    u.write(outLine + os.linesep)
                    if not failedFilters:
                        outLine = "\t".join([chrom, str(position + 1), ".", stats.ref, ",".join(stats.passedAltAlleles),
                                             str(stats.maxPassedAltQual), filter, ";".join(infoCol)])
                        o.write(outLine + os.linesep)

        sys.stderr.write(
            "\t".join([self._printPrefix, time.strftime('%X'), "Variant Calling Complete\n"]))

    def _indelFromRealignment(self, origRead, realignRead):
        """
        Obtains INDEL positions a read realigned via SSW

        :param origRead: A pysam.AlignedSegment object
        :param realignRead: A skbio.alignment.AlignmentStrcture generated using the origRead.query_sequence
        :return:
        """
        insPos = {}
        delPos = {}
        length = ""
        currentPos = realignRead.query_begin
        seqPos = realignRead.target_begin
        for char in realignRead.cigar:
            if char.isdigit():
                length = length + char
            else:
                opLength = int(length)
                if char == "M":
                    currentPos += opLength
                    seqPos += opLength
                elif char == "D":  # This is actually an insertion in terms of the read (the cigar is relative to the reference)
                    insPos[currentPos] = realignRead.target_sequence[seqPos: seqPos + opLength + 1]
                    seqPos += opLength
                elif char == "I":  # Reversed due to cigar being relative to the reference
                    delPos[currentPos] = opLength
                    currentPos += opLength
                else:
                    raise TypeError("Unknown operator returned by skbio.alignment.AlignmentStructure: \"%s\". "
                                    "It is likely that skbio's API has been updated since this code was written "
                                    "(Generated for skbio=0.5.1)" % char)
                length = ""

        """
        outCigar = []
        insPos = {}
        delPos = {}
        
        
        # Where was this read expected to start?
        expectedStart = origRead.reference_start - self._refStart # Start of this read, relative to the reference window
        if origRead.cigartuples[0][0] == 4:
            expectedStart -= origRead.cigartuples[0][1]

        assert realignRead.query_start >= expectedStart
        outCigar.append((4, realignRead.target_begin))
        # Convert the alignmentStructure cigar string to a cigar tuple

        length = ""
        foundIndel = False
        currentPos = realignRead.query_start
        for char in realignRead.cigar:
            if char.isdigit():
                length = length + char
            else:
                opLength = int(length)
                if char == "M":
                    currentPos += opLength
                    outCigar.append((0, opLength))
                elif char == "D":
                    delPos[currentPos] = currentPos
                    outCigar.append((1, opLength))  # Invert operator, since this is the cigar of the reference sequence
                    currentPos += opLength
                    foundIndel = True
                elif char == "I":
                    outCigar.append((2, opLength))  # Invert operator, since this is the cigar of the reference sequence
                    foundIndel = True
                else:
                    raise TypeError("Unknown operator returned by skbio.alignment.AlignmentStructure: \"%s\". "
                                    "It is likely that skbio's API has been updated since this code was written "
                                    "(Generated for skbio=0.5.1)" % char)
                length = ""

        # Add any trailing soft-clipping
        expectedEnd = len(origRead.query_sequence)
        assert expectedEnd >= realignRead.target_end_optimal
        outCigar.append((4, expectedEnd - realignRead.target_end_optimal))

        return outCigar, foundIndel
        """
        return delPos, insPos

    def cigarToTuple(self, cigarTuples):
        """
        Converts a pysam-style cigar tuples into a regular tuple, where 1 operator = 1 base

        Ignores soft-clipping, and flags a
        :param cigarTuples: A list of tuples
        :return: A tuple listing expanded cigar operators, as well as a boolean indicating if there is an INDEL in this read
        """

        cigarList = []
        hasIndel = False
        for cigarElement in cigarTuples:
            cigOp = cigarElement[0]
            # Ignore soft-clipped bases
            if cigOp == 4:
                continue
            # If there is an insertion or deletion in this read, indicate that
            if cigOp == 1 or cigOp == 2:
                hasIndel = True
            cigarList.extend([cigOp] * cigarElement[1])

        return tuple(cigarList), hasIndel

    def _combineIndel(self, events):
        """
        Attempt to combine multiple INDEL events

        :param events: A distionary containing a dictionary of reference sequences supporting the event
        :return:
        """

        # Generate a baseline k-mer index
        # refKmers = IndelUtils.SeqIndexSet(self._refWindow)

        refAlign = alignment.StripedSmithWaterman(self._refWindow)

        distToRef = sortedcontainers.SortedDict()

        # For each event, obtain a baseline score relative to the reference
        # We will use this to perform a comparison later with all other events
        for pos, eventsAtPos in events.items():
            for eventKey, event in eventsAtPos.items():
                mapping = refAlign(event)

                # Is there already an event with the same mapping score?
                # If so, see if they are smiliar enough that they can be combined
                refScore = mapping.optimal_alignment_score
                if refScore in distToRef:
                    readAlign = alignment.StripedSmithWaterman(event)
                    canMerge = False
                    for altEvent in distToRef[refScore]:
                        altSeq = events[altEvent[0]][altEvent[1]]
                        altMapping = readAlign(altSeq)

                        # If this current event is closer to this event than the reference, we can just use the
                        # existing alternate seq as a template, and thus merge these events
                        if altMapping.optimal_alignment_score < refScore:
                            canMerge = True
                            break

                    if canMerge:
                        continue
                    else:
                        distToRef[refScore].append((pos, eventKey))
                else:
                    distToRef[refScore] = [(pos, eventKey)]

        collapsedEvents = {}
        # Now that we have the distance from the reference for each event, work outwards
        # compare the event closest to the ref with all other events, and see if they be merged
        for distance, altEvents in distToRef.items():
            for location in altEvents:
                # Create a mapping for this event
                altSeq = events[location[0]][location[1]]
                mapping = alignment.StripedSmithWaterman(altSeq)

                distToDelete = []
                # Compare this event to all other events, and see if they can be merged
                for otherDistance, otherAltEvents in distToRef.items():
                    i = 0
                    posToDelete = sortedcontainers.SortedList()
                    for otherLocation in otherAltEvents:
                        # If this event is currently being processed, ignore it
                        if distance == otherDistance and location == otherLocation:
                            continue

                        otherAltSeq = events[otherLocation[0]][otherLocation[1]]
                        altMap = mapping(otherAltSeq)
                        if altMap.optimal_alignment_score < distance:
                            posToDelete.add(i)
                        i += 1

                    # Delete any events that are no longer needed
                    for i in reversed(posToDelete):
                        otherAltEvents.pop(i)
                        if len(otherAltEvents) == 0:
                            distToDelete.append(otherDistance)

                # Delete and reference distances which are to be collapsed
                for i in distToDelete:
                    del distToRef[i]

                # Now that all possible events have been collapsed into this one, store this event
                if location[0] not in collapsedEvents:
                    collapsedEvents[location[0]] = {}
                collapsedEvents[location[0]][location[1]] = events[location[0]][location[1]]

        return collapsedEvents

    def _getAltSeq(self):
        """
        Identify the sequence indicated by each INDEL event.
        If possible, combine multiple sequences into a consensus

        :return:
        """

        indels = {}
        refRealign = None

        # In the easiest case, the aligner has already marked each event for us, and we can just use those
        for read in self._indelReads:
            eventStart = read.reference_start - self._refStart
            readStart = 0
            delEvents = {}
            insEvents = {}
            numEvents = 0
            for operation, length in read.cigartuples:

                if operation == 0:  # Add this offset to the event start
                    eventStart += length
                    readStart += length
                elif operation == 4:  # Soft clipping is not included in the offset
                    readStart += length
                elif operation == 2:  # This event is a deletion
                    delEvents[eventStart] = length
                    eventStart += length
                    numEvents += 1
                elif operation == 1:  # This event is an insertion
                    insEvents[eventStart] = read.query_sequence[readStart: readStart + length + 1]
                    numEvents += 1
                    readStart += length

            # If there were several events, or the aligner didn't find any, perform local realignment
            if numEvents != 1:  # We could possibly find a simplier solution if more than one indel was called
                if refRealign is None:  # Generate a SW object from the reference, if one does not exist already
                    refRealign = alignment.StripedSmithWaterman(self._refWindow)
                readRealign = refRealign(read.query_sequence)
                delEvents, insEvents = self._indelFromRealignment(read, readRealign)

            # Store these events
            for position, event in delEvents.items():  # Store deletions
                delKey = "-" * event
                if position in indels:
                    # See if this deletion is the same as any other deletions at this position
                    if event not in indels[position]:  # Add a reference which contains this event
                        indels[position][delKey] = self._refWindow[:position] + self._refWindow[position + event:]
                else:
                    indels[position] = {}
                    indels[position][delKey] = self._refWindow[:position] + self._refWindow[position + event:]

            for position, event in insEvents.items():  # Store insertions
                if position in indels:
                    # See if this insertion is the same as any other insertion
                    if event not in indels[position]:
                        indels[position][event] = self._refWindow[:position] + event + self._refWindow[position:]
                else:
                    indels[position] = {}
                    indels[position][event] = self._refWindow[:position] + event + self._refWindow[position:]

        return indels

    def _cigarFromRealignment(self, realignStructure, read, eventPos, event):
        """
        Obtains a new pysam style cigar tuple from a SSW alignment structure

        :param realignStructure: A SSW alignment structure
        :param read: A pysam.AlignedSegment()
        """

        cigarToOp = {"M": 0, "D": 1, "I": 2}  # Reversed, as the SSW cigar is actually relative to the reference
        outCigar = []
        # If there was no soft-clipping, where would this read start?
        startPos = read.reference_start
        # Compensate for leading soft clipping in the original alignment
        if read.cigartuples[0][0] == 4:
            startPos -= read.cigartuples[0][1]

        # How much soft clipping is in the original alignment?
        startSoftClip = realignStructure.target_begin
        if startSoftClip != 0:
            outCigar.append((4, startSoftClip))
            startPos += startSoftClip
            # If the event is inside the soft-clipped region, then something catastophic happened...
            # In this case, preserve the validity of the original alignment

        cigarPos = startPos - self._refStart  # Use this to add the event back into the coresponding position
        if cigarPos >= eventPos:
            raise TypeError("ERROR during INDEL realignment")
        # Process the cigar
        length = ""
        for char in realignStructure.cigar:
            if char.isdigit():
                length = length + char
            else:
                opLength = int(length)
                if cigarPos < eventPos and cigarPos + opLength > eventPos:  # Aka the event itself falls within this region
                    x = eventPos - cigarPos
                    outCigar.append((cigarToOp[char], x))
                    # Add this event to the cigar
                    if set(event) == {"-"}:  # This event is a deletion
                        outCigar.append((2, len(event)))
                        x = cigarPos + opLength - eventPos
                    else:
                        outCigar.append((1, len(event)))
                        x = cigarPos + opLength - eventPos - len(event)

                        outCigar.append((cigarToOp[char], x))
                elif opLength > 0:
                    outCigar.append((cigarToOp[char], opLength))
                length = ""
                if char != "D":  # Insertions do not consume reference positions, so do not add to the counter
                    cigarPos += opLength

        # Compensate for any trailing soft-clipping
        endSoftClip = len(realignStructure.target_sequence) - realignStructure.target_end_optimal - 1
        if endSoftClip != 0:
            outCigar.append((4, endSoftClip))

        return tuple(outCigar), startPos

    def _mapReadsToIndels(self, events, realignWindow=500):
        """
        Map all reads which may contain an INDEL to the collapsed INDELs

        TODO: Handle reads which support more than one event
        :param events: A dictionary storing positions: events
        :param realignWindow: How close to the start of a read an event must be before it is considered for ealignment
        :return:
        """
        sswStructures = {}
        # Create a SSW alignment structure for each alternate allele to improve performance
        for position, eventAtPos in events.items():
            sswStructures[position] = {}
            for event, seq in eventAtPos.items():
                sswStructures[position][event] = alignment.StripedSmithWaterman(seq)

        refAlignment = alignment.StripedSmithWaterman(self._refWindow)

        indelVar = {}

        for read in self._indelReads:
            # Obtain a baseline measure of how well this read compares to the reference
            mapping = refAlignment(read.query_sequence)
            bestScore = mapping.optimal_alignment_score
            bestReadCigar = read.cigartuples
            bestStart = read.reference_start
            bestEvent = None
            bestPosition = None


            # Map this read to each INDEL, and see which one it aligns to best
            for position, eventAtPos in events.items():
                if position - realignWindow > read.reference_start or position + realignWindow < read.reference_start:
                    continue
                for event in eventAtPos.keys():
                    altAlignment = sswStructures[position][event](read.query_sequence)

                    if altAlignment.optimal_alignment_score > bestScore:
                        # This event is now the best event
                        bestScore = altAlignment.optimal_alignment_score
                        # bestReadCigar, bestStart = self._cigarFromRealignment(altAlignment, read, position, event)
                        bestEvent = event
                        bestPosition = position
            # Now that we have the best event for this read, add this read to the variant, and assign the new cigar to the read
            if bestEvent is not None:  # This read supports a non-reference event
                isDel = set(bestEvent) == {"-"}
                refPos = bestPosition + self._refStart
                # Modify the read
                # read.reference_start = bestStart
                # read.cigartuples = bestReadCigar

                if self._chrom not in self._candidateIndels:
                    self._candidateIndels[self._chrom] = sortedcontainers.SortedDict()

                if refPos not in self._candidateIndels[self._chrom]:  # This event has not been flagged as a variant yet
                    # Create a position object which will store this INDEL
                    if isDel:  # Obtain the reference window for this deletion
                        ref = self._refWindow[bestPosition:bestPosition + len(bestEvent) + 1]
                    else:
                        ref = self._refWindow[bestPosition]
                    self._candidateIndels[self._chrom][refPos] = Position(ref)

                if isDel:
                    alt = "-"
                    qual = 0
                else:
                    alt = self._candidateIndels[self._chrom][refPos].ref + bestEvent
                    readIndelPos = refPos - read.reference_start
                    qual = int(sum(x for x in read.query_alignment_qualities[readIndelPos: readIndelPos + len(alt)]) / len(alt))

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
                # If clipoverlap was run on this BAM file, the originating strand of each base will be stored in a tag
                try:  # Just in case this read does not contain that tag
                    rMapStrand = read.get_tag("co")
                except KeyError:  # i.e. The required tag is not present. Either these reads do not overlap at all, or clipoverlap was never run on the input BAM file
                    if read.is_reverse:
                        rMapStrand = ["R"] * len(read.query_sequence)
                    else:
                        rMapStrand = ["F"] * len(read.query_sequence)

                rMapStrand = rMapStrand[refPos - read.reference_start]
                rMappingQual = read.mapping_quality

                # Read stats for this event
                dist = min(refPos - read.reference_start, read.reference_end - refPos)
                self._candidateIndels[self._chrom][refPos].add(alt, qual, rFSize, rPosParent, rMapStrand, dist, rMappingQual, counter)

    def _findIndels(self):

        # We need to find what Insertion and Deletions are represented by these reads
        """
        events = self._getAltSeq()  # Obtain INDELs
        if len(events) != 0:
            events = self._combineIndel(events)  # Combine INDELs

            # Perform local realignment of reads, and identify what INDELs are supported
            self._mapReadsToIndels(events)
        # Finally, map the rest of these reads back to the reference
        for read in self._indelReads:
            rCigar, hasIndel = self.cigarToTuple(read.cigartuples)
            self.processRead(read, rCigar, mismatchMax=1.1)
        """
        self._indelReads = []

    def processRead(self, read, rCigar, mismatchMax=0.02, discardBad=False):
        """
        Add all positions covered by this read to the pileup

        :param read: A pysam.AlignedSegment() object
        :param rCigar: A list containing cigar operators
        :param mismatchMax: The lead length multiplied this is the maximum number of mismatched permitted before the read is flagged as supporting an indel
        """

        maxMismatch = len(read.query_alignment_sequence) * mismatchMax + 1
        if maxMismatch < 5:
            maxMismatch = 5

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
        # If clipoverlap was run on this BAM file, the originating strand of each base will be stored in a tag
        try:  # Just in case this read does not contain that tag
            rMapStrand = read.get_tag("co")
        except KeyError:  # i.e. The required tag is not present. Either these reads do not overlap at all, or clipoverlap was never run on the input BAM file
            if read.is_reverse:
                rMapStrand = ["R"] * len(read.query_alignment_sequence)
            else:
                rMapStrand = ["F"] * len(read.query_alignment_sequence)

        rMappingQual = read.mapping_quality

        # All positions covered by this read
        positions = []

        # The number of variants supported by this read
        numMismatch = 0

        # Iterate through all bases and cigar indexes
        cigarIndex = 0  # Keep this seperate to account for deletions
        readIterator = zip(read.query_alignment_sequence, read.query_alignment_qualities,
                           rMapStrand)  # TODO: Fix this
        position = read.reference_start

        try:
            while True:
                # for base, cigar, qual, mapStrand in zip(read.query_alignment_sequence, rCigar, read.query_alignment_qualities, rMapStrand):
                cigar = rCigar[cigarIndex]
                cigarIndex += 1
                if cigar == 1:
                    next(readIterator)
                    continue
                elif cigar == 2:
                    position += 1
                    continue

                # If this position has not been covered before, we need to generate a new pileup at this position
                if position not in self.pileup[self._chrom]:
                    refBase = self._refWindow[position - self._refStart]
                    self.pileup[self._chrom][position] = Position(refBase)

                pileupPos = self.pileup[self._chrom][position]
                position += 1

                base, qual, mapStrand = next(readIterator)

                distFromEnd = min(position - read.reference_start, read.reference_end - position)

                # Ignore bases that are "N"
                if base == "N":
                    continue
                isAlt = pileupPos.add(base, qual, rFSize, rPosParent, mapStrand, distFromEnd, rMappingQual, counter)

                # Check to see if this position now has evidence of a variant
                if isAlt:
                    numMismatch += 1
                positions.append(pileupPos)

        except IndexError:
            pass
        except StopIteration as e:
            print(read)
            raise TypeError("Malformed Read Detected: %s" % read.query_name) from e

        # If this read has too many mismatches, it might have an indel which was not called.
        # Remove this read from the pileup
        if numMismatch > maxMismatch:
            for pos in positions:
                pos.removeLast()
            self._indelReads.append(read)
        else:
            # Store the number of total mismatches for this read at each variant
            for pos in positions:
                pos.addMismatchNum(numMismatch)

    def generatePileup(self, windowBuffer = 200):

        def processWindow(pos, coordinate):

            # If this position is an INDEL, add the statistics from this position to that INDEL
            if self._chrom in self._candidateIndels and coordinate in self._candidateIndels[self._chrom]:
                self._candidateIndels[self._chrom][coordinate].addPos(pos)
            # Ignore non-variant positions
            if pos.alt:

                # Is this position within the capture space?
                if self._captureSpace is not None:
                    if self._chrom not in self._captureSpace or bisect.bisect_right(self._captureSpace[self._chrom], coordinate) % 2 != 1:
                        return  # This position does not fall within the capture space. Ignore it

                # If there is a variant, we need to obtain some general characteristics to be used for filtering
                # First, store the sequence of the surrounding bases. This will be used to identify errors
                # due to homopolymers
                homoWindow = []

                windowPos = coordinate - self._refStart # Where, relative to the reference window, are we located?
                for j in range(windowPos - self._homopolymerWindow, windowPos + self._homopolymerWindow):
                    if windowPos == j:  # i.e. We are at the position we are currently examining, flag it so we have a reference point
                        homoWindow.append("M")  # It's me!
                    else:
                        try:
                            homoWindow.append(self._refWindow[j])
                        except IndexError:  # In the rare case where we overflow outside the window
                            homoWindow.append(self._refGenome[self._chrom][j + self._refStart])

                # Next, examine how "noisy" the neighbouring region is (i.e. how many candidiate variants there are)
                # Flag all nearby candidate variants. To avoid flagging the same variant twice, add this
                # candidate variant position to all variants behind it, and vise-versa
                if self._chrom not in self.candidateVar:
                    self.candidateVar[self._chrom] = sortedcontainers.SortedDict()

                nearbyVar = []
                for j in self.candidateVar[self._chrom].irange(coordinate - self._noiseWindow):
                    self.candidateVar[self._chrom][j].nearbyVar.append(pos)
                    nearbyVar.append(self.candidateVar[self._chrom][j])

                # Finally, save all these stats we have just generated inside the variant, to be examined during
                # filtering
                pos.processVariant(homoWindow, nearbyVar, self._noiseWindow)
                self.candidateVar[self._chrom][coordinate] = pos


        # Print status messages
        sys.stderr.write("\t".join([self._printPrefix, time.strftime('%X'), "Starting...\n"]))
        sys.stderr.write("\t".join([self._printPrefix, time.strftime('%X'), "Finding Candidate Variants\n"]))
        readsProcessed = 0
        lastPos = -1

        for read in self._inFile:

            # Prior to adding each base in this read onto the pileup, we need to analyze the current read,
            # and obtain general characteristics (family size, mapping quality etc)

            # Print out status messages
            readsProcessed += 1
            if readsProcessed % 100000 == 0:
                sys.stderr.write("\t".join([self._printPrefix, time.strftime('%X'), "Reads Processed:%s\n" % readsProcessed]))

            # Ignore unmapped reads
            if read.is_unmapped:
                continue

            # Ignore reads without a cigar string, as they are malformed or empty
            if not read.cigartuples:
                continue

            # Have we advanced positions? If so, we may need to change the loaded reference window, and process
            # positions which no more reads will be mapped to
            if read.reference_name != self._chrom or read.reference_start != lastPos:

                # Sanity check that the input BAM is sorted
                if read.reference_name == self._chrom and read.reference_start < lastPos:
                    sys.stderr.write("ERROR: The input BAM file does not appear to be sorted\n")
                    exit(1)

                # Process positions at which no more reads are expected to map to (assuming the input is actually sorted)
                # We do this to drastically reduce the memory footprint of the pileup, since we don't care about non-variant
                # positions very much
                try:
                    if read.reference_name != self._chrom and self._chrom is not None:

                        # We need to process all remaining reads which support an INDEL
                        self._findIndels()
                        # If we have switched chromosomes, process all variants stored on the previous chromosome
                        posToProcess = self.pileup[self._chrom].keys()
                        for coordinate in posToProcess:
                            processWindow(self.pileup[self._chrom][coordinate], coordinate)
                        del self.pileup[self._chrom]
                        # Since this is the first read we are processing from this chromosome, we need to create the dictionary which
                        # will store all bases on this chromosome
                        self.pileup[read.reference_name] = sortedcontainers.SortedDict()

                    # Otherwise, check to see if previous positions fall outside the buffer window. If so, process them
                    elif self._chrom is not None:

                        # We need to process all remaining reads which support an INDEL
                        self._findIndels()
                        posToProcess = tuple(self.pileup[self._chrom].irange(minimum=None, maximum=read.reference_start - self._pileupWindow))
                        for coordinate in posToProcess:
                            processWindow(self.pileup[self._chrom][coordinate], coordinate)
                        for coordinate in posToProcess:
                            del self.pileup[self._chrom][coordinate]
                    elif self._chrom is None:
                        self.pileup[read.reference_name] = sortedcontainers.SortedDict()

                except KeyError:  # i.e. All reads which mapped to this contig or previous positions failed QC. There are no positions to process
                    pass

                # Load the current reference window and it's position
                # If we have switched contigs, or moved outside the window that is currently buffered,
                # we need to re-buffer the reference
                readWindowStart = read.reference_start - windowBuffer
                if readWindowStart < 0:
                    readWindowStart = 0
                readWindowEnd = read.reference_end + windowBuffer  #TODO: Don't estimate soft-clipping, but instead calculate it

                if read.reference_name != self._chrom or readWindowStart < self._refStart or readWindowEnd > refEnd:
                    # To reduce the number of times this needs to be performed, obtain a pretty large buffer
                    self._refStart = readWindowStart - self._pileupWindow - self._homopolymerWindow
                    if self._refStart < 0:
                        self._refStart = 0
                    refEnd = readWindowStart + 5000
                    self._chrom = read.reference_name
                    self._refWindow = self._refGenome[self._chrom][self._refStart:refEnd].seq

                lastPos = read.reference_start

            # Does this read overlap the capture space? If not, ignore it
            # Note that this will keep reads which partially overlap the capture space
            # We will perform a more precise restricting of the regions which overlap the capture space later
            if not self._insideCaptureSpace(read):
                continue


            # Convert the cigar to a more user-friendly format
            rCigar, rHasIndel = self.cigarToTuple(read.cigartuples)

            # If this read supports an INDEL, we need to set it aside and process it separately.
            if rHasIndel:
                self._indelReads.append(read)
                continue

            self.processRead(read, rCigar)

        # Now that we have finished processing all reads in the input file, process all remaining positions
        try:
            posToProcess = self.pileup[self._chrom].keys()
            for coordinate in posToProcess:
                processWindow(self.pileup[self._chrom][coordinate], coordinate)
        except KeyError:
            pass

        sys.stderr.write("\t".join([self._printPrefix, time.strftime('%X'), "Reads Processed:" + str(readsProcessed) + "\n"]))
        sys.stderr.write(
            "\t".join([self._printPrefix, time.strftime('%X'), "Candidate Variants Identified\n"]))


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
        raise parser.error("Unable to locate %s. Please ensure the file exists, and try again" % file)


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
        if parameter == "True":
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
    parser.add_argument("-o", "--output", metavar="VCF", required=True, help="Output VCF file")
    parser.add_argument("-u", "--unfiltered", metavar="VCF",
                        help="An additional output VCF file, which lists all candidate variants, including those that failed filters")
    parser.add_argument("-r", "--reference", metavar="FASTA", required=True, type=lambda x: isValidFile(x, parser),
                        help="Reference Genome, in FASTA format")
    parser.add_argument("-t", "--targets", metavar="BED", type=lambda x: isValidFile(x, parser),
                        help="A BED file containing regions in which to restrict variant calling")
    parser.add_argument("-f", "--filter", metavar="PICKLE", required=True, type=lambda x: isValidFile(x, parser),
                        help="A pickle file containing the trained filter. Can be generated using \'produse train\'")
    parser.add_argument("--strand_bias_threshold", metavar="FLOAT", type=float, default=0.05,
                        help="Any variants with a strand bias above this threshold will be filtered out")
    validatedArgs = parser.parse_args(listArgs)
    return vars(validatedArgs)


parser = argparse.ArgumentParser(description="Identifies and calls variants")
parser.add_argument("-c", "--config", metavar="INI", type=lambda x: isValidFile(x, parser), help="An optional configuration file which can provide one or more arguments")
parser.add_argument("-i", "--input", metavar="BAM", type=lambda x: isValidFile(x, parser), help="Input post-collapse or post-clipoverlap BAM file")
parser.add_argument("-o", "--output", metavar="VCF", help="Output VCF file, listing all variants which passed filters")
parser.add_argument("-u", "--unfiltered", metavar="VCF", help="An additional output VCF file, which lists all candidate variants, including those that failed filters")
parser.add_argument("-r", "--reference", metavar="FASTA", type=lambda x: isValidFile(x, parser), help="Reference Genome, in FASTA format")
parser.add_argument("-t", "--targets", metavar="BED", type=lambda x: isValidFile(x, parser), help="A BED file containing regions in which to restrict variant calling")
parser.add_argument("-f", "--filter", metavar="PICKLE", type=lambda x: isValidFile(x, parser),
                    help="A pickle file containing a trained filter. Can be generated using \'produse train\'")
parser.add_argument("--strand_bias_threshold", metavar="FLOAT", type=float, help="Any variants with a strand bias above this threshold will be filtered out")


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
                if argument in args and not args[argument]:  # Aka this argument is used by call, and a parameter was not provided at the command line
                    args[argument] = parameter
        except KeyError:  # i.e. there is no section named "call" in the config file
            sys.stderr.write(
                "ERROR: Unable to locate a section named \'call\' in the config file \'%s\'\n" % (args["config"]))
            exit(1)

    args = validateArgs(args)
    pileup = PileupEngine(args["input"], args["reference"], args["targets"], printPrefix=printPrefix)

    # Find candidate variants
    pileup.generatePileup()

    # Filter variants
    with open(args["filter"], "r+b") as o:
        filterModel = pickle.load(o)
    pileup.filterAndWriteVariants(args["output"], filterModel, args["unfiltered"], args["strand_bias_threshold"])


if __name__ == "__main__":
    main()
