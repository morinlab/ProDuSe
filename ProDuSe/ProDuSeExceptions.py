#! /usr/bin/env python

"""
Contains all custom exceptions thrown by ProDuSe

"""

# Thrown if a FASTA or BAM index is missing
class IndexNotFound(Exception):
    pass


# Thrown if a read (pysam.AlignedSegment) doesn't appear valid
# Usually this is thrown if the cigar, sequence, and qualities scores are not compatible
class MalformedReadException(Exception):
    pass


# Thrown if the input BAM/SAM/BED file is unsorted
class UnsortedInputException(Exception):
    pass


# Thrown if the input file is not in the required format
class InvalidInputException(Exception):
    pass
