The ProDuSe Pipeline
====================

The ProDuSe pipeline consists of multiple steps, each of which is described breifly here

Config
^^^^^^

Sets up output directories, file symbolic links, and configuration files for each step in the ProDuSe pipeline

.. seealso:: `configure_produse.py`_

.. _configure_produse.py: Config.html


Trim
^^^^

Trims the adapter sequence off of each read, and prepends the adapter sequence to the read name. In addition, any reads which
do not match the barcode sequence (and are outside the mismatch range) are discarded

    - Input: Paired fastq files
    - Output: Trimmed paired fastq files

.. note:: Trim can be used to demultiplex samples, assuming the barcodes corresponding to each sample are sufficiently distinct

.. seealso:: `trim.py`_

.. _trim.py: Trim.html

bwa
^^^

Maps provided reads to a reference genome using the Burrows-Wheeler Aligner.

    - Input: Trimmed paired fastq files
    - Output: BAM file consisting of trimmed reads

.. seealso:: `bwa.py`_

.. _bwa.py: Bwa.html

Collapse
^^^^^^^^

Collapses duplicate reads into a consensus sequence.

    - Input: Trimmed BAM file
    - Output: Paired collapsed (consensus) fastq files

.. seealso:: `collapse.py`_

.. _collapse.py: Collapse.html

bwa
^^^

Maps provided reads to a reference genome using the Burrows-Wheeler Aligner.

    - Input: Collapsed paired fastq files
    - Output: BAM file consisting of Collapsed reads

.. seealso:: `bwa.py`_

.. _bwa.py: Bwa.html

Stitcher
^^^^^^^^

Merges forward and reverse reads into a consensus sequence if they overlap. Used to correct errors in overlapping bases.

    - Input: Collapsed BAM file
    - Output: Stitched (forward and reverse reads are merged) BAM file

.. seealso:: `Stitcher Documentation`_

.. _Stitcher Documentation: https://github.com/Illumina/Pisces/wiki/Stitcher-5.1.6-Design-Document

SplitMerge
^^^^^^^^^^

Splits reads merged by Stitcher back into forward and reverse reads. The shared bases are assigned to one of the two reads.

    - Input: Stitched BAM file
    - Output: BAM file with merged reads split. Unsorted

.. seealso:: `SplitMerge.py`_

.. _SplitMerge.py: SplitMerge.html

SNV
^^^

Identifies positions whereby one or more bases support an alternate allele.

    - Input: A BAM file containing de-stitched reads
    - Output: A VCF file listing all positions with non-reference bases, and the count of these bases

.. seealso:: `snv.py`

.. _snv.py: snv.html

Filter
^^^^^^

Filters candidate variant calls based upon overall capture space and locus characteristics.

    - Input: A VCF file listing candidate variants
    - Output: A VCF file listing filtered variants

.. seealso:: `filter_produse.py`_

.. _filter_produse.py: Filter.html


