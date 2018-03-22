The ProDuSe Pipeline
====================

The ProDuSe pipeline consists of multiple stages, each of which is described breifly here.

Trim
^^^^

Trims the barcode sequence off each read, and stores it in a FASTQ comment. Any reads where the barcode deviates significantly
from the expected degenerate range are discarded

    - Input: Paired fastq files
    - Output: Trimmed paired fastq files

.. note:: Trim can be used to demultiplex samples, assuming the barcodes used in each sample are sufficiently distinct

.. seealso:: `Trim`_

.. _Trim: Trim.html

bwa
^^^

Maps provided reads to a reference genome using the Burrows-Wheeler Aligner (mem algorithm). The resulting SAM file is converted into
a BAM file and sorted, with the FASTQ comment stored as a read tag

    - Command: bwa mem <reference> <trimmed_fastq.R1.fastq> <trimmed_fastq.R2.fastq> | samtools view -b | samtools sort > out.trimmed.bam

Collapse
^^^^^^^^

Collapses duplicate reads into a consensus sequence. In addition, reads which are in "duplex" (i.e. originate from the same parental
molecule) are flagged here

    - Input: Trimmed BAM file
    - Output: Collapsed BAM File

.. seealso:: `Collapse`_

.. _Collapse: Collapse.html

ClipOverlap
^^^^^^^^^^^

Idenfies bases that overlap between each read pair, and generates a consensus from the overlap. This
consensus is then assigned to only one read in the read pair, thus removing overlapping bases.

    - Input: Collapsed BAM file
    - Output: Clipped BAM file

.. seealso:: `ClipOverlap`_

.. _ClipOverlap: ClipOverlap.html

Call
^^^^^^

Identifies all possible variants in the specified file, then filters those variants
based upon capture-space and locus-specific characteristics.

    - Input: Clipped BAM file
    - Output: Two VCF files, listing raw (all) and filtered variants

.. seealso:: `Call`_

.. _Call: Call.html


