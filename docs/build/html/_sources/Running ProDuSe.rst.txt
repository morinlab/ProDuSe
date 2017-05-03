Running ProDuSe
===============

The ProDuSe pipeline consists of multiple steps, each of which is described breifly here

ProDuSe Config
^^^^^^^^^^^^^^

Sets up output directories, file symlinks, and configuration files for each step in the ProDuSe pipeline

.. seealso:: `ProDuSe Config`_

.. _ProDuSe Config: Config.html


ProDuSe Trim
^^^^^^^^^^^^

Trims the adapter sequence off of each read, and prepends the adapter sequence to the read name. In addition, any reads which
do not match the barcode sequence (and are outside the mismatch range) are discarded

    - Input: Forward and reverse fastq files
    - Output: Trimmed forward and reverse fastq files

.. note:: Trim can be used to demultiplex samples, assuming the barcodes coresponding to each sample are sufficiently distinct

.. seealso:: `ProDuSe Trim`_

.. _ProDuSe Trim: Trim.html
