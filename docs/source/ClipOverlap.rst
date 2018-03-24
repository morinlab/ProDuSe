ClipOverlap
===========

Description
^^^^^^^^^^^

Identifies all positions which overlap inside of a given read pair. If any bases overlap, a consensus is generated
and that consensus is assigned to only one read in the read pair.

Run Using
^^^^^^^^^

::

	produse clip

or

::

	/path/to/ProDuSe/Clone/ProDuSe/ClipOverlap.py

Parameters
^^^^^^^^^^

	:-c --config:
		A configuration file which can provide any of the following arguments. See the `config page`_ for more info.
	:-i --input:
		An input SAM/BAM file (use "-" to read from standard in). Does not need to be sorted.
	:-o --output:
		An output SAM/BAM file in which to write clipped reads (use "-" to write to stdout). Will be UNSORTED. The file type is determined from the file extension, or the input file type if stdout is specified.
	:--tag_origin:
		Add a read tag indicating which read a consensus base originated. S=Both reads agree.

	.. _config page: Config_Files.html

.. warning:: The output BAM file will be unsorted, even if the input BAM file is sorted.


Helpful Tips
^^^^^^^^^^^^

The output BAM file can be sorted easily by pipeline the output to "samtools sort". For example::

	produse clip -i inFile.bam -o - | samtools sort > outFile.bam

Additional Info
^^^^^^^^^^^^^^^

Any overlap between read pairs is determined solely using the alignments in the BAM file: No realignment is performed.
If the bases at a position disagrees between a read pair, the base with the highest mapping quality is used. In the case of
a tie, the base from the read which starts later is used.

When obtaining a consensus between INDELs, the read with the lowest number of INDELs in the overlapping region is used as the consensus.
No realignment is performed.

The consensus overlap is assigned to either read 1 or read 2. For the other read, the positions coresponding to the consensus are soft clipped.

.. warning:: Overlap between read pairs is removed via soft-clipping. This may cause problems for some structural variant callers which examine soft-clipped bases

