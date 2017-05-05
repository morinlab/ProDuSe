SplitMerge
==========

Processes each read in the Stitched BAM file. Unstitched reads are simply saved to the output BAM file. If a read was stitched, the bases unique to the original forward and reverse read are
identified and assigned to a replacement forward and reverse read. Bases shared between the two sequences are assigned to only one of the paired reads, to remove overlap.
In addition, if Stitcher misaligned a portion of the Stitched sequence, the correct start position of each read is determined, and an INDEL is added as necessary.

Run Using
^^^^^^^^^

::

    produse split_merge

or

::

    python /path/to/ProDuSe/ProDuSe/ProDuSe/SplitMerge.py

Parameters
^^^^^^^^^^

    :-c --config:
        An optional configuration file which can provide any of the arguments described below.
    :-i --input:
        A BAM file containing reads merged (stitched) by Stitcher. This file must be sorted by read name.
    :-u --unstitched_input:
        A BAM file containing reads immediately before they were merged by stitcher. If running the ProDuSe pipeline, this will be the collapsed BAM file. This file must
        be sorted by read name.
    :-o --output:
        Path and name of the output BAM file.

.. note:: Input files must be sorted by read name, not position.


Additional Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^

While the vast majority (>99.99%) of stitched reads will be successfully separated, some reads with very complicated stitching will be discarded. Note that while this step
will identify some additional INDELs that bwa did not identify due to their size, ProDuSe does not currently call INDELs.


