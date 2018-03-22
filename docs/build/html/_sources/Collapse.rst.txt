Collapse
========

Identifies the start position, barcode sequence, and mapping strand for each read pair in the supplied BAM file. If both reads in one or more reads share the same start position, mapping strand, and barcode sequence (within mismatch tolerance), they are flagged as a "family", and merged into a single consensus sequence.
If family members disagree at a given position, the most common base is used as a consensus.  In the case of a tie, the base with the highest aggregated quality score across all family members is used. The quality of each base set to the highest quality base at that position.

Run Using
^^^^^^^^^

::

    produse collapse

or

::

    python /path/to/ProDuSe/ProDuSe/Collapse.py


Parameters
^^^^^^^^^^

    :-c --config:
        A configuration file which can supply any of the arguments below. See the `config page`_ for more details.
    :-i --input:
        Input SAM/CRAM/BAM file. Each read must contain a read tag which stores adapter sequences
    :-o --output:
        Output SAM/CRAM/BAM file containing collapsed reads (Use "-" for stdout). The output file format will be chosen based upon the supplied file extension, or if "-" is used, will be the same as the input file format. Will be unsorted.
    :-fm --family_mask:
        | Positions in the barcode sequence to use when comparing barcode sequences between reads which originate from the same parental strand.
        | 1=Use this position, 0=Do not use this position.
    :-dm --duplex_mask:
        | The positions in the adapter sequence to use when comparing adapter sequences for reads of opposing types (i.e. forward vs reverse reads).
        | 1=Use this position, 0=Do not use this position.
    :-fmm --family_max_mismatch:
        The maximum number of mismatches allowed between the barcode sequence of two read pairs before two read pairs are considered members of different families (See -fm).
    :-dmm --duplex_max_mismatch:
        The maximum number of mismatches allowed between the barcode sequence of two families before they are considered as not in duplex (See -dm).
    :-r --reference:
        Reference genome, in FASTA format. Must be the same genome version that the reads were aligned against.
    :-t --targets:
        A BED3 file or better listing regions of interest. Any read pairs which fall entirely outside these regions will be discarded
    :--tag_family_members:
    	Store the original name of all reads which were incorporated into a family in the read tag "Zm"

.. _config page: Config_Files.html

Additional Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^

The runtime of Collapse depends not only on the absolute number of reads, but the proportion of reads which are duplicates. BAM files with high duplicate rates will take significantly longer than BAM files with a lower duplicate rate.

Currently, this version of Collapse does not perform local realignment of soft-clipped regions.

