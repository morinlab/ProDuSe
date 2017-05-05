Collapse
========

Identifies the start position, adapter sequence, and mapping strand for each read in the supplied BAM file. If two or more reads share the same start position, mapping strand, and adapter sequence (within mismatch tolerance), they are merged into a single consensus sequence.
If there is a mismatch at a given position, the most common base is used as a consensus. The quality of each base set to the highest quality base at that position. If an individual read contains too many mismatches, it is discarded prior to collapsing.

Run Using
^^^^^^^^^

::

    produse collapse

or

::

    python /path/to/ProDuSe/ProDuSe/collapse.py


Parameters
^^^^^^^^^^

    :-c --config:
        A configuration file which can provide any of the parameters below
    :-i --input:
        Input BAM file. The name of each read must contain the adapter sequence.
    :-o --output:
        Path and name of output collapsed fastq files. This parameter must be specified exactly twice.
        The output fastqs can be gzipped automatically by appending ".gz" to the output name.
    :-sp --strand_position:
        The positions in the adapter sequence to use when comparison adapter  sequences for reads of the same type (i.e. between forward reads, or between reverse reads).
        1=Use this position, 0=Do not use this position.
    :-dp --duplex_position:
        The positions in the adapter sequence to use when comparing adapter sequences for reads of opposing types (i.e. forward vs reverse reads).
        1=Use this position, 0=Do not use this position.
    :-amm --adapter_max_mismatch:
        The maximum number of mismatches allowed between the expected and actual adapter sequences when comparing reads of the same type (See -sp).
    :-dmm --duplex_max_mismatch:
        The maximum number of mismatches allowed between the expected and actual adapter sequences when comparing molecules of different types (See -dp).
    :-smm --sequence_max_mismatch:
        The maximum number of mismatches allowed in an individual read before it is discarded. This threshold should be adjusted based upon read length.
    :-oo --original_output:
        Path and name of fastq files to write original (i.e. pre-collapse) reads. Reads exceeding mismatch thresholds will still be discarded.
        This option must be be specified exactly twice, or not at all. These fastqs can be gzipped automatically by appending ".gz" to the output name.
    :-sf --stats_file:
        Path and name of a text file to store collapsing statistics.

Additional Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^

The runtime of this script depends not only on the absolute number of reads, but the proportion of reads which are duplicates. BAM files with high duplicate rates will take significantly longer than BAM files with a lower duplicate rate.


