Trim
============

Purpose
^^^^^^^

Trims adapter sequence from barcoded reads, and saves this information in the read name.
Reads which do not fall within barcode range and mismatch threshold are discarded.

Run Using
^^^^^^^^^

::

    produse trim

or

::

    python /path/to/ProDuSe/ProDuSe/trim.py

Parameters
^^^^^^^^^^

    :-c, --config:
        A configuration file which can provide any of the following arguments
    :-i --input:
        A fastq file to be trimmed. This argument must be specified twice. This file may be gzipped.
    :-o --output:
        | An output file for writing trimmed reads. This argument must be specified exactly twice.
        | These files can automatically be gzipped by appending '.gz' to the output file name.
        | The output file order corresponds with the input file order.
    :-as --adapter_sequence:
        The barcode sequence range, described using IUPAC bases.
    :-ap --adapter_position:
        | The positions in the adapter sequence to consider when comparing reference (i.e. -as) and actual adapter sequences
        | 0 = Do not consider this position
        | 1 = Consider this position
        | Note that the entire adapter sequence will be trimmed, even if a position is labeled 0.
    :--mm --max_mismatch:
        The maximum number of mismatches allowed between the reference and actual adapter sequences.
    :-v:
        | Instead of discarding reads with mismatching adapter sequences and saving all other read, do the opposite (save mismatching reads, discard all others).
        | Useful for debugging.
    :-u:
        Do not trim the adapter sequence.


Helpful tips
^^^^^^^^^^^^

If the read discard rate is extremely high, check the supplied adapter sequence.

Additional Notes
^^^^^^^^^^^^^^^^

If the supplied fastqs are multiplexed, a single sample can be extracted if no other samples use barcoded adapters, or if the barcoded adapter sequences between are distinct enough (i.e. the difference between the two barcode sequences exceeds the maximum mismatch threshold).
Note that there may be some spillover if non-barcoded reads start with a sequence that falls within the barcode sequence range by chance, or if the differences between barcoded sequences is only slightly higher than the maximum mismatch threshold.
