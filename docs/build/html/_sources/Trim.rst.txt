Trim
============

Purpose
^^^^^^^

Trims the barcode sequence from each read, and stores the barcode in a FASTQ tag
Reads where the barcode does not fall within specified barcode range and mismatch threshold are discarded.

Run Using
^^^^^^^^^

::

    produse trim

or

::

    python /path/to/ProDuSe/ProDuSe/Trim.py

Parameters
^^^^^^^^^^

    :-c, --config:
        A configuration file which can provide any of the following arguments. See the `config page`_ for more details.
    :-i --input:
        Two FASTQ files to be trimmed. These files may be gzipped.
    :-o --output:
        | Two output FASTQ files
        | These files can automatically be gzipped by appending '.gz' to the output file name.
        | The output file order corresponds with the input file order (i.e -i read1.fastq read2.fastq -> -o read1.out.fastq read2.out.fastq)
    :-b --barcode_sequence:
        The barcode sequence range, described using IUPAC bases. Can be determined using `adapter_predict`_
    :-p --barcode_position:
        | The positions in the adapter sequence to consider when comparing reference (i.e. -b) and actual adapter sequences
        | 0 = Do not consider this position
        | 1 = Consider this position
        | Note that the entire barcode sequence will be trimmed, even if a position is labeled 0.
    :--mm --max_mismatch:
        The maximum number of mismatches allowed between the reference and actual adapter sequences before a read pair is discarded
    :---reverse:
        | Instead of discarding reads with mismatching barcode sequences, do the opposite (save mismatching reads, discard all others).
        | Useful for debugging.
    :--no_trim:
        Do not trim the adapter sequence. Only store the adapter sequence a read tag.
    :--trim_other_end:
        | Examine the other end of the read for barcode sequences as well.
        | Will not remove partial barcodes
        | Useful when read lengths are significantly longer than the median fragment length

    .. _config page: Config_Files.html
    .. _adapter_predict: adapter_predict.html

Helpful Tips
^^^^^^^^^^^^

If the read discard rate is extremely high, check the supplied barcode sequence.

Additional Notes
^^^^^^^^^^^^^^^^

If the supplied fastqs are multiplexed, a single sample can be extracted if no other samples use barcoded adapters, or if the barcoded adapter sequences between are distinct enough (i.e. the difference between the two barcode sequences exceeds the maximum mismatch threshold).
Note that there may be some spillover if non-barcoded reads start with a sequence that falls within the barcode sequence range by chance, or if the differences between barcoded sequences is only slightly higher than the maximum mismatch threshold.
