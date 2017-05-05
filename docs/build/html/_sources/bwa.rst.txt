bwa
===

Maps supplied reads to a reference genome using the Burrows-Wheeler aligner. Sorts and converts the resulting SAM file into a sorted BAM file. Note that, for all intents and purposes, this script is simply a python wrapper around bwa, with some functionality restrictions.

Run Using
^^^^^^^^^

::

    produse align

or

::

    python /path/to/ProDuSe/ProDuSe/bwa.py

Parameters
^^^^^^^^^^

    :-c --config:
        A configuration file which can supply any of the parameters described below
    :-i --input:
        A fastq file containing the reads to be mapped. This option must be provided exactly twice. These files may be gzipped.
    :-o --output:
        Path and name of the output BAM file.
    :-t --threads:
        Number of threads to use while running BWA.
    :-r --reference:
        Reference genome, in fasta format. BWA indexes should exist in the same parent directory. If these do not exist, use `bwa index` to generate them.
    :-R --readgroup:
        Read group information to add to the BAM file header. Useful for running some exterior tools, such as Picard.

Additional Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^

To keep the terminal as clean as possible, many bwa output messages are suppressed. This can inadvertently include some error messages. If bwa.py exits abruptly, try running bwa directly to determine the cause of the failure.

