Call
====

Identifies all possible positions which support an alternate allele (no matter how weakly) across the entire capture space. These variants are then
filtered using a random forest filter, based upon numerous characteristics


Run Using
^^^^^^^^^

::

    produse call

or

::

    python /path/to/ProDuSe/ProDuSe/Call.py

Parameters
^^^^^^^^^^

    :-c --config:
        An optional configuration file which can specify any of the arguments described below
    :-i --input:
        An input BAM containing collapsed (and ideally clipped) reads. This file must be sorted by position.
    :-o --output:
        Output VCF file containing filtered variants
    :-u, --unfilted:
        Output VCF file containing ALL possible variants
    :-r --reference:
        Reference genome, in FASTA format. A reference index should be present in the same directory
    :-t --target_bed:
        A BED3 file specifying a capture space to restrict variant calling
    :-f, --filter:
        A piclke containing a tranined Random Forest Classifier



