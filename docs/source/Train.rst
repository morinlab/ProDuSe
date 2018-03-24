Train
=====

Description
^^^^^^^^^^^

Trains ProDuSe's variant calling filter using one or more sets of validated variants. A random forest
classifier will be trained using this data

Run Using
^^^^^^^^^

::

	produse train

or

::

	/path/to/ProDuSe/Clone/ProDuSe/Train.py

Parameters
^^^^^^^^^^

	:-c --config:
		A configuration file which can provide any of the following arguments. See the `config page`_ for more details.
	:-b --bam:
		One or more post-clipoverlap BAM files, which were used for variant filtering and validations. Must be specified in an order that coresponds to the order specified by -v/--validations
	:-v --validations:
		One or more VCF files listing validated variants. The number and order of file must corespond to files specified in -b/--bam.
	:-o --output:
		Output file which will store the random forest classifier
	:-r --reference:
		Reference genome, in FASTA format. An index should also be present in the same directory.
	:-t --targets:
		Optional. One or more BED files specifying regions in which to restrict candidate variant identification. Must be specified in an order coresponding to the order specified in -b/--bam.

.. _config page: Config_Files.html

Additional Information
^^^^^^^^^^^^^^^^^^^^^^

For each sample provided, ProDuSe will identify all candidate variants in a manner identical to Call.
If a candidate variant is flagged in the validation VCF file with VALIDATED=TRUE, it will be considered
a true variant, while variants flagged with VALIDATED=UNKNOWN will be ignored completely. All remaining
variants will be flaged as artifacts. False variants will be randomly subset so there are the same number
of false and true variants prior to training the variant filter.
