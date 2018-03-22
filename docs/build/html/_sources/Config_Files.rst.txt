Configuration Files
===================

Most ProDuSe commands allow arguments to be specified using a configuration INI file. In most cases, this is a much more convinient
method of specifying arguments if ProDuSe is to be run multiple times. While you can view more detailed info on `config files here`_, 
a brief overview will be described here.

.. _config files here: http://www.voidspace.org.uk/python/configobj.html#the-config-file-format


Standard Config Files
^^^^^^^^^^^^^^^^^^^^^

All config files provided to produse commands (using -c/--config) use the following format::

	[Script_Name]
	argument=parameter
	argument2=parameter2

Where Script_Name is the name of the tool which will use those arguments, with the arguments listed in the section below. Short form argument
names are not supported.

Note that config file sections and arguments not relevent to a given tool will safely be ignored, so it is possible to merge config files::

	[Trim]
	barcode_sequence=NNNWTWYYT
	max_mismatch=3

	[Pipeline]
	reference=genome.fa

A demo sample configuration file can be obtained here_

.. _here: https://github.com/morinlab/ProDuSe/blob/master/etc/produse_config.ini

Sample Configuration Files
^^^^^^^^^^^^^^^^^^^^^^^^^^

When running multiple samples uing `run_produse`_, a special sample configuration file is specified instead::

	[Sample_Name]
	fastqs=/some/path/read.R1.fastq.gz,/some/path/read.R2.fastq.gz

	[Sample_Name_2]
	fastqs=/some/other/path/read.R1.fastq.gz,/some/other/path/read.R2.fastq.gz

Where Sample_Name is the name used for that sample (this is prepended to all intermediate and output files). Any parameters supported by `run_produse`_
can be specified, and will be used for that sample **only**. For instance::

	[Sample_Name_1]
	fastqs=/some/path/read.R1.fastq.gz,/some/path/read.R2.fastq.gz
	reference=GRCh37.fa

	[Sample_Name_2]
	fastqs=/some/other/path/read.R1.fastq.gz,/some/other/path/read.R2.fastq.gz
	reference=GRCh38.fa

GRCh37 will be used as the reference genome for Sample 1, while GRCh38 will be used for sample 2. `Click here to view a demo sample config file`_

.. _Click here to view a demo sample config file: https://github.com/morinlab/ProDuSe/blob/master/etc/sample_config.ini
.. _run_produse: run_produse.html