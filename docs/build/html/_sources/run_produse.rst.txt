Run_ProDuSe
===========

Purpose
^^^^^^^

Wrapper script which runs each step of the ProDuSe Pipeline sequentially on each sample. Can process multiple samples simultaneously.

Run Using
^^^^^^^^^

::

	produse run_produse

or

::

	/path/to/ProDuSe/ProdusePipeline.py

Parameters
^^^^^^^^^^

General Parameters for running analysis

    :-c, --config:
        An INI configuration file, which can specify any of the following 
        arguments. Any arguments passed at the command line override those specified in the configuration file, which override defaults.
    :--fastqs:
    	A pair of filepaths to a set of paired-end FASTQ files. Mutually 
    	exclusive with -sc\--sample_config.
    :-sc, --sample_config:
    	An INI configuration file, specifying one or more samples to analyze. This file can also provide sample-specific arguments. 
    :-d, --outdir:
    	Base output directory for intermediate files and results. 
    	--directory_name will be created inside this directory. Default is the 
    	current working directory.
    :-r, --reference:
    	A filepath to a reference genome FASTA file. Ideally, BWA and .fai 
    	indexes should be present in the same directory. If they are not, they 
    	will be generated automatically.
    :-j, --jobs:
    	Number of samples to analyze in parallel. Use 0 or a negative number to 
    	process as many samples as possible simultaneously. Default is 1.

Additional Analysis Parameters

	:--bwa:
		A filepath to a bwa executable. Default is $(which bwa).
	:--samtools:
		A filepath to a samtools executable. Default is $(which samtools)
	:--directory_name:
		Name of the directory to create inside -d/--outdir to store intermediate files and results. Default is "produse_analysis_directory".
	:--append_to_directory:
		If --directory_name already exists inside -d/--outdir, place the intermediate files and results for this analysis inside this directory. If any samples have the same name as those inside --directory_name, they will not be analyzed.

Barcode Trimming Parameters

	:-b, --barcode_sequence:
		The sequence of the degenerate barcode, specified in IUPAC bases. 
	:-p, --barcode_position:
		Positions in the -b/--barcode sequence to use when comparing expected and actual barcode sequences (1=Yes, 0=No). The entire barcode will still be trimmed, regardless of which positions are flagged with 0.
	:-mm, --max_mismatch:
		Maximum number of positions (specified using -p/--barcode_position) which can fall outside the expected degeneracy range before the read pair is discarded. This is the total between both the forward and reverse read.
	:--trim_other_end:
		Examine the end of the read for the presence of a barcode (for example, in the case of read-through). Will not remove partial barcodes

Family Collapsing Parameters

	:-fm, --family_mask:
		Barcode positions to consider when determining if two reads are members of the same family (1=Yes, 0=No).
	:-fmm, --family_mismatch:
		Maximum number of mismatched positions (specified using -fm/--family_mask) permitted before two reads are considered members of different families.
	:-dm, --duplex_mask:
		Barcode positions to consider when determining if two families are in duplex (1=Yes, 0=No).
	:-dmm, --duplex_max_mismatch:
		Maximum number of mismatched positions (specified using -dm/--duplex_mask) permitted before two families are not considered a duplex
	:-t, --targets:
		A filepath to a BED file listing the capture regions of interest. Read pairs that do not overlap these positions will be discarded
	:--tag_family_members:
		Store the name of each read incorporated into a family in the read tag "Zm"

Filtering Parameters

	:-f, --filter:
		A filepath to a pickled Random Forest Classifier, used to filter variants. Can be generated using 'produse train'


.. note:: If no -b/--barcode is specified for one or more samples, adapter_predict will be run on those samples automatically

Pipeline Stages
^^^^^^^^^^^^^^^

This command will run the following stages of the ProDuSe Pipeline sequentially on each sample:

	- adapter_predict (Estimates the degenerate barcode sequence, only run if no barcode is specified)
	- Trim (Trim Barcodes)
	- Align reads to reference (Burrows-Wheeler Aligner)
	- Collapse (Identify and merge duplicate reads into a consensus)
	- ClipOverlap (Identifies positions which overlap between read pairs, and generates a consensus)
	- Call (Flag and filter variants)

Configuration Files
^^^^^^^^^^^^^^^^^^^

In lieu of specifying all arguments at the command line, arguments can be specified in **either** the sample configuration file,
or the main produse configuration file. This approach is recommended, as it improves reproducibility between runs. 

If a single argument is specified multiple times, the following orders of precedence apply: 

1. -sc/--sample_config file
2. Command line
3. -c/--config file

Arguments specified in the -sc/--sample_config file only apply to the specified sample.
More information on configuration files can be viewed here.

Directory Layout
^^^^^^^^^^^^^^^^

When run_produse is called, it creates a directory structure for each sample inside of --directory_name, as follows::

	directory_name
		Sample_Name
			config
			tmp
			results
		ProDuSe_Task.log

The contents of each folder are as follows:

	:config:
		Stores configuration files for each step of the pipeline, and files indicating when a pipeline stage is complete
	:tmp:
		Stores intermediate files (ex. Trimmed FASTQ files, raw BAM files)
	:results:
		Stores the final BAM file and variant calls

All parameters used to run a given instance, as well as software versions, are specified in ProDuSe_Task.log


