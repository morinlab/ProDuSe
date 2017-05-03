Config
==============

Configures output directories, input and output file paths for each step of the ProDuSe pipeline, and configuration files for each script

produse_config.py performs the following tasks:

    - Creates a main output directory (produse_analysis_direcetory) in the specified directory
    - Creates a log file listing the supplied arguments in the main output directory
    - Creates a Makefile in the main output directory (Depreciated, has no purpose)
    - (If necessary) Creates normal and bwa indexes for the designated reference genome

For each sample supplied:

    - Creates a sample-specific directory and subdirectories in the main output directory using the sample's name
    - Symbolically Links input files into this sample-specific directory
    - Creates a config file for each step in the ProDuSe pipeline, containing sample-specific parameters

Run Using
^^^^^^^^^

::

    produse configure_produse
    /path/to/ProDuSe/ProDuSe/configure_produse.py

Parameters
^^^^^^^^^^

    :-f --fastqs:
        | If running a single sample, the fastq files for the sample to be analyzed.
        | Two arguments must be supplied.
        | A single sample directory names 'Sample' will be created.
        | Mutually exclusive with -sc --sample_config.
    :-sc --sample_config:
        | A configuration file listing sample names, fastq locations, and sample-specific parameters.
        | A directory will be created for each sample, using the name supplied.
        | Mutually exclusive with -f --fastqs.
    :-d --output_directory:
        The directory to output ProDuSe results. Default is the current working directory.
    :-r --reference:
        Reference genome build, in fasta format.
    :-c --config:
        A configuration file listing paramters to be supplied to each stage of the ProDuSe pipeline

Additional Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^

While indexes will be automatically generated if none are present in the reference genome directory, this make take a significant amount of time
