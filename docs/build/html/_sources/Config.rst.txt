Config
======

Configures output directories, input and output file paths for each step of the ProDuSe pipeline, and configuration files for each script.

produse_config.py performs the following tasks:

    - Creates a main output directory (produse_analysis_directory) in the specified directory
        - If a directory named 'produse_analysis_directory' is specified as the base output directory, it will be used directly as the main output directory
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

or

::

    python /path/to/ProDuSe/ProDuSe/configure_produse.py

Parameters
^^^^^^^^^^

    :-f --fastqs:
        | If running a single sample, the fastq files for the sample to be analyzed.
        | Two arguments must be supplied.
        | A sample-specific output directory will be created, using the base name of the input files.
        | Mutually exclusive with -sc --sample_config.
    :-sc --sample_config:
        | A configuration file listing sample names, fastq locations, and sample-specific parameters.
        | A directory will be created for each sample, using the name supplied.
        | Mutually exclusive with -f --fastqs, -n --normal.
    :-d --output_directory:
        The directory to output ProDuSe results. Default is the current working directory.
    :-r --reference:
        | Reference genome build, in fasta format.
        | If bwa indexes do not exist in the same directory, they will be generated automatically.
    :-c --config:
        A configuration file listing parameters to be supplied to each stage of the ProDuSe pipeline.
    :-n --normal:
        | A BAM file or two fastq files originating from a matched normal for -f --fastqs. 
        | The input type is automatically determined by the file extension.
        | Mutually exclusive with -sc --sample_config.
    :-ns --normal_adapter_sequence:
        | An optional parameter representing the barcode sequence of -n --normal FASTQ files.
        | Only valid if fastq files are provided for -n --normal.
        | If this parameter is provided, -np --normal_adapter_position must be provided as well.
        | If not provided, no trim configuration file will be created for the -n --normal fastq files, and ProdusePipeline.py will not trim these fastqs.
    :-np --normal_adapter_position:
        | An optional parameter representing the barcode position of -n --normal FASTQ files.
        | Only valid if fastq files are provided for -n --normal.
        | If this parameter is provided, -ns --normal_adapter_sequence must be provided as well.
        | If not provided, no trim configuration file will be created for the -n --normal fastq files, and ProdusePipeline.py will not trim these fastqs.
    :--append_to_directory:
        | Configure sample-specific output directories inside an existing produse_analysis_directory.
        | If a sample-specific output directory already exists inside produse_analysis_directory, that sample will be skipped, and the conflicting directory will remain unmodified.


Additional Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^

While indexes will be automatically generated if none are present in the reference genome directory, this make take a significant amount of time.
