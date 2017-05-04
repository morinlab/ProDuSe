Quick Start
===========

Check and install dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To run ProDuSe successfully, ensure the following non-python dependencies are installed on your system:

    - Samtools_ (>=1.3.1)
    - `Burrows-Wheeler Aligner`_ (>=0.7.12)
    - Python_ (>=2.7)
    - Stitcher_ (>=5.1.3)

    .. _Samtools: http://samtools.sourceforge.net/
    .. _Burrows-Wheeler Aligner: http://bio-bwa.sourceforge.net/
    .. _Python: https://www.python.org/
    .. _Stitcher: https://github.com/Illumina/Pisces

(Optional) Install ProDuSe
^^^^^^^^^^^^^^^^^^^^^^^^^^

To run ProDuSe as a command line command, install ProDuSe as follows:

From the command line::

    cd /path/to/ProDuSe/clone/
    python setup.py install

This will automatically install all python dependencies

Setup ProDuSe configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run ProDuSe, you will need the following:

 1. Sequence (fastq) files from your sample of interest
 2. The barcode sequence of your sample
 3. `ProDuse configuration file`_

 .. _ProDuse configuration file: https://github.com/morinlab/ProDuSe/blob/master/etc/produse_config.ini

To start, modify the ProDuSe configuration file with your sample's information::

    [snv]
    min_molecules = 2
    mutant_molecules = 3
    variant_allele_fraction_threshold = 0.01
    min_reads_per_uid = 1

    [adapter_qc]
    adapter_sequence = NNNWSMRWSYWKMWWT
    adapter_position = 0001111111111110
    max_mismatch = 3

    [collapse_bwa]

    [collapse]
    strand_position = 0001111111111110
    adapter_max_mismatch = 3
    duplex_position = 0000000001111110
    duplex_max_mismatch = 2

    [trim_bwa]

    [trim]
    adapter_sequence = NNNWSMRWSYWKMWWT
    adapter_position = 0001111111111110
    max_mismatch = 3

    [filter_produse]

See the configuration help page for information on what each of these options specify.

Run ProDuSe
^^^^^^^^^^^

Finally, run ProDuSe with the specified information

If you installed ProDuSe::

    produse run_produse -x /path/to/Stitcher.exe -r /path/to/reference/genome/build -c /path/to/produse/config/file -f /path/to/fastq.R1 /path/to/fastq.R2

If you did not install ProDuSe::

    /path/to/ProDuSe/ProDuSe/ProdusePipleine.py -x /path/to/Stitcher.exe -r /path/to/reference/genome/build -c /path/to/produse/config/file -f /path/to/fastq.R1 /path/to/fastq.R2

All results will be placed in the current working directory (this can be changed using `-d`), in the folder `produse_analysis_directory`

Running multiple samples at once
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ProDuSe can process multiple samples sequentially with the same command. You will need a `sample configuration` file to use this feature

.. _sample configuration: https://github.com/morinlab/ProDuSe/blob/master/etc/sample_config.ini

Specify each sample in this configuration file::

    [Sample1]
    fastqs=/path/to/Sample1/fastq.R1,/path/to/Sample1/fastq.R2

    [Sample2]
    fastqs=/path/to/Sample2/fastq.R1,/path/to/Sample2/fastq.R2

.. note:: There is no space after the comma

If you would like to specify different parameters for each sample, they can be specified in the sample_config file::

    [Sample1]
    fastqs=/path/to/Sample1/fastq.R1,/path/to/Sample1/fastq.R2
    adapter_sequence=NNNWSMRWSYWKMWW
    adapter_position=000111111111111

    [Sample2]
    fastqs=/path/to/Sample2/fastq.R1,/path/to/Sample2/fastq.R2
    max_mismatch=2

These parameters will be used for that sample alone

.. note:: Parameters specified in the sample configuration file override ProDuSe configuration file paramters for that sample, while ProDuSe configuration file parameters override command line parameters for ALL samples

When running ProDuSe, this sample configuration file can be specified instead of fastq files

If ProDuSe is installed::

   produse run_produse -x /path/to/Stitcher.exe -r /path/to/reference/genome/build -c /path/to/produse/config/file *-sc /path/to/sample/configuration/file*

If ProDuSe is not installed::

    /path/to/ProDuSe/ProDuSe/ProdusePipleine.py -x /path/to/Stitcher.exe -r /path/to/reference/genome/build -c /path/to/produse/config/file *-sc /path/to/sample/configuration/file*

All results will be outputted in individual directories under 'produse_analysis_directory'



