Quick Start
===========

Check and install dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To run ProDuSe successfully, ensure the following non-python dependencies are installed on your system:

    - Samtools_ (>=1.3.1)
    - `Burrows-Wheeler Aligner`_ (>=0.7.0)
    - Python_ (>=3.4)

    .. _Samtools: http://samtools.sourceforge.net/
    .. _Burrows-Wheeler Aligner: http://bio-bwa.sourceforge.net/
    .. _Python: https://www.python.org/

(Optional) Install ProDuSe
^^^^^^^^^^^^^^^^^^^^^^^^^^

To run ProDuSe as a command line command, install ProDuSe as follows:

From the command line::

    cd /path/to/ProDuSe/clone/
    python setup.py install

This will automatically install all python dependencies.

Setup ProDuSe configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run ProDuSe, you will need the following:

 1. Sequence (fastq) files from your sample of interest
 3. `ProDuse configuration file`_

 .. _ProDuse configuration file: https://github.com/morinlab/ProDuSe/blob/master/etc/produse_config.ini

If you are unsure of the barcode sequence, it can be estimated using `adapter_predict`_

.. _adapter_predict: Adapter Predict.rst

To start, modify the ProDuSe configuration file with your sample's information::

    [Pipeline]
    barcode_sequence=NNNWSMRWSYWKMWWT
    barcode_position=0001111111111110
    max_mismatch=3
    family_mask=0001111111111110
    family_mismatch=3
    duplex_mask=0000000001111110
    duplex_mismatch=2
    reference=test/reference/GRCh38_TNFSF14.fa
    sample_config=etc/sample_config.ini
    filter=etc/default_filter.p

See run_produse_ for a detailed description of each argument

.. _run_produse: run_produse.html

Run ProDuSe
^^^^^^^^^^^

Finally, run ProDuSe with the specified information.

If you installed ProDuSe::

    produse run_produse -c /path/to/produse/config/file --fastqs /path/to/fastq.R1 /path/to/fastq.R2

If you did not install ProDuSe::

    /path/to/ProDuSe/ProDuSe/ProdusePipleine.py -c /path/to/produse/config/file --fastqs /path/to/fastq.R1 /path/to/fastq.R2

All results will be placed in the current working directory (this can be changed using `-d`), under `produse_analysis_directory`.

Running multiple samples at once
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ProDuSe can process multiple samples sequentially with the same command. You will need a `sample configuration` file to use this feature.

.. _sample configuration: https://github.com/morinlab/ProDuSe/blob/master/etc/sample_config.ini

Specify each sample in this configuration file::

    [Sample1]
    fastqs=/path/to/Sample1/fastq.R1,/path/to/Sample1/fastq.R2

    [Sample2]
    fastqs=/path/to/Sample2/fastq.R1,/path/to/Sample2/fastq.R2


If you would like to specify different parameters for each sample, they can be specified in the sample_config file::

    [Sample1]
    fastqs=/path/to/Sample1/fastq.R1,/path/to/Sample1/fastq.R2
    barcode_sequence=NNNWSMRWSYWKMWW
    barcode_position=000111111111111

    [Sample2]
    fastqs=/path/to/Sample2/fastq.R1,/path/to/Sample2/fastq.R2
    max_mismatch=2

Arguments specified in the sample config file will only be used for that sample

.. note:: Parameters specified in the sample configuration file take precidence for that sample

When running ProDuSe, this sample configuration file will be specified instead of fastq files.

If ProDuSe is installed::

   produse run_produse -c /path/to/produse/config/file *-sc /path/to/sample/configuration/file*

If ProDuSe is not installed::

    /path/to/ProDuSe/ProDuSe/ProdusePipleine.py -c /path/to/produse/config/file *-sc /path/to/sample/configuration/file*

All results will be outputted in individual sample directories under 'produse_analysis_directory'.
