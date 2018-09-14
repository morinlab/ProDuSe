[![Build Status](https://travis-ci.org/morinlab/ProDuSe.svg?branch=master)](https://travis-ci.org/morinlab/ProDuSe)

# ProDuSe
A variant caller designed for use with libraries generated using semi-degenerate barcoded adapters.

## Description

See the full wiki page for more information: http://produse.readthedocs.io/en/latest/

## Installation 

### Dependencies

You will need to install the following tools before installing the ProDuSe package:

* `python>=3.4`
* `bwa>=0.7.0`
* `samtools>=1.3.1`

ProDuSe will check that a valid version of these tools are installed prior to running the pipeline

To install the ProDuSe package run the following command:

```bash
cd path/to/github/clone/ProDuSe
python setup.py install
```
All required python dependencies will be installed during this step

## Running ProDuSe

### The Analysis Pipeline: Very Quick Start

You can view more detailed instructions on the [wiki](http://produse.readthedocs.io/en/latest/)

All parameters required to run ProDuSe can be viewed using the following:
```
    produse run_produse -h
```

Alternatively, if you wish to run ProDuSe without installing it, you can run `ProdusePipeline.py` manually in a similar manner:
```
    /path/to/ProDuSe/ProdusePipeline.py -h
```

While these parameters can be specified individually, they can also be provided using a configuration file

To run the analysis pipeline you simply need to run the following command:
```
    produse run_produse
    -c /path/to/github/clone/etc/produse_config.ini
```

Alternatively:
```
    /path/to/ProDuSe/ProdusePipeline.py 
    -c /path/to/github/clone/etc/produse_config.ini
```

This will run the entire ProDuSe pipeline on all samples specified in the sample_config.ini file, which can be found in 
etc/sample_config.ini

Results will be located in the following directory:

```bash
ls ./produse_analysis_directory
```

### Helper Scripts

The ProDuSe package is comprised of several stages to aid in the analysis of duplex sequencing data.

These stages can be be viewed by running the following:

```bash
produse -h
```

#### produse adapter_predict

If you need to confirm the expected adapter sequence of a sample you should run the following command:

```bash
produse adapter_predict -i input1.fastq input2.fastq
```

This tool will print a predicted adapter sequence based off of ACGT abundances at each position. It uses these observed abundances and finds the closest expected abundance for an IUPAC unambiguous or ambiguous base.

