# ProDuSe
An analysis pipeline, helper scripts and Python classes to **Pro**cess **Du**plex **Se**quence data

## Description


## Installation 

### Dependencies

You will need to install the following tools before installing the ProDuSe package:

* `python==2.7`
* `bwa==0.7.12`
* `samtools==1.3.1`

To install the ProDuSe package run the following command:
### Currently being updated

```bash
cd path/to/github/clone/ProDuSe
python setup.py install
```

## Running ProDuSe

### The Analysis Pipeline

You will first need to retrieve two configuration files:

#### `produse_config.ini`
 * command line arguments for each stage in the analysis pipeline
 * retrieve a sample config.ini file [here](https://github.com/morinlab/ProDuSe/blob/master/etc/produse_config.ini)

#### `sample_config.ini`
 * paired fastq files for all samples you wish to run the analysis pipeline on
 * retrieve a sample sample_config.ini file [here](https://github.com/morinlab/ProDuSe/blob/master/etc/sample_config.ini)

To run the analysis pipeline you simply need to run the following command:

```produse run_produse
    -c config.ini
    -sc sample_config.ini
    -r /path/to/ref.fa
    -o /path/to/output/directory
```

This will run the entire ProDuSe pipeline on all samples specified in the sample_config.ini file
Results will be located in the following directory:

```bash
cd /path/to/output/produse_analysis_directory
```

### Helper Scripts

The ProDuSe package includes a variety of helper scripts to aid in the analysis of duplex sequencing data.

All scripts included in the current package can be found by running the following:

```bash
produse -h
```

#### produse adapter_predict

If you need to confirm the expected adapter sequence of a sample you should run the following command:

```bash
produse adapter_predict -i input1.fastq input2.fastq
```

This tool will print a predicted adapter sequence based off of ACGT abundances at each position. It uses these observed abundances and finds the closest expected abundance for an IUPAC unambiguous or ambiguous base.

### Python Classes

Two major python classes are included with ProDuSe. 

#### The Alignment Class

The first is the alignment class. This linearly processes reads from a BAM file until both read pairs have been identified, at which point the first yield to the developer occurs.

#### The Position Class

The second if the position class. This class aims to create a duplex sequencing ready mpileup class.

Full descriptions of two python classes can be retrieved here
