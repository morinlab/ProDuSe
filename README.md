# ProDuSe
An analysis pipeline, helper scripts and Python classes to **Pro**cess **Du**plex **Se**quence data

## Description


## Installation 

### Dependencies

You will need to install the following tools before installing the ProDuSe package:

* `python==2.7`
* `bwa==0.7.12`
* `samtools==0.1.19`

To install the ProDuSe package run the following command:

```bash
pip install ProDuSe
```

## Running ProDuSe

### The Analysis Pipeline

To run the analysis pipeline you simply need to run the following command:

```bash
produse analysis_pipeline 
    -c config.ini
    -sc sample_config.ini
    -r /path/to/ref.fa
    -o /path/to/output
```

The command required two configuration files:

#### `config.ini`
 * command line arguments for each stage in the analysis pipeline
 * retrieve a sample config.ini file [here](https://github.com/morinlab/ProDuSe/blob/master/etc/produse_config.ini)

#### `sample_config.ini`
 * paired fastq files for all samples you wish to run the analysis pipeline on
 * retrieve a sample sample_config.ini file [here](https://github.com/morinlab/ProDuSe/blob/master/etc/sample_config.ini)

