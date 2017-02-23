#! /usr/bin/env python

# USAGE:
#   See configure_produse.py -h for details
#
# DESCRIPTION:
#   Creates directories and config files for running the ProDuSe pipeline
#
# AUTHORS:
#   Creator: Marco Albuquerque
#   Christopher Rushton (ckrushto@sfu.ca)


import argparse
import configparser
import os
import sys
import platform


def make_directory(sample_dir, fastqs, sampleConfig, pconfig, reference):
    """
        Creates subdirectories, symlinks, and config files for ProDuSe

        Creates subdirectories to store intermediate and final results, symlinks for required
        input and files for each script in the ProDuSe pipeline, and sets up config files for each
        script. Sample config paramaters superseed produce config parameters

        Args:
            sample_dir: Path to the output directory for ProDuSe
            fasqs: Array listing paired fastq files, either compressed or uncompressed
            sampleConfig: Dictionary listing the parameters specified in the sample config file
            pconfig: Path to ProDuSe configurations file
            reference: Reference genome, in FASTA format


    """

    # Creates ProDuSe subdirectories
    tmp_dir = os.sep.join([sample_dir, "tmp"])
    os.makedirs(tmp_dir)

    config_dir = os.sep.join([sample_dir, "config"])
    os.makedirs(config_dir)

    figures_dir = os.sep.join([sample_dir, "figures"])
    os.makedirs(figures_dir)

    data_dir = os.sep.join([sample_dir, "data"])
    os.makedirs(data_dir)

    results_dir = os.sep.join([sample_dir, "results"])
    os.makedirs(results_dir)

    log_dir = os.sep.join([sample_dir, "logs"])
    os.makedirs(log_dir)

    # Import produse config file
    cparser = configparser.RawConfigParser()
    cparser.read(pconfig)

    # Checks to ensure fastq files exist
    if not os.path.isfile(fastqs[0]):
        sys.stderr.write("ERROR: Unable to locate %s\n" % (fastqs[0]))
        sys.exit(1)

    if not os.path.isfile(fastqs[1]):
        sys.stderr.write("ERROR: Unable to locate %s\n" % (fastqs[1]))
        sys.exit(1)

    # Sets up fastq files
    raw_fastq1 = os.sep.join([data_dir, "raw_R1.fastq"])
    raw_fastq2 = os.sep.join([data_dir, "raw_R2.fastq"])
    if os.path.splitext(fastqs[0])[1] == ".gz" and os.path.splitext(fastqs[1])[1] == ".gz":
        raw_fastq1 = '.'.join([raw_fastq1, "gz"])
        raw_fastq2 = '.'.join([raw_fastq2, "gz"])

    elif os.path.splitext(fastqs[0])[1] == ".gz" or os.path.splitext(fastqs[1])[1] == ".gz":
        sys.stderr.write("ERROR: Fastq files must be either both gziped or both decompressed\n")
        sys.exit(1)

    os.symlink(os.path.abspath(fastqs[0]), raw_fastq1)
    os.symlink(os.path.abspath(fastqs[1]), raw_fastq2)

    # Sets up symlinks for intermediate  ProDuSe files
    trim_one_tmp = os.sep.join([tmp_dir, "trim_R1.fastq.gz"])
    trim_two_tmp = os.sep.join([tmp_dir, "trim_R2.fastq.gz"])
    trim_one_data = os.sep.join([data_dir, "trim_R1.fastq.gz"])
    trim_two_data = os.sep.join([data_dir, "trim_R2.fastq.gz"])

    os.symlink(trim_one_tmp, trim_one_data)
    os.symlink(trim_two_tmp, trim_two_data)

    # Creates Trim config file
    new_config = configparser.RawConfigParser()
    new_config.add_section("config")
    new_config.set("config", "input", ''.join(["[", ','.join([raw_fastq1, raw_fastq2]), "]"]))
    new_config.set("config", "output", ''.join(["[", ','.join([trim_one_tmp, trim_two_tmp]), "]"]))
    for (key, val) in cparser.items("trim"):
        if key in sampleConfig:
            val = sampleConfig[key]
        new_config.set("config", key, val)
    output = open(os.sep.join([config_dir, "trim_task.ini"]), 'w')
    new_config.write(output)
    output.close()

    # Creates BWA output symlinks
    trim_tmp = os.sep.join([tmp_dir, "trim.bam"])
    trim_data = os.sep.join([data_dir, "trim.bam"])
    os.symlink(trim_tmp, trim_data)

    # Creates BWA trim config file
    new_config = configparser.RawConfigParser()
    new_config.add_section("config")
    new_config.set("config", "input", ''.join(["[", ','.join([trim_one_data, trim_two_data]), "]"]))
    new_config.set("config", "output", trim_tmp)
    new_config.set("config", "reference", reference)
    for (key, val) in cparser.items("trim_bwa"):
        new_config.set("config", key, val)
        if key in sampleConfig:
            val = sampleConfig[key]
    output = open(os.sep.join([config_dir, "trim_bwa_task.ini"]), 'w')
    new_config.write(output)
    output.close()

    # Creates collapse output symlinks
    collapse_one_tmp = os.sep.join([tmp_dir, "collapse_R1.fastq.gz"])
    collapse_two_tmp = os.sep.join([tmp_dir, "collapse_R2.fastq.gz"])
    collapse_one_data = os.sep.join([data_dir, "collapse_R1.fastq.gz"])
    collapse_two_data = os.sep.join([data_dir, "collapse_R2.fastq.gz"])

    os.symlink(collapse_one_tmp, collapse_one_data)
    os.symlink(collapse_two_tmp, collapse_two_data)

    # Creates collapse config file
    new_config = configparser.RawConfigParser()
    new_config.add_section("config")
    new_config.set("config", "input", trim_data)
    new_config.set("config", "output", ''.join(["[", ','.join([collapse_one_tmp, collapse_two_tmp]), "]"]))
    for (key, val) in cparser.items("collapse"):
        new_config.set("config", key, val)
        if key in sampleConfig:
            val = sampleConfig[key]
    output = open(os.sep.join([config_dir, "collapse_task.ini"]), 'w')
    new_config.write(output)
    output.close()

    # Creates bwa collapse output symlinks
    collapse_tmp = os.sep.join([tmp_dir, "collapse.bam"])
    collapse_data = os.sep.join([data_dir, "collapse.bam"])

    os.symlink(collapse_tmp, collapse_data)

    # Creates bwa collapse config file
    new_config = configparser.RawConfigParser()
    new_config.add_section("config")
    new_config.set("config", "input", ''.join(["[", ','.join([collapse_one_data, collapse_two_data]), "]"]))
    new_config.set("config", "output", collapse_tmp)
    new_config.set("config", "reference", reference)
    for (key, val) in cparser.items("collapse_bwa"):
        new_config.set("config", key, val)
        if key in sampleConfig:
            val = sampleConfig[key]
    output = open(os.sep.join([config_dir, "collapse_bwa_task.ini"]), 'w')
    new_config.write(output)
    output.close()

    # Creates SNV calling config file
    new_config = configparser.RawConfigParser()
    new_config.add_section("config")
    new_config.set("config", "input", results_dir + os.sep + "SplitMerge.sorted.bam")
    new_config.set("config", "output", os.sep.join([results_dir, "variants.vcf"]))
    new_config.set("config", "reference", reference)
    for (key, val) in cparser.items("snv"):
        new_config.set("config", key, val)
        if key in sampleConfig:
            val = sampleConfig[key]
    output = open(os.sep.join([sample_dir, "config", os.sep, "snv_task.ini"]), 'w')
    new_config.write(output)


def make_makefile(produse_directory, analysis_dir, samples):
    """
        Creates a Make file, listing the location of the various scripts that will be run, their order, and parameters

        Args:
            produse_directory: Directory containing ProDuSe scripts
            analyis_directory: Output directory for intermediate and final results of the produse pipeline
            samples: Samples listed in the optional sample_config file
    """
    with open(os.sep.join([analysis_dir, "Makefile"]), "w") as make:
        make.write(''.join(["produse_dir :=", produse_directory, '\n']))
        make.write(''.join(["analysis_dir :=", analysis_dir, '\n']))
        make.write("\n")
        make.write("trim_script := $(produse_dir)/trim.py\n")
        make.write("bwa_script := $(produse_dir)/bwa.py\n")
        make.write("collapse_script := $(produse_dir)/collapse.py\n")
        make.write("adapter_qc_script := $(produse_dir)/adapter_qc.py\n")
        make.write("snv_script := $(produse_dir)/snv.py\n")
        make.write("\n")
        make.write("".join(["SAMPLES = ", " ".join(samples)]))
        make.write("\n")
        make.write("SNV_TASK = $(addprefix snv_task-, $(SAMPLES))")
        make.write("\n\n")
        make.write("all: $(SNV_TASK)\n")
        make.write("\n")
        make.write("snv_task-%: collapse_index_bwa_task-%\n")
        make.write("\tpython $(snv_script) --config=$(analysis_dir)/$*/config/snv_task.ini && touch $@\n")
        make.write("\n")
        make.write("collapse_index_bwa_task-%: collapse_bwa_task-%\n")
        make.write("\tsamtools index $(analysis_dir)/$*/data/collapse.bam && touch %@\n")
        make.write("\n")
        make.write("collapse_bwa_task-%: collapse_task-%\n")
        make.write("\tpython $(bwa_script) --config=$(analysis_dir)/$*/config/collapse_bwa_task.ini && touch $@\n")
        make.write("\n")
        make.write("adapter_qc_task-%: collapse_task-%\n")
        make.write("\tpython $(adapter_qc_script) --config=$(analysis_dir)/$*/config/adapter_qc_task.ini && touch $@\n")
        make.write("\n")
        make.write("collapse_task-%: trim_bwa_task-%\n")
        make.write("\tpython $(collapse_script) --config=$(analysis_dir)/$*/config/collapse_task.ini && touch $@\n")
        make.write("\n")
        make.write("trim_bwa_task-%: trim_task-%\n")
        make.write("\tpython $(bwa_script) --config=$(analysis_dir)/$*/config/trim_bwa_task.ini && touch $@\n")
        make.write("\n")
        make.write("trim_task-%:\n")
        make.write("\tpython $(trim_script) --config=$(analysis_dir)/$*/config/trim_task.ini && touch $@\n")
        make.write("\n")
        make.close()


"""
Processes command line arguments

Returns:
    args: A namespace object listing the provided command line arguments
Raises:
    parser.error(): THrown if neither fastq files nor sampleconfig files are provided
"""

# Process command line arguments
desc = "Process Duplex Sequencing Data"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument(
    "-f", "--fastqs",
    nargs=2,
    type=lambda x: is_valid_file(x, parser),
    help="Forward and reverse fastq files from paired end sequencing"
    )
parser.add_argument(
    "-o", "--output_directory",
    default="." + os.sep,
        help="All ProDuSe results and subdirectories will be placed in this directory [Default: %(default)s]"
    )
parser.add_argument(
    "-r", "--reference",
    required=True,
    type=lambda x: is_valid_file(x, parser),
    help="Reference genome, in FASTA format"
    )
parser.add_argument(
    "-c", "--config",
    required=True,
    help="ProDuSe config file, listing adapter sequences, positions, and other paramters to be passed to pipeline scripts"
    )
parser.add_argument(
    "-sc", "--sample_config",
    required=False,
    help="Config file listing sample names and FASTQ locations, used to analyze multiple samples"
    )

def is_valid_file(file, parser):
    """
    Checks to ensure the specified file exists, and throws an error if it does not

    Args:
        file: A filepath
        parser: An argparse.ArgumentParser() object. Used to throw an exception if the file does not exist

    Returns:
        type: The file itself

    Raises:
        parser.error(): An ArgumentParser.error() object, thrown if the file does not exist
    """
    if os.path.isfile(file):
        return file
    else:
        parser.error("The file %s does not exist" % (file))


def main(args=None):

    """
    Performs setup of ProDuSe config files and subdirectories

    Args:
        args: A namespace object listing ProDuSe parameters. See get_args() for a list of supported options

    """

    if args is None:
        args = parser.parse_args()
    # If neither -f nor -sc was specified, throw an error
    if not args.sample_config and not args.fastqs:
        raise parser.error("Either --fastqs or --sample_config must be provided")

    # Sets up ProDuSe output directory
    output_directory = os.path.abspath(os.sep.join([args.output_directory, "produse_analysis_directory"]))
    produse_directory = os.path.dirname(os.path.realpath(__file__))

    # If the output directory already exists, throw an error
    if os.path.isdir(output_directory):
        sys.stderr.write("ERROR: %s already exists\n" % output_directory)
        sys.exit(1)
    else:
        os.makedirs(output_directory)

    # Reads sections and key-value pairs from config file
    sample_config = configparser.RawConfigParser()
    sample_config.read(args.sample_config)

    # Creates a MakeFile
    make_makefile(produse_directory, output_directory, sample_config.sections())

    # Create a sample direectory and config files for each sample listed i the sample config file
    for sample in sample_config.sections():

        sample_dir = os.sep.join([output_directory, sample])
        os.makedirs(sample_dir)

        sampleDict = dict(sample_config.items(sample))

        # Obtains fastq files from either the command line or sampleconfig file
        if "fastqs" in sampleDict:

            fastqs = sampleDict["fastqs"].split(",")
            if len(fastqs) != 2:
                sys.stderr.write("ERROR: Exactly two fastq files must be list in the sample config file, comma seperated\n")
                sys.stderr.write("Example: fastqs=path/foward.fq.gz,path/reverse.fq.gz")
                exit(1)

        elif not args.fastqs:
            sys.stderr.write("ERROR: \'fastqs\' are not specified in %s, nor were they provided in the arguments\n" % (args.sample_config))
            sys.stderr.write("Fastqs can be specified using \'-f\'\n")
            exit(1)
        else:
            fastqs = args.fastqs

        # Setup sample directory
        make_directory(sample_dir, fastqs, sampleDict, args.config, args.reference)

    print("Configuration complete")


if __name__ == '__main__':

    main()
