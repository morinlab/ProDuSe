import argparse
import ConfigParser
import os
import sys


def make_directory(sample_dir, fastqs, config, reference):


    configparser = ConfigParser.RawConfigParser()
    configparser.read(config)

    tmp_dir = '/'.join([sample_dir, "tmp"])
    os.makedirs(tmp_dir)

    config_dir = '/'.join([sample_dir, "config"])
    os.makedirs(config_dir)

    figures_dir = '/'.join([sample_dir, "figures"])
    os.makedirs(figures_dir)

    data_dir = '/'.join([sample_dir, "data"])
    os.makedirs(data_dir)

    results_dir = '/'.join([sample_dir, "results"])
    os.makedirs(results_dir)

    log_dir = '/'.join([sample_dir, "logs"])
    os.makedirs(log_dir)

    if not os.path.isfile(fastqs[0]):
        print >> sys.stderr, ''.join(["This file does not exist: ", fastqs[0]])
        sys.exit(1)

    if not os.path.isfile(fastqs[1]):
        print >> sys.stderr, ''.join(["This file does not exist: ", fastqs[1]])
        sys.exit(1)
    
    raw_one = '/'.join([data_dir, "raw_R1.fastq"])
    raw_two = '/'.join([data_dir, "raw_R2.fastq"])
    if fastqs[0][-3:] == ".gz" and fastqs[1][-3:] == ".gz":
        raw_one = '.'.join([raw_one, "gz"])
        raw_two = '.'.join([raw_two, "gz"])

    elif fastqs[0][-3:] == ".gz" or fastqs[1][-3:] == ".gz":
        print >> sys.stderr, ''.join(["Fastq files must be both gziped or both decompressed: ", sample])
        sys.exit(1)

    os.symlink(fastqs[0], raw_one)
    os.symlink(fastqs[1], raw_two)

    trim_one_tmp = '/'.join([tmp_dir, "trim_R1.fastq.gz"])
    trim_two_tmp = '/'.join([tmp_dir, "trim_R2.fastq.gz"])
    trim_one_data = '/'.join([data_dir, "trim_R1.fastq.gz"])
    trim_two_data = '/'.join([data_dir, "trim_R2.fastq.gz"])

    os.symlink(trim_one_tmp, trim_one_data)
    os.symlink(trim_two_tmp, trim_two_data)

    ### MAKE TRIM CONFIG
  
    new_config = ConfigParser.RawConfigParser()
    new_config.add_section("config")
    new_config.set("config", "input", ''.join(["[",','.join([raw_one, raw_two]),"]"]))
    new_config.set("config", "output", ''.join(["[",','.join([trim_one_tmp, trim_two_tmp]),"]"]))
    for (key,val) in configparser.items("trim"):
        new_config.set("config", key, val)
    output = open('/'.join([config_dir, "trim_task.ini"]) , 'w')
    new_config.write(output)
    output.close()

    ### MAKE TRIM BWA I/O

    trim_tmp = '/'.join([tmp_dir, "trim.bam"])
    trim_data = '/'.join([data_dir, "trim.bam"])

    os.symlink(trim_tmp, trim_data)

    ### MAKE TRIM BWA CONFIG

    new_config = ConfigParser.RawConfigParser()
    new_config.add_section("config")
    new_config.set("config", "input", ''.join(["[",','.join([trim_one_data, trim_two_data]),"]"]))
    new_config.set("config", "output", trim_tmp)
    new_config.set("config", "reference", reference)
    for (key,val) in configparser.items("trim_bwa"):
        new_config.set("config", key, val)
    output = open('/'.join([config_dir, "trim_bwa_task.ini"]) , 'w')
    new_config.write(output)
    output.close()

    ### MAKE COLLAPSE I/O

    collapse_one_tmp = '/'.join([tmp_dir, "collapse_R1.fastq.gz"])
    collapse_two_tmp = '/'.join([tmp_dir, "collapse_R2.fastq.gz"])
    collapse_one_data = '/'.join([data_dir, "collapse_R1.fastq.gz"])
    collapse_two_data = '/'.join([data_dir, "collapse_R2.fastq.gz"])

    os.symlink(collapse_one_tmp, collapse_one_data)
    os.symlink(collapse_two_tmp, collapse_two_data)

    ### MAKE COLLAPSE CONFIG

    new_config = ConfigParser.RawConfigParser()
    new_config.add_section("config")
    new_config.set("config", "input", trim_data)
    new_config.set("config", "output", ''.join(["[",','.join([collapse_one_tmp, collapse_two_tmp]),"]"]))
    for (key,val) in configparser.items("collapse"):
        new_config.set("config", key, val)
    output = open('/'.join([config_dir, "collapse_task.ini"]) , 'w')
    new_config.write(output)
    output.close()    

    ### MAKE COLLAPSE BWA I/O

    collapse_tmp = '/'.join([tmp_dir, "collapse.bam"])
    collapse_data = '/'.join([data_dir, "collapse.bam"])

    os.symlink(collapse_tmp, collapse_data)

    ### MAKE COLLAPSE BWA CONFIG

    new_config = ConfigParser.RawConfigParser()
    new_config.add_section("config")
    new_config.set("config", "input", ''.join(["[",','.join([collapse_one_data, collapse_two_data]),"]"]))
    new_config.set("config", "output", collapse_tmp)
    new_config.set("config", "reference", reference)
    for (key,val) in configparser.items("collapse_bwa"):
        new_config.set("config", key, val)
    output = open('/'.join([config_dir, "collapse_bwa_task.ini"]) , 'w')
    new_config.write(output)
    output.close()

    ### MAKE SNV CONFIG

    new_config = ConfigParser.RawConfigParser()
    new_config.add_section("config")
    new_config.set("config", "input", collapse_data )
    new_config.set("config", "output", "/".join([results_dir, "variants.txt"]) )
    new_config.set("config", "reference", reference)
    for (key,val) in configparser.items("snv"):
        new_config.set("config", key, val)
    output = open('/'.join([sample_dir, "config/snv_task.ini"]), 'w')
    new_config.write(output)   

    return 0

def make_makefile(produse_directory, analysis_dir, samples):
    make = open('/'.join([analysis_dir, "Makefile"]), "w")
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
    make.write("snv_task-%: collapse_bwa_task-%\n")
    make.write("\tpython $(snv_script) --config=$(analysis_dir)/$*/config/snv_task.ini && touch $@\n")
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
    return 0

if __name__ == '__main__':

    desc = "Process Duplex Sequencing Data"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-f", "--fastqs",
        nargs=2,
        required=False,
        help="Takes the file locations of the forward and reverse fastq files from paired end sequencing"
        )
    parser.add_argument(
        "-o", "--output_directory",
        required=True,
        help="Takes the file locations to store the trimmed forward and reverse fastq files"
        )
    parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Takes the file location of the reference genome"
        )
    parser.add_argument(
        "-c", "--config",
        required=True,
        help="Optional config if other inputs not defined"
        )
    parser.add_argument(
        "-sc", "--sample_config",
        required=False,
        help="Optional for setting up directory for several samples"
        )
    args = parser.parse_args()

    output_directory = os.path.abspath('/'.join([args.output_directory, "produse_analysis_directory"]))
    produse_directory = os.path.dirname(os.path.realpath(__file__))

    if os.path.isdir(output_directory):
        print >> sys.stderr, ''.join(["This directory already exists: ", output_directory])
        sys.exit(1)
    else:
        os.makedirs(output_directory)

    if not args.sample_config:
        if not args.fastqs:
            parser.error('without --sample_config, you must specify --fastqs')
    
    else:

        sample_config = ConfigParser.RawConfigParser()
        sample_config.read(args.sample_config)

        make_makefile(produse_directory, output_directory, sample_config.sections())

        for sample in sample_config.sections():

            sample_dir = '/'.join([output_directory, sample])
            os.makedirs(sample_dir)

            for (key,val) in sample_config.items(sample):

                ### MAKE SAMPLE DIRECTORY

                make_directory(sample_dir, val.split(","), args.config, args.reference)