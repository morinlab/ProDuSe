#!/bin/usr/env python

#
#
#
#
#

import argparse
import os
import nucleotide
import fastq
import qsub
import random
import subprocess
import sys

if __name__ == '__main__':

    desc = "Convert paired-end fastq files with variable adapters into sub-fastq's containing same variable adapter sequence"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-pe", "--paired_end",
        nargs=2,
        required=True,
        help="Takes the forward and reverse fastq files from paired end sequencing"
        )
    parser.add_argument(
        "-as", "--adapter_sequence",
        default=0,
        help="Specify the contiguous arrangement of nucleotides in the adapter sequence using any symbol defined by the IUPAC"
        )
    parser.add_argument(
        "-mm", "--max_mismatch",
        type=int, 
        required=False,
        help="The maximum number of mismatching nucleotides to see in the adapter and reads. Default is the length of the adapter (i.e. do not discard any reads)"
        )
    parser.add_argument(
        "-fp", "--file_prefix",
        required=False,
        default='',
        help="The prefix to include on all the sub-fastq's, if not specified the files will be named 1.NNNNN.fastq or 2.NNNNN.fastq, where NNNNN is the respected variable adapter sequence"
        )
    parser.add_argument(
        "-o", "--output_directory",
        default="./",
        help="Specify and output directory, otherwise the current working directory will be used"
        )
    parser.add_argument(
        "-rpa", "--record_parse_amount",
        type=int,
        default=100000,
        help="The number of records to read before writing to sub-fastq files"
        )
    parser.add_argument(
        "-t", "--number_of_threads",
        type=int,
        default=1,
        help="Specify the number of threads to each alignment, mark duplicates, and the final merge"
        )
    parser.add_argument(
        "-r", "--reference",
        required=False,
        help="Specify the reference genome to align each set of sub-fastqs to"
        )
    parser.add_argument(
        "-rd", "--remove_duplicates", 
        action="store_true", 
        default=False,
        help="After alignment, should PCR duplicates be removed or flagged"
        )
    parser.add_argument(
        "-fb", "--fix_barcode",
        action="store_true",
        default=True,
        help="Attempt to determine if a unique sequence still exists even though there are mismatches in barcode"
        )
    parser.add_argument(
        "-s", "--sge",
        action="store_true", 
        default=False,
        help="If cluster services are available, jobs will be submitted using qsub"
        )
    parser.add_argument(
        "-c", "--clean", 
        action="store_true", 
        default=False,
        help="Remove any temporary files (e.g. sub-bams and sub-fastqs)"
        )
    parser.add_argument(
        "-F", "--fastq_preprocess",
        action="store_true",
        default=False,
        help="If specified, only trim barcode sequence from fastq (skip other stages in pipeline)"
        )
    args = parser.parse_args()


    ### INPUT CHECK ################################################################################

    #if not os.path.isfile(args.paired_end[0]):
        #sys.exit("The forward fastq file does not exist")

    #if not os.path.isfile(args.paired_end[1]):
        #sys.exit("The reverse fastq file does not exist")

    #if not os.path.isfile(args.reference):
        #sys.exit("The reference file does not exist")

    #if not os.path.exists(args.output_directory):
        #sys.exit("The output directory specified does not exist")

    if len(args.adapter_sequence) == 0:
        sys.exit("The adapter sequence must contain at least a single IUPAC nucleotide symbol")

    if not nucleotide.is_iupac(args.adapter_sequence):
        sys.exit("The adapter sequence must be defined using IUPAC uppercase symbols")

    if nucleotide.is_unambiguous(args.adapter_sequence):
        sys.exit("The adapter sequence you provided is unambiguous, there is no need for this method")
    if not args.reference:
        if not args.fastq_preprocess:
            sys.exit("The reference sequence must be provided")

    ### DETERMINE POSSIBLE ADAPTERS ###
    ref_adapter = ''.join([args.adapter_sequence, args.adapter_sequence])
    good_adapters = nucleotide.make_unambiguous(ref_adapter)


    ### CREATE TEMPORARY DIRECTORY AND FILES FOR THIS RUN ##########################################
    
    # Define the file names for all possible output/input
    forward_output_file = 'forward.fastq'
    reverse_output_file = 'reverse.fastq'
    bwa_sam_file = 'bwa.sam'
    read_group_bam_file = 'read_group.bam'
    picard_bam_file = 'picard.bam'
    picard_metric_file = 'metric_file.txt'
    filtered_bam_file = 'filtered.bam'

    # Append the output prefix to these file names if specified as a 
    # command line argument
    if not args.file_prefix == '':
        forward_output_file = '.'.join([args.file_prefix, forward_output_file])
        reverse_output_file = '.'.join([args.file_prefix, reverse_output_file])
        bwa_sam_file = '.'.join([args.file_prefix, bwa_sam_file])
        read_group_bam_file = '.'.join([args.file_prefix, read_group_bam_file])
        picard_bam_file = '.'.join([args.file_prefix, picard_bam_file])
        picard_metric_file = '.'.join([args.file_prefix, picard_metric_file])
        filtered_bam_file ='.'.join([args.file_prefix, filtered_bam_file])

    # Make a new unique temporary directory
    prev_tmp_dir = '/'.join([args.output_directory, 'tmp'])
    tmp_dir = prev_tmp_dir
    while os.path.exists(tmp_dir):
        tmp_dir = ''.join([prev_tmp_dir, str(random.randint(1000000000, 9999999999))])
    os.makedirs(tmp_dir)
    print ':'.join(['PROCESS-FASTQ\tCreated Temporary Directory', tmp_dir[3:]])


    ### PROCESS FASTQ INPUT / OUTPUT ###############################################################
    print ':'.join(['PROCESS-FASTQ\tTrimming Adapter Sequence', args.adapter_sequence])


    # Open Fastq files for reading and writing
    forward_input = fastq.FastqOpen(args.paired_end[0], 'r')
    reverse_input = fastq.FastqOpen(args.paired_end[1], 'r')
    forward_output = fastq.FastqOpen('/'.join([tmp_dir, forward_output_file]), 'w')
    reverse_output = fastq.FastqOpen('/'.join([tmp_dir, reverse_output_file]), 'w')

    # For each read in the forward and reverse fastq files, trim the adapter, at it to the read id
    # and output this new fastq record to the temprary output fastq
    count = 0
    for forward_read in forward_input:
        count = count + 1
        reverse_read = reverse_input.next()
        forward_adapter = forward_read.trim(len(args.adapter_sequence))
        reverse_adapter = reverse_read.trim(len(args.adapter_sequence))
        forward_read.id = ''.join(['@', forward_adapter.seq, reverse_adapter.seq, forward_read.id[1:]])
        reverse_read.id = ''.join(['@', forward_adapter.seq, reverse_adapter.seq, reverse_read.id[1:]])
        forward_output.next(forward_read)
        reverse_output.next(reverse_read)
        if count % 100000 == 0:
            print ':'.join(['PROCESS-FASTQ\tProcessed Reads', str(count)])

    print ':'.join(['PROCESS-FASTQ\tProcessed Reads', str(count)])
    # Close the Fastq Files (especially since these files buffer the records), want to make sure
    # BWA has all the records
    forward_input.close()
    reverse_input.close()
    forward_output.close()
    reverse_output.close()

    if args.fastq_preprocess:
        sys.exit(0)

    ### SETUP SGE ##################################################################################

    # These variables will only be used if args.sge is true

    sge = qsub.Qsub()

    ### FASTQ PRE-PROCESSING #######################################################################





    ### BWA ########################################################################################

    commands = [
        'bwa', 'mem',
        '-M', 
        '-t', str(args.number_of_threads),
        args.reference,
        '/'.join([tmp_dir, forward_output_file]),
        '/'.join([tmp_dir, reverse_output_file]),
        '>', '/'.join([tmp_dir, bwa_sam_file])
        ]

    if args.sge:
        sge.change(job_name='PROCESS-FASTQ:BWA')
        sge.run(commands)
    else:
        subprocess.call(commands)


    ### ADD READGROUPS #############################################################################

    commands = [
        'python', 'add_readgroups.py',
        '--input', '/'.join([tmp_dir, bwa_sam_file]),
        '--output', '/'.join([tmp_dir, read_group_bam_file]),
        '--adapter_sequence', args.adapter_sequence
        ]

    if args.sge:
        sge.change(job_name='PROCESS-FASTQ:ADD_READ_GROUP')
        sge.run(commands)
    else:
        subprocess.call(commands)


    ### PICARD #####################################################################################

    commands = [
        'java', '-Xmx4g',
        '-jar', '$PICARDROOT/MarkDuplicates.jar',
        ''.join(["INPUT=", '/'.join([tmp_dir, bwa_sam_file])]),
        ''.join(["OUTPUT=", '/'.join([tmp_dir, picard_bam_file])]),
        ''.join(["METRIC_FILE=", '/'.join([tmp_dir, picard_metric_file])]),
        'REMOVE_DUPLICATES=True'
        'ASSUME_SORTED=False',
        'SORT_ORDER=coordinate'
        ]

    if args.sge:
        sge.change(job_name='PROCESS-FASTQ:PICARD')
        sge.run(commands)
    else:
        subprocess.call(commands)


    ### REMARK DUPLICATES ##########################################################################






    ### CLEAN UP ###################################################################################









