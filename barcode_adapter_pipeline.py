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
        required=True,
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
        default="True",
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
    args = parser.parse_args()


    ### INPUT CHECK ###
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


    ### DETERMINE POSSIBLE ADAPTERS ###
    ref_adapter = ''.join([args.adapter_sequence, args.adapter_sequence])
    good_adapters = nucleotide.make_unambiguous(ref_adapter)
    
    forward_files = {}
    reverse_files = {}
    sam_files = {}
    bwa_job_id = {}
    forward_outputs = {}
    reverse_outputs = {}

    good_adapters.append('MISC')

    # Create file names for all good adapters as well as a miscellaneous file for mismatching sequences
    for adapter in good_adapters:
        forward_files[adapter] = '.'.join(['1',adapter,'fastq'])
        reverse_files[adapter] = '.'.join(['2',adapter,'fastq'])
        sam_files[adapter] = '.'.join([adapter, 'sam'])
        if not args.file_prefix == '':
            forward_files[adapter] = '.'.join([args.file_prefix, forward_files[adapter]])
            reverse_files[adapter] = '.'.join([args.file_prefix, reverse_files[adapter]])
            sam_files[adapter] = '.'.join([args.file_prefix, sam_files[adapter]])
        if args.sge:
            bwa_job_id[adapter] = '-'.join(['BWA', sam_files[adapter]])

    ### PROCESS FASTQ INPUT / OUTPUT ###

    forward_input = open(args.paired_end[0], 'r')
    reverse_input = open(args.paired_end[1], 'r')

    end_of_file = False
    while not end_of_file:

        forward_records = {}
        reverse_records = {}  

        for adapter in good_adapters:
            forward_records[adapter] = []
            reverse_records[adapter] = []

        forward_read = list(itertools.islice(forward_input, 4 * args.record_parse_amount))
        reverse_read = list(itertools.islice(reverse_input, 4 * args.record_parse_amount))
 
        if len(forward_read) < args.record_parse_amount:
            end_of_file = True

        for i in range(len(forward_read)):
            
            if i%4 == 0:
                adapter = ''.join([forward_read[i+1][:len(args.adapter_sequence)],reverse_read[i+1][:len(args.adapter_sequence)]])

                if nucleotide.is_match(adapter,ref_adapter):
                    forward_records[adapter].append(forward_read[i])
                    reverse_records[adapter].append(reverse_read[i])
                    forward_records[adapter].append(forward_read[i+1][len(args.adapter_sequence):])
                    reverse_records[adapter].append(reverse_read[i+1][len(args.adapter_sequence):])
                    forward_records[adapter].append(forward_read[i+2])
                    reverse_records[adapter].append(reverse_read[i+2])    
                    forward_records[adapter].append(forward_read[i+3][len(args.adapter_sequence):])
                    reverse_records[adapter].append(reverse_read[i+3][len(args.adapter_sequence):])                                   
                else:
                    forward_records['MISC'].append(forward_read[i])
                    reverse_records['MISC'].append(reverse_read[i])
                    forward_records['MISC'].append(forward_read[i+1][len(args.adapter_sequence):])
                    reverse_records['MISC'].append(reverse_read[i+1][len(args.adapter_sequence):])
                    forward_records['MISC'].append(forward_read[i+2])
                    reverse_records['MISC'].append(reverse_read[i+2])    
                    forward_records['MISC'].append(forward_read[i+3][len(args.adapter_sequence):])
                    reverse_records['MISC'].append(reverse_read[i+3][len(args.adapter_sequence):])  

            else:
                continue

        for adapter in good_adapters:
            handler = open('/'.join([args.output_directory, forward_files[adapter]]), 'a')
            handler.write(''.join(forward_records[adapter]))
            handler.close()
            handler = open('/'.join([args.output_directory, reverse_files[adapter]]), 'a')
            handler.write(''.join(reverse_records[adapter]))
            handler.close()

    ### BWA ###

    for adapter in good_adapters:
        
        print ':'.join(['BWA', adapter])

        commands = ' '.join([
            'bwa', 'mem',
            '-M', 
            '-t', str(args.number_of_threads),
            args.reference,
            '/'.join([args.output_directory, forward_files[adapter]]),
            '/'.join([args.output_directory, reverse_files[adapter]]),
            '>', '/'.join([args.output_directory, sam_files[adapter]])
            ])

        if args.sge:
            qsub_commands = ' '.join([
                'qsub', '-cwd',
                '-b', 'y',
                '-N', bwa_id[adapter],
                commands
                ])
        else:
            os.popen(commands, 'r')

    ### MARK DUPLICATES ###

    


    ### MERGE BAMS ####



    ### FILTER MISC BAM AND MERGE ###



    ### PROCESS FINAL BAM AND REMOVE READS AT SAME POSITION with different adapters that have an edit distance greater than default ###



    ### CLEAN UP ###

    #for i in 0:(len(forward_file)-1):
