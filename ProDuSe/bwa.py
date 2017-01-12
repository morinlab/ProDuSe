import configargparse
import sys
import os
import subprocess
import time
import itertools
from distutils.version import StrictVersion

desc = "Align Paired Reads for the PRODUSE Pipeline"
parser = configargparse.ArgParser(description=desc)
parser.add(
    "-c", "--config",
    required=False,
    is_config_file=True,
    help="An optional configuration file for any of the input arguments."
    )
parser.add(
        "-t", "--threads",
        required=False,
        default=1,
        help="how many threads to allow bwa to use"
        )
parser.add(
    "-i", "--input",
    required=True,
    action="append",
    type=str,
    help="A set of fastq files for reading processed by trim.py. The first two should be un-merged R1 and R2 (in that order) that did not get merged by FLASH and the third file should be the FLASH-merged reads"
    )
parser.add(
    "-o", "--output",
    required=True,
    type=str,
    help="A prefix for your bam files for writing initial alignments and merged bam"
    )
parser.add(
    "-r", "--reference",
    required=True,
    type=str,
    help="A genome reference file with BWA indexes"
    )


def main(args=None):
    
    if args == None:
        args = parser.parse_args()

    if not len(args.input) == 3:
        parser.error('--input must be trio of fastq files from FLASH')

    print_prefix = "PRODUSE-BWA        " ;
    sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Starting...\n")    

    if not os.path.isfile(args.input[0]) or not os.path.isfile(args.input[1]) or not os.path.isfile(args.input[2]):
        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Error: one or more of input fastq files does not exist - " + args.input[0] + args.input[1] + args.input[2]+ "\n")
        sys.exit(1)


    output_format = "BAM"
    #this should not be optional
    

    if os.path.isfile(args.output) :
        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Error: output file already exists - " + args.output + "\n")
        sys.exit(1)

    # Determine samtools version and make sure it is 1.3.1 or higher
    version_samtools_command1 = [
        'samtools'
        ]
    
    ps1 = subprocess.Popen(version_samtools_command1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps1_iterator = iter(ps1.stderr.readline, b"")
    for line in ps1_iterator:
        if line.startswith("Version"):
            if StrictVersion(line.split(" ")[1]) < StrictVersion("1.3.1"):
                sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Error: samtools version must be 1.3.1 or higher - " + line.split(" ")[1] + "\n")
                sys.exit(1)
            break
    #run BWA for merged reads
    threads = "-t " + str(args.threads)
    command1 = [
        'bwa', 'mem', threads,
        args.reference,
        args.input[2]
        ]

    command2 = [
        'samtools', 'sort',
        '--output-fmt', output_format,
        '-'
        ]
    out1 = args.output + ".flash_reads.bam"
    out2 = args.output + ".unmerged_reads.bam"
    out_merge = args.output + ".merge.bam"
    out1_fh = open(out1,"w")
    out2_fh = open(out2,"w")
    #output_fh = open(args.output, 'w')
    print "running"
    print command1
    ps1 = subprocess.Popen(command1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps2 = subprocess.Popen(command2, stdout=out1_fh, stdin=ps1.stdout, stderr=subprocess.PIPE)

    ps1_iterator = iter(ps1.stderr.readline, b"")
    ps2_iterator = iter(ps2.stderr.readline, b"")

    counter = 0 
    for line in itertools.chain(ps1_iterator, ps2_iterator):
        if line.startswith("[M::mem_process_seqs]"):
            counter += int(line.split(" ")[2])
            sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Reads Processed:" + str(counter) + "\n")

    ps1.stdout.close()
    ps2.wait()
    command3 = [
        'bwa','mem',threads,
        args.reference,
        args.input[0],
        args.input[1]
        ]
    ps1 = subprocess.Popen(command3, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps2 = subprocess.Popen(command2, stdout=out2_fh, stdin=ps1.stdout, stderr=subprocess.PIPE)
    ps1_iterator = iter(ps1.stderr.readline, b"")
    ps2_iterator = iter(ps2.stderr.readline, b"")

    counter = 0
    for line in itertools.chain(ps1_iterator, ps2_iterator):
        if line.startswith("[M::mem_process_seqs]"):
            counter += int(line.split(" ")[2])
            sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Reads Processed:" + str(counter) + "\n")

    ps1.stdout.close()
    ps2.wait()
    out1_fh.close()
    out2_fh.close()
    
    merge_command = [
        'samtools', 'merge', out_merge, 
        out1, out2]
    ps1 = subprocess.Popen(merge_command)
    ps1.wait()
    index_command = ['samtools', 'index', out_merge]
    ps2 = subprocess.Popen(index_command)
    ps2.wait()
if __name__ == '__main__':
    main()
