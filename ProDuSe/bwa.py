import configargparse
import sys
import os
import subprocess
import time
import itertools

desc = "Align Paired Reads for the PRODUSE Pipeline"
parser = configargparse.ArgParser(description=desc)
parser.add(
    "-c", "--config",
    required=False,
    is_config_file=True,
    help="An optional configuration file for any of the input arguments."
    )
parser.add(
    "-i", "--input",
    required=True,
    action="append",
    type=str,
    help="A pair of fastq files for reading processed by trim.py"
    )
parser.add(
    "-o", "--output",
    required=True,
    type=str,
    help="An output bam file for writing"
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

    if not len(args.input) == 2:
        parser.error('--input must be a pair (i.e. a sized two list) of fastq files')

    print_prefix = "PRODUSE-BWA        " ;
    sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Starting...\n")    

    if not os.path.isfile(args.input[0]):
        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Error: input fastq files does not exist - " + args.input[0] + "\n")
        sys.exit(1)

    if not os.path.isfile(args.input[1]):
        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Error: input fastq files does not exist - " + args.input[1] + "\n")
        sys.exit(1)       

    output_prefix = ""
    resulting_file = ""
    if args.output[-4:] == ".bam":
        output_prefix = args.output[:-4]
        resulting_file = args.output
    else:
        output_prefix = args.output
        resulting_file = args.output + ".bam"

    if os.path.isfile(resulting_file) :
        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Error: output bam file already exists - " + resulting_file + "\n")
        sys.exit(1)

    command1 = [
        'bwa', 'mem',
        args.reference,
        args.input[0],
        args.input[1]
        ]

    command2 = [
        'samtools', 'view',
        '-Sb', "-"
        ]

    command3 = [
        'samtools', 'sort',
        "-", output_prefix
        ]

    ps1 = subprocess.Popen(command1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ps2 = subprocess.Popen(command2, stdout=subprocess.PIPE, stdin=ps1.stdout, stderr=subprocess.PIPE)
    ps3 = subprocess.Popen(command3, stdin=ps2.stdout, stderr=subprocess.PIPE)

    ps1_iterator = iter(ps1.stderr.readline, b"")
    ps2_iterator = iter(ps2.stderr.readline, b"")
    ps3_iterator = iter(ps3.stderr.readline, b"")

    counter = 0 
    for line in itertools.chain(ps1_iterator, ps2_iterator, ps3_iterator):
        if line.startswith("[M::mem_process_seqs]"):
            counter += int(line.split(" ")[2])
            sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Reads Processed:" + str(counter) + "\n")

    ps1.stdout.close()
    ps2.stdout.close()
    ps3.wait()

if __name__ == '__main__':
    main()
