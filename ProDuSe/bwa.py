import argparse
import printer
import sys
import select
import os
import signal
import subprocess

if __name__ == '__main__':

    printer.general('BWA')

    desc = "Align Paired Reads for the PRODUSE Pipeline"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-f", "--fastqs",
        nargs=2,
        required=True,
        help="Input Fastqs"
        )
    parser.add_argument(
        "-b", "--bam",
        nargs=1,
        required=True,
        help="Output Sorted Bam"
        )
    parser.add_argument(
        "-r", "--reference",
        nargs=1,
        required=True,
        help="BWA Reference"
        )
    args = parser.parse_args()

    output_prefix = ""
    if args.bam[0][-4:] == ".bam":
        output_prefix = args.bam[0][:-4]
    else:
        output_prefix = args.bam[0]

    command1 = [
        'bwa', 'mem',
        args.reference[0],
        args.fastqs[0],
        args.fastqs[1]
        ]

    command2 = [
        'samtools', 'view',
        '-Sb', "-"
        ]

    command3 = [
        'samtools', 'sort',
        "-", output_prefix
        ]

    ps1 = subprocess.Popen(command1, stdout=subprocess.PIPE) #, stderr=subprocess.PIPE)
    ps2 = subprocess.Popen(command2, stdout=subprocess.PIPE, stdin=ps1.stdout) #, stderr=subprocess.PIPE)
    ps3 = subprocess.Popen(command3, stdin=ps2.stdout) #, stderr=subprocess.PIPE)

    #ps1_iterator = iter(ps1.stderr.readline, b"")
    #ps2_iterator = iter(ps2.stderr.readline, b"")
    #ps3_iterator = iter(ps3.stderr.readline, b"")

    #counter = 0
    #printer.general(" ".join(["Reads Aligned:", str(counter)]))

    #for line in ps1_iterator:
 
    #    if line.startswith("[M::main_mem]"):
    #        counter += int(line.split()[2])
    #        printer.general(" ".join(["Reads Aligned:", str(counter)]))

    ps1.stdout.close()
    ps2.stdout.close()
    ps3.wait()
