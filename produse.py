#/bin/usr/env python
import qsub
import argument_parser
import printer
import setup
import subprocess
import sys
import select
import os
import signal

if __name__ == '__main__':

    printer.general('START')

    ### GET ARGUMENTS ##############################################################################

    printer.general('ARGS')
    args = argument_parser.parse()
    argument_parser.check(args)

    ### SETUP ######################################################################################

    printer.general('SETUP')
    data = setup.Setup(args)

    ### TRIM FASTQ #################################################################################

    printer.general('TRIM')

    if data.trim:    
        
        commands = [
            'python', 'trim.py',
            '-i', data.forward, data.reverse,
            '-o', data.trimmed_forward, data.trimmed_reverse,
            '-mm', str(args.max_mismatch),
            '-as', args.adapter_sequence
            ]
   
        printer.command(commands)

        if args.sge:
            printer.general('Running on SGE')
            sge.change()
            sge.run(commands)
        else:
            printer.general('Running on Command Line')
            subprocess.call(commands)
    else:
        printer.general('Trimming Ignored')

    ### BWA MEM ####################################################################################

    printer.general('BWA')

    if data.trim_bam:
  
        commands = [
            'python', 'bwa.py',
            '-f', data.trimmed_forward, data.trimmed_reverse,
            '-b', data.trimmed_bam,
            '-r', args.reference
            ]

        printer.command(commands)

        if args.sge:
            printer.general('Running on SGE')
            sge.change(job_name='PRODUSE:TRIM-ALIGNMENT')
            sge.run(commands)

        else:
            printer.general('Running on Command Line')
            subprocess.call(commands)

    else:
        printer.general('Alignment Ignored')


    ### COLLAPSE DUPS ##############################################################################

    printer.general('COLLAPSE')

    if data.collapse:
       
        commands = [
            'python', 'collapse.py',
            '-i', '.'.join([data.trimmed_bam, "bam"]),
            '-o', data.collapsed_forward, data.collapsed_reverse,
            '-mm', str(args.max_mismatch),
            '-as', args.adapter_sequence
            ]

        printer.command(commands)

        if args.sge:
            printer.general('Running on SGE')
            sge.change(job_name='PROCESS-FASTQ:COLLAPSE')
            sge.run(commands)

        else:
            printer.general('Running on Command Line')
            subprocess.call(commands)

    else:
        printer.general('Remark Duplicates Ignored')

    ### COLLAPSE BWA ###############################################################################

    printer.general('BWA')

    if data.collapse_bam:

        printer.general('BWA')
  
        commands = [
            'python', 'bwa.py',
            '-f', data.collapsed_forward, data.collapsed_reverse,
            '-b', data.collapsed_bam,
            '-r', args.reference
            ]

        printer.command(commands)

        if args.sge:
            printer.general('Running on SGE')
            sge.change(job_name='PRODUSE:COLLAPSED-ALIGNMENT')
            sge.run(commands)

        else:
            printer.general('Running on Command Line')
            subprocess.call(commands)

    else:
        printer.general('Alignment Ignored')

    ### CLEAN ######################################################################################

    printer.general('CLEAN')

    ### END ########################################################################################

    printer.general('END')
