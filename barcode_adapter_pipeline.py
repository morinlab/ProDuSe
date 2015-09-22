#/bin/usr/env python
import qsub
import argument_parser
import printer
import setup
import subprocess

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
            'python', 'trim_fastq.py',
            '-i', data.forward, data.reverse,
            '-o', data.trimmed_forward, data.trimmed_reverse,
            '-mm', args.max_mismatch,
            '-as', args.adapter_sequence
            ]
    
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

    if data.bwa:
        
        commands = [
            'bwa', 'mem',
            '-M',
            '-t', str(args.number_of_threads),
            args.reference,
            data.trimmed_forward,
            data.trimmed_reverse,
            '>', data.bam
            ]
    
        if args.sge:
            sge.change(job_name='PROCESS-FASTQ:BWA')
            sge.run(commands)
        else:
            subprocess.call(commands)
    
    else:
        printer.general('Alignment Ignored')

    ### PICARD MARKDUPS ############################################################################

    printer.general('PICARD')

    if data.picard:
        
        commands = [
            'java', '-Xmx4g',
            '-jar', args.picard_jar, "MarkDuplicates",
            ''.join(["INPUT=", data.bam]),
            ''.join(["OUTPUT=", data.marked_bam]),
            ''.join(["METRIC_FILE=", data.metric_file]),
            'REMOVE_DUPLICATES=True',
            'ASSUME_SORTED=False',
            'SORT_ORDER=coordinate'
            ]
    
        if args.sge:
            sge.change(job_name='PROCESS-FASTQ:PICARD')
            sge.run(commands)
        else:
            subprocess.call(commands)

    else:
        printer.general('Mark Duplicates Ignored')


    ### REMARK DUPS ################################################################################

    printer.general('REMARK')

    if data.remark:
        
        commands = [
            'python', 'remark_adapter_dups.py',
            '-i', data.marked_bam,
            '-o', data.remarked.bam,
            '-mm', args.max_mismatch,
            '-as', args.adapter_sequence
            ]

        if args.sge:
            printer.general('Running on SGE')
            sge.change()
            sge.run(commands)
        else:
            printer.general('Running on Command Line')
            subprocess.call(commands)

    else:
        printer.general('Remark Duplicates Ignored')

    ### CLEAN ######################################################################################

    printer.general('CLEAN')

    ### END ########################################################################################

    printer.general('END')
