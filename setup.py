import printer
import os
import random
import subprocess

class Setup:

    def __init__(self, args):

        self.forward = None
        self.reverse = None
        self.trimmed_forward = None
        self.trimmed_reverse = None
        self.collapsed_forward = None
        self.collapsed_reverse = None
        self.trimmed_sam = None
        self.trimmed_bam = None
        self.trimmed_sort = None
        self.collapsed_sam = None
        self.collapsed_bam = None
        self.collapsed_sort = None
        
        self.trim = False
        self.trim_bam = False
        self.collapse = False
        self.collapse_bam = False

        printer.general('Creating Temporary Output Directory')
        
        prev_tmp_dir = '/'.join([args.output_directory, args.prefix])
        tmp_dir = prev_tmp_dir
        
        while os.path.exists(tmp_dir):
            tmp_dir = ''.join([prev_tmp_dir, str(random.randint(1000000000, 9999999999))])
        os.makedirs(tmp_dir)
        printer.general(':'.join(['Temporary output directory is', tmp_dir]))

        printer.general('Creating All Output File Names')


        if not args.fastqs == None:
            self.forward = args.fastqs[0]
            self.reverse = args.fastqs[1]

        elif not args.trimmed_fastqs == None:
            self.trimmed_forward = args.trimmed_fastqs[0]
            self.trimmed_reverse = args.trimmed_fastqs[1]

        elif not args.trimmed_bam == None:
            self.trimmed_bam = args.trimmed_bam

        elif not args.collapsed_fastqs == None:
            self.collapsed_forward = args.collapsed_fastqs[0]
            self.collapsed_reverse = args.collapsed_fastqs[1]

        if args.collapse_bam:
            self.collapsed_bam = '/'.join([tmp_dir, 'collapsed'])
            self.collapse_bam = True

        if args.collapse_fastq:
            if args.collapsed_fastqs == None:
                self.collapsed_forward = '/'.join([tmp_dir, 'collapsed_forward.fastq'])
                self.collapsed_reverse = '/'.join([tmp_dir, 'collapsed_reverse.fastq'])
                self.collapse = True
                if args.trim_bam:
                    self.trimmed_bam = '/'.join([tmp_dir, 'trimmed'])
                    self.trim_bam = True
                    if args.trim_fastq:
                        self.trimmed_forward = '/'.join([tmp_dir, 'trimmed_forward.fastq'])
                        self.trimmed_reverse = '/'.join([tmp_dir, 'trimmed_reverse.fastq'])
                        self.trim = True
                    else:
                        pass
                else:
                    pass
            else:
                pass

        
