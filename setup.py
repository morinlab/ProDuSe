import printer
import os
import random

class Setup:

    def __init__(self, args):

        self.forward = None
        self.reverse = None
        self.trimmed_forward = None
        self.trimmed_reverse = None
        self.bam = None
        self.marked_bam = None
        self.remarked_bam = None
        self.trim = False
        self.bwa = False
        self.picard = False
        self.remark = False

        printer.general('Creating Temporary Output Directory')
        prev_tmp_dir = '/'.join([args.output_directory, 'tmp'])
        tmp_dir = prev_tmp_dir
        while os.path.exists(tmp_dir):
            tmp_dir = ''.join([prev_tmp_dir, str(random.randint(1000000000, 9999999999))])
        os.makedirs(tmp_dir)
        printer.general(':'.join(['Temporary output directory is', tmp_dir]))

        printer.general('Creating All Output File Names')
        self.remarked_bam = '/'.join([tmp_dir, 'remarked.bam'])
        self.remark = True
        if args.marked_bam != None:
            self.marked_bam = args.marked_bam
        else:
            self.picard = True
            self.marked_bam = '/'.join([tmp_dir, 'marked.bam'])
            if args.bam != None:
                self.bam = args.bam
            else:
                self.bwa = True
                self.bam = '/'.join([tmp_dir, 'bwa.bam'])
                if args.trimmed_paired_end != None:
                    self.trimmed_forward = args.trimmed_paired_end[0]
                    self.trimmed_reverse = args.trimmed_paired_end[1]
                else:
                    self.trim = True
                    self.trimmed_forward = '/'.join([tmp_dir, 'forward_trimmed.fastq'])
                    self.trimmed_forward = '/'.join([tmp_dir, 'reverse_trimmed.fastq'])
                    self.forward = args.paired_end[0]
                    self.reverse = args.paired_end[1]
