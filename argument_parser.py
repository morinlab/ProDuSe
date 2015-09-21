import argparse
import printer

def parse():
    desc = "Convert paired-end fastq files with variable adapters into sub-fastq's containing same variable adapter sequence"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-pe", "--paired_end",
        nargs=2,
        help="Takes the forward and reverse fastq files from paired end sequencing"
        )
    parser.add_argument(
        "-tpe", "--trimmed_paired_end",
        nargs=2,
        help="Takes the forward and reverse fastq files from paired end sequencing whose adapters have already been trimmed"
        )
    parser.add_argument(
        "-b", "--bam",
        nargs=1,
        help="Takes the sorted, aligned bam file"
        )
    parser.add_argument(
        "-mb", "--marked_bam",
        nargs=1,
        help="Takes the marked, sorted, aligned bam file"
        )
    parser.add_argument(
        "-mm", "--max_mismatch",
        type=int,
        help="The maximum number of mismatching nucleotides to see in the adapter and reads. Default is the length of the adapter (i.e. do not discard any reads)"
        )
    parser.add_argument(
        "-fp", "--file_prefix",
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
        help="Specify the reference genome to align each set of sub-fastqs to"
        )
    parser.add_argument(
        "-pj", "--picard_jar",
        help="Specify the picard jar file to run MarkDuplicates"
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
        "-F", "--trim_fastq",
        action="store_true",
        default=False,
        help="Only trim Barcodes"
        )
    parser.add_argument(
        "-B", "--bwa",
        action="store_true",
        default=False,
        help="Only run BWA"
        )
    parser.add_argument(
        "-P", "--picard",
        action="store_true",
        default=False,
        help="Only run Picard"
        )
    parser.add_argument(
       "-R", "--remark_bam",
       action="store_true",
       default="True",
       help="Only remark Bam"
       )

    args = parser.parse_args()  

    return args

def check( args ):
    if sum([args.paired_end != None, args.trimmed_paired_end != None, args.bam != None, args.marked_bam != None]) != 1:
        printer.issue('You must specify exactly one of the paired_end, trimmed_paired_end, bam or marked_bam options')





