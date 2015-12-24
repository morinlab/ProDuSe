import argparse
import printer

def parse():
    desc = "Convert paired-end fastq files with variable adapters into sub-fastq's containing same variable adapter sequence"
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument(
        "-f", "--fastqs",
        nargs=2,
        help="The untrimmed, raw fastqs generated from paired end sequencing"
        )

    parser.add_argument(
        "-tf", "--trimmed_fastqs",
        nargs=2,
        help="The trimmed fastqs generated from the first stage of the pipeline"
        )

    parser.add_argument(
        "-tb", "--trimmed_bam",
        help="The aligned bam generated from the second stage of the pipeline"
        )

    parser.add_argument(
        "-cf", "--collapsed_fastqs",
        nargs=2,
        help="The collapsed fastqs generated from the third stage of the pipeline"
        )

    parser.add_argument(
        "-cb", "--collapsed_bam",
        help="The collapsed bam generated from the fourth stage of the pipeline"
        )

    parser.add_argument(
        "-mm", "--max_mismatch",
        type=int,
        default=3,
        help="The maximum number of mismatching nucleotides to see in the adapter and reads. Default is the length of the adapter (i.e. do not discard any reads)"
        )

    parser.add_argument(
        '-as', "--adapter_sequence",
        type=str,
        default="WSWSWGACT"
        )

    parser.add_argument(
        "-o", "--output_directory",
        default=".",
        help="Specify and output directory, otherwise the current working directory will be used"
        )

    parser.add_argument(
        "-p", "--prefix",
        default="tmp",
        help="Specify a prefix instead of tmp"
        )

    parser.add_argument(
        "-r", "--reference",
        help="Specify the reference genome to align each set of sub-fastqs to, must contain bwa indexes"
        )

    parser.add_argument(
        "-s", "--sge",
        action="store_true",
        default=False,
        help="If cluster services are available, jobs will be submitted using qsub"
        )

    parser.add_argument(
        "-sd", "--source_directory",
        default="",
        help="A directory to Source before each qsub submition"
        )

    parser.add_argument(
        "-1", "--trim_fastq",
        action="store_true",
        default=False,
        help="Only trim Barcodes"
        )

    parser.add_argument(
        "-2", "--trim_bam",
        action="store_true",
        default=False,
        help="Only run BWA"
        )

    parser.add_argument(
        "-3", "--collapse_fastq",
        action="store_true",
        default="True",
        help="Only collapse Fastqs"
        )

    parser.add_argument(
        "-4", "--collapse_bam",
        action="store_true",
        default="True",
        help="Only run BWA on collapsed fastqs"
        )

    parser.add_argument(
        "-5", "--call_snv",
        action="store_true",
        default="True",
        help="Only call variants"
        )

    parser.add_argument(
        "-6", "--call_sv",
        action="store_true",
        default="True",
        help="Only call variants"
        )

    parser.add_argument(
        "-7", "--call_cnv",
        action="store_true",
        default="True",
        help="Only call variants"
        )

    parser.add_argument(
        "-8", "--clean",
        action = "store_true",
        default="False",
        help="Only Clean tmp Directory"
        )

    parser.add_argument(
        "--batch",
        type=str,
        help="Submit Multiple Jobs At Once"
        )

    args = parser.parse_args()  

    return args

def check( args ):
    if not sum([not args.fastqs == None, not args.trimmed_fastqs == None, not args.trimmed_bam == None, not args.collapsed_fastqs == None, not args.batch == None]) == 1:
        printer.issue('You must specify exactly one of the fastqs, trimmed_fastqs, trimmed_bam, collapsed_fastqs or batch options')


        





