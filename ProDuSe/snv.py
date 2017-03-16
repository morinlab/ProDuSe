#! /usr/bin/env python

# If not installed, or running in python2, this works fine
try:
    import position
    import bed
except ImportError:
    # If installed and running in python3
    from ProDuSe import position
    from ProDuSe import bed

import configargparse
import pysam
import os

"""
Processes command line arguments

Returns:
    args: A namespace object listing command line paramters
"""
desc = "Call SNVs on collapsed BAMs containing adapter sequence information"
parser = configargparse.ArgParser(description=desc)
parser.add(
    "-c", "--config",
    required=False,
    is_config_file=True,
    help="An configuration file which can list one or more input arguments"
    )
parser.add(
    "-i", "--input",
    required=True,
    type=lambda x: is_valid_file(x, parser),
    help="An input bam file, generated using \'collapse.py\' and \'bwa.py\'"
    )
parser.add(
    "-o", "--output",
    required=True,
    type=str,
    help="Output Variant Call Format (vcf) file"
    )
parser.add(
    "-r", "--reference",
    required=True,
    type=str,
    help="Reference genome, in FASTA format"
    )
parser.add_argument(
    "-tb", "--target_bed",
    required=False,
    help="A tab-delinated file listing regions on which variant calling will be restricted to"
    )
parser.add_argument(
    "-vaft", "--variant_allele_fraction_threshold",
    default=0.01,
    type=float,
    help="Minimum variant frequency threshold, above which variants will be called as being real for each strant [Default: %(default)s]"
    )
parser.add_argument(
    "-mo", "--min_molecules",
    default=400,
    type=int,
    help="Number of unique supporting molecules required to call a variant as real for each strand. Reduce this if you are running only on positions you expect to be mutated [Default: %(default)s]"
    )
parser.add_argument(
    "-mum", "--mutant_molecules",
    default=3,
    required=False,
    type=int,
    help="Number of TOTAL molecules (i.e. on the forward and reverse strand) required to call a variant as real (set to 0 if you are forcing variant calling at known sites)"
    )
parser.add_argument(

    "-mrpu", "--min_reads_per_uid",
    default=2,
    required=True,
    type=int,
    help="Bases with support between MRPU and SSBT will be classified as a weak supported base"
    )
parser.add_argument(
    "-ssbt", "--strong_supported_base_threshold",
    default=3,
    type=int,
    help="Bases with support equal to or greater then SSBT, will be classified as a strong supported base."
    )

parser.add_argument(
    "-eds", "--enforce_dual_strand",
    action='store_true',
    help="require at least one molecule to be read in each of the two directions"
    )
parser.add_argument(
    "-mq", "--min_qual",
    default=3,
    type=int,
    help="Minimum base quality threshold, below which a base will be treated as \'N\'")
# For backwards compatability only. These arguments do nothing
parser.add_argument(
    "-abct", "--alt_base_count_threshold",
    default=5,
    type=int,
    help=configargparse.SUPPRESS  #  "The minimum number of alternative bases to identify separately in the positive and negative read collections"
    )
parser.add_argument(
    "-sbt", "--strand_bias_threshold",
    help=configargparse.SUPPRESS
    )
parser.add_argument(
    "--adapter_sequence",
    help=configargparse.SUPPRESS  #  "The minimum number of alternative bases to identify separately in the positive and negative read collections"
    )
parser.add_argument(
    "--adapter_position",
    help=configargparse.SUPPRESS
    )
parser.add_argument(
    "--max_mismatch",
    help=configargparse.SUPPRESS
    )

def is_valid_file(file, parser):
    """
    Checks to ensure the specified file exists, and throws an error if it does not

    Args:
        file: A filepath
        parser: An argparse.ArgumentParser() object. Used to throw an exception if the file does not exist

    Returns:
        type: The file itself

    Raises:
        parser.error(): An ArgumentParser.error() object, thrown if the file does not exist
    """
    if os.path.isfile(file):
        return file
    else:
        parser.error("The file %s does not exist" % (file))


def main(args=None):

    if args is None:
        args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.input, 'rb')
    fastafile = pysam.FastaFile(args.reference)
    targetbed = None
    if not args.target_bed == None:
        targetbed = bed.BedOpen(args.target_bed, 'r');
    #else:
    #    sys.exit("You need to specify a targetbed")

    posCollection = position.PosCollectionCreate(bamfile, fastafile, filter_overlapping_reads=True, target_bed=targetbed, min_reads_per_uid=args.min_reads_per_uid, min_base_qual=args.min_qual);

    output = None;
    if not args.output == "-":
        output = open(args.output, 'w');

    # if args.mode == 'validation':
    #     for pos in posCollection:
    #         if pos.coords2() in targetbed:
    #             is_variant = pos.is_variant( \
    #                 adapter_sequence = args.adapter_sequence, \
    #                 max_mismatch = args.max_mismatch, \
    #                 alt_base_count_threshold = args.alt_base_count_threshold, \
    #                 strand_bias_threshold = args.strand_bias_threshold, \
    #                 variant_allele_fraction_threshold = args.variant_allele_fraction_threshold );

    ##             if args.output == "-":
    #                 print(str(pos));
    #             else:
    #                 output.write(str(pos));
    #                 output.write('\n');

    # elif args.mode == 'discovery':

    m = 1
    first = True
    for pos in posCollection:

        pos.calc_base_stats(
            min_reads_per_uid=args.strong_supported_base_threshold
            )
        # print "done %s  positions in for pos in posCollection at %s" % (m,pos.coords())

        m += 1

        if pos.is_variant(args.variant_allele_fraction_threshold, args.min_molecules, args.enforce_dual_strand, args.mutant_molecules):

            if first:
                pos.write_header(output)
                first = False
            pos.write_variant(output)
            # output.write(pos.coords() + "\n")
            # output.write(pos.ref + " > " + ''.join(pos.alt) + "\n")
            # output.write(str(pos))

        # if pos.coords2() in targetbed:

        #     if args.output == "-":
        #         print(str(pos));

        #     else:
        #         output.write(str(pos));
        #         output.write('\n');

    # else:

    #     for pos in posCollection:

    #         is_variant = pos.is_variant( \
    #             adapter_sequence = args.adapter_sequence, \
    #             max_mismatch = args.max_mismatch, \
    #             alt_base_count_threshold = args.alt_base_count_threshold, \
    #             strand_bias_threshold = args.strand_bias_threshold, \
    #             variant_allele_fraction_threshold = args.variant_allele_fraction_threshold );

    #         if is_variant or pos.coords() in bedfile:

    #             if args.output == "-":
    #                 print(str(pos));
    #             else:
    #                 output.write(str(pos));
    #                 output.write('\n');

    if not args.output == "-":
        output.close()


if __name__ == '__main__':

    main()
