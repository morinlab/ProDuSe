import argparse
import printer
import nucleotide
import position
import pysam
import bed

if __name__ == '__main__':

    desc = "Call SNVs on collapsed BAMs containing adapter sequence information"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="File location collapsed BAM"
        )
    parser.add_argument(
        "-as", "--adapter_sequence",
        type=str,
        required=True,
        help="The adapter sequence following IUPAC DNA naming"
        )
    parser.add_argument(
        "-mm", "--max_mismatch",
        type=int,
        default=3,
        help="The maximum number of mismatches to accept in reads when pairing positive and negative reads"
        )

    parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Reference file with samtools faidx index"
        )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output variant file"
        )

    parser.add_argument(
        "-m", "--mode",
        default='discovery',
        choices=['discovery', 'validation', 'both'],
        help="select one of the modes"
        );

    parser.add_argument(
        "-pb", "--positional_bed",
        required=False,
        help="Print out snv information for positions in the bed file in validation mode"
        )

    parser.add_argument(
        "-tb", "--target_bed",
        required=False,
        help="Restrict read fetching to intervals defined in this file"
        )

    parser.add_argument(
        "-abct", "--alt_base_count_threshold",
        default=5,
        type=int,
        help="The minimum number of alternative bases to identify separately in the positive and negative read collections"
        );

    parser.add_argument(
        "-sbt", "--strand_bias_threshold",
        default=0.2,
        type=float,
        help=""
        );

    parser.add_argument(
        "-vaft", "--variant_allele_fraction_threshold", 
        default=0.01,
        type=float,
        help=""
        );

    args = parser.parse_args()

    if (args.mode == 'validation' or args.mode == 'both') and not args.target_bed:
        sys.exit("You must specify a positional bed in validation mode.");

    if args.mode == 'discovery' and args.positional_bed:
        sys.warning("Mode is set to discovery, the positional bed file will be ignored.");

    bamfile = pysam.AlignmentFile(args.input, 'rb');
    fastafile = pysam.FastaFile(args.reference);
    targetbed = None
    if args.target_bed:
        targetbed = bed.BedOpen(args.target_bed, 'r');

    posCollection = position.PosCollectionCreate(bamfile, fastafile, filter_overlapping_reads = True, target_bed = targetbed);

    bedfile = None;
    if args.positional_bed:
        bedfile = bed.BedOpen(args.positional_bed, 'r');

    output = None;
    if not args.output == "-":
        output = open(args.output, 'w');

    if args.mode == 'validation':
        for pos in posCollection:
            if pos.coords2() in targetbed:
                is_variant = pos.is_variant( \
                    adapter_sequence = args.adapter_sequence, \
                    max_mismatch = args.max_mismatch, \
                    alt_base_count_threshold = args.alt_base_count_threshold, \
                    strand_bias_threshold = args.strand_bias_threshold, \
                    variant_allele_fraction_threshold = args.variant_allele_fraction_threshold );

                if args.output == "-":
                    print(str(pos));
                else:
                    output.write(str(pos));
                    output.write('\n');

    elif args.mode == 'discovery':

        for pos in posCollection:

            is_variant = pos.is_variant( \
                adapter_sequence = args.adapter_sequence, \
                max_mismatch = args.max_mismatch, \
                alt_base_count_threshold = args.alt_base_count_threshold, \
                strand_bias_threshold = args.strand_bias_threshold, \
                variant_allele_fraction_threshold = args.variant_allele_fraction_threshold );

            if is_variant:

                if args.output == "-":
                    print(str(pos));
                else:
                    output.write(str(pos));
                    output.write('\n');

    else:

        for pos in posCollection:

            is_variant = pos.is_variant( \
                adapter_sequence = args.adapter_sequence, \
                max_mismatch = args.max_mismatch, \
                alt_base_count_threshold = args.alt_base_count_threshold, \
                strand_bias_threshold = args.strand_bias_threshold, \
                variant_allele_fraction_threshold = args.variant_allele_fraction_threshold );

            if is_variant or pos.coords() in bedfile:

                if args.output == "-":
                    print(str(pos));
                else:
                    output.write(str(pos));
                    output.write('\n');
    
    if not args.output == "-":
        output.close();