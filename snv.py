import argparse
import printer
import nucleotide
import position
import pysam

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
        default=0.1,
        type=float,
        help=""
        );
    args = parser.parse_args()
    bamfile = pysam.AlignmentFile(args.input, 'r');
    fastafile = pysam.FastaFile(args.reference);
    posCollection = position.PosCollectionCreate(bamfile, fastafile, filter_overlapping_reads = True);
    
    output = None;
    if not args.output == "-":
        output = open(args.output, 'w');

    for pos in posCollection:

        if pos.is_variant( \
            adapter_sequence = args.adapter_sequence, \
            max_mismatch = args.max_mismatch, \
            alt_base_count_threshold = args.alt_base_count_threshold, \
            strand_bias_threshold = args.strand_bias_threshold, \
            variant_allele_fraction_threshold = args.variant_allele_fraction_threshold ):

            if args.output == "-":
                print(str(pos));
            else:
                output.write(str(pos));
                output.write('\n');
    
    if not args.output == "-":
        output.close();