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
        default="WSWSWGACT",
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
    args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.input, 'r');
    fastafile = pysam.FastaFile(args.reference);
    posCollection = position.PosCollectionCreate(bamfile, fastafile);
    output = open(args.output, 'w');

    for pos in posCollection:
        if pos.is_variant(adapter_sequence=args.adapter_sequence, max_mismatch=args.max_mismatch):
            output.write(str(pos));
            output.write('\n');

    output.close();
