import argparse
import ConfigParser
import random
import numpy
import nucleotide
import pysam
import fastq
import re

class dsMolecule:

    def __init__(self, seq1, seq2):
        self.forward = seq1
        self.reverse = seq2

    def denature(self):
        return [ ssMolecule(self.forward, "+"), ssMolecule(nucleotide.reverse(self.reverse), "+")]

class ssMolecule:

    def __init__(self, seq, sense):
        self.seq = seq
        self.sense = sense

    def polymerase(self, error_rate):
        if self.sense == "+":
            return dsMolecule(self.seq, nucleotide.random_mismatch(nucleotide.complement(self.seq), error_rate))
        else:
            return dsMolecule(nucleotide.random_mismatch(nucleotide.complement(self.seq), error_rate), self.seq)

if __name__ == '__main__':

    desc = "Create test paired_end Duplex Sequenced Fastq Files"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-as", "--adapter_sequence",
        #required=True,
        help="Specify the contiguous arrangement of nucleotides in the adapter sequence using any symbol defined by the IUPAC"
        )

    parser.add_argument(
        "-c", "--config",
        required=True,
        help="Algoritm Parameter Config File"        
        )

    parser.add_argument(
        "-o", "--output",
        nargs=2,
        #required=True,
        help="Output Fastq Files"
        )

    parser.add_argument(
        "-r", "--reference",
        #required=True,
        help="Reference To Create sequences from"
        )

    args = parser.parse_args()

    config = ConfigParser.RawConfigParser()
    config.read(args.config)
    pcr_cycle_count = config.getint("pcr parameters", "pcr_cycle_count")
    pcr_cycle_error_rate = config.getfloat("pcr parameters", "pcr_cycle_error_rate")
    seq_read_length = config.getint("sequencing parameters", "seq_read_length")
    seq_error_rate = config.getfloat("sequencing parameters", "seq_error_rate")
    seq_molecule_count = config.getint("sequencing parameters", "seq_molecule_count")
    lib_template_length_avg = config.getint("library parameters", "lib_template_length_avg")
    lib_template_length_sd = config.getint("library parameters", "lib_template_length_sd")
    lib_molecule_count = config.getint("library parameters", "lib_molecule_count")

    write = 'w'
    is_output_one_gzipped = not not re.search('.*\.gz', args.output[0])
    is_output_two_gzipped = not not re.search('.*\.gz', args.output[1])
    if is_output_one_gzipped and is_output_two_gzipped:
        write = ''.join([write, 'g']);
    elif is_output_one_gzipped or is_output_two_gzipped:
        print 'Output files must both be gzipped or both uncompressed\n'
        sys.exit()

    forward_output = fastq.FastqOpen(args.output[0], write)
    reverse_output = fastq.FastqOpen(args.output[1], write)

    print "Extracting DNA from Cells..."

    fastafile = pysam.FastaFile(args.reference)
    seq = fastafile.fetch("chr1", 1, 1000);
    indexes_start = list(numpy.random.choice( range(len(seq))[ 0 : len(seq) - lib_template_length_avg ] , lib_molecule_count, replace=True))
    template_lengths = [ int(a) for a in list(numpy.random.normal(lib_template_length_avg, lib_template_length_sd, lib_molecule_count))]
    indexes_end = [None] * len(indexes_start)
    for i in range(len(indexes_start)):
        indexes_end[i] = indexes_start[i] + template_lengths[i]
    for i in range(len(indexes_end)):
        if indexes_end[i] >= len(seq):
            indexes_end[i] = len(seq) - 1

    print "Ligating Randomized Duplex Tags..."

    list_of_dna = [None] * lib_molecule_count
    for i in range(len(indexes_start)):
        cur_seq = seq[indexes_start[i]:indexes_end[i]]
        cur_seq_complement = nucleotide.complement(cur_seq)
        alpha = nucleotide.random_unambiguous(args.adapter_sequence)
        alpha_complement = nucleotide.complement(alpha)
        beta_reverse = nucleotide.reverse(nucleotide.random_unambiguous(args.adapter_sequence))
        beta_reverse_complement = nucleotide.complement(beta_reverse)
        list_of_dna[i] = dsMolecule(alpha + cur_seq + beta_reverse_complement, alpha_complement + cur_seq_complement + beta_reverse)

    print "Running Polymerase Chain Reaction..."

    for i in range(pcr_cycle_count):
        list_of_rna = [ item for dna in list_of_dna for item in dna.denature() ]
        list_of_dna = [ rna.polymerase(pcr_cycle_error_rate) for rna in list_of_rna ] 

    print "Sequencing..."

    list_of_dna = numpy.random.choice(list_of_dna, seq_molecule_count, replace=False)

    for i in range(len(list_of_dna)):
        read_id = "".join(["@R",str(i)])
        forward = nucleotide.random_mismatch(list_of_dna[i].forward[0:seq_read_length], seq_error_rate)
        reverse = nucleotide.random_mismatch(nucleotide.reverse(list_of_dna[i].reverse[len(list_of_dna[i].reverse)-seq_read_length:len(list_of_dna[i].reverse)]), seq_error_rate)
        forward_output.next(fastq.Fastq(read_id, forward, '+', ''.join(numpy.random.choice(['A','B','C','D','E','F','G','H'], seq_read_length, replace=True))))
        reverse_output.next(fastq.Fastq(read_id, reverse, '+', ''.join(numpy.random.choice(['A','B','C','D','E','F','G','H'], seq_read_length, replace=True))))