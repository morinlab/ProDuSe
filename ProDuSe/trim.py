import configargparse
import fastq
import nucleotide
import re
import time
import sys
import os

desc = "Trim paired-end fastq files that contain an adapter sequence (paste this sequence in QNAME)"
parser = configargparse.ArgParser( description=desc )
parser.add(
    "-c", "--config",
    required=False,
    is_config_file=True,
    help="An optional configuration file for any of the input arguments."
    )
parser.add(
    "-i", "--input",
    required=True,
    action="append",
    type=str,
    help="A pair of fastq files for reading which contain flanked randomized adapter sequences"
    )
parser.add(
    "-o", "--output",
    required=True,
    action="append",
    type=str,
    help="A pair of empty fastq files for writing"
    )
parser.add(
    "-as", "--adapter_sequence",
    type=str,
    required=True,
    help="The randomized adapter sequence flanked in input fastq files described using IUPAC bases"
    )
parser.add(
    "-ap", "--adapter_position",
    type=str,
    required=True,
    help="The positions in the adapter sequence to include in distance calculations, 0 for no, 1 for yes"
    )
parser.add(
    "-mm", "--max_mismatch",
    type=int,
    required=True,
    help="The maximum number of mismatches allowed between the expected and actual adapter sequences",
    )
parser.add(
    "-v",
    action="store_true",
    help="Instead, output entries that are distant from the adapter sequence"
    )
parser.add(
    "-t",
    action="store_true",
    help="Instead, output entries without trimming the adapter sequence"
    )
if __name__ == '__main__':

    args = parser.parse_args()

    if not len(args.input) == 2:
        parser.error('--input must be a pair (i.e. a sized two list) of fastq files')

    if not len(args.output) == 2:
        parser.error('--output must be a pair (i.e. a sized two list) of fastq files')

    print_prefix = "PRODUSE-TRIM       " ;
    sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Starting...\n")

    if not len(args.adapter_position) == len(args.adapter_sequence):
        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Error: adapter_position and adapter_sequence must have same length\n")
        sys.exit(1)

    if not os.path.isfile(args.input[0]):
        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Error: input fastq files does not exist - " + args.input[0] + "\n")
        sys.exit(1)

    if not os.path.isfile(args.input[1]):
        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Error: input fastq files does not exist - " + args.input[1] + "\n")
        sys.exit(1)       

    if os.path.isfile(args.output[0]) or os.path.isfile(args.output[1]) :
        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Warning: output fastq file(s) already exist, appending to file(s)\n")

    # Determine Possible Adapters
    ref_adapter = ''.join([args.adapter_sequence, args.adapter_sequence])
    ref_indexes = list(''.join([args.adapter_position, args.adapter_position]))
    ref_indexes = [ i for i in range(len(ref_indexes)) if ref_indexes[i] == "1" ]

    # Check If Input is Gzip and call appropariate FastqOpen
    read = 'r'
    is_input_one_gzipped = not not re.search('.*\.gz', args.input[0])
    is_input_two_gzipped = not not re.search('.*\.gz', args.input[1])
    if is_input_one_gzipped and is_input_two_gzipped:
        read = ''.join([read, 'g']);
    elif is_input_one_gzipped or is_input_two_gzipped:
        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Error: input files must both be gzipped or both uncompressed\n")
        sys.exit(1)

    # Check If Output is Gzip and call appropariate FastqOpen
    write = 'w'
    is_output_one_gzipped = not not re.search('.*\.gz', args.output[0])
    is_output_two_gzipped = not not re.search('.*\.gz', args.output[1])
    if is_output_one_gzipped and is_output_two_gzipped:
        write = ''.join([write, 'g']);
    elif is_output_one_gzipped or is_output_two_gzipped:
        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Error: output files must both be gzipped or both uncompressed\n")
        sys.exit(1)

    # Open Fastq files for reading and writing
    forward_input = fastq.FastqOpen(args.input[0], read)
    reverse_input = fastq.FastqOpen(args.input[1], read)
    forward_output = fastq.FastqOpen(args.output[0], write)
    reverse_output = fastq.FastqOpen(args.output[1], write)

    # For each read in the forward and reverse fastq files, trim the adapter, at it to the read id
    # and output this new fastq record to the temprary output fastq
    count , discard = 0 , 0
    #sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Discard Rate:" + str(discard / count) + "\tCount:" + str(count) + "\n")
    # Loop over every fastq entry in file
    for forward_read in forward_input:

        count += 1

        if count % 100000 == 0:
            sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Discard Rate:" + str(round(float(discard) / float(count), 3) * 100) + "%    Count:" + str(count) + "\n")

        # Fetch associated reverse read
        reverse_read = reverse_input.next()

        # Get the actual adapter sequence from both fastq files
        # As well as trimming the seq and qual bases from these regions
        forward_adapter = forward_read.trim(len(args.adapter_sequence))
        reverse_adapter = reverse_read.trim(len(args.adapter_sequence))

        # Glue these together
        sequenced_adapter = ''.join([forward_adapter.seq, reverse_adapter.seq]);
        
        # Determine the distance between the actual and expected adapter sequence
        if args.v and nucleotide.distance(sequenced_adapter, ref_adapter, ref_indexes) <= args.max_mismatch:
            
            # Do not include in output fastq files if distance is greater then max_mismatch
            discard += 1
            continue

        elif not args.v and nucleotide.distance(sequenced_adapter,ref_adapter, ref_indexes) > args.max_mismatch:

            discard += 1
            continue

        if args.t:

            # Glue back the trimmed sequences and quality scores            
            forward_read.prepend(forward_adapter)
            reverse_read.prepend(reverse_adapter)
         
        else:

            # Update read names to include adapter sequence prefixes
            forward_read.id = ''.join(['@', forward_adapter.seq, reverse_adapter.seq, ':', forward_read.id[1:]])
            reverse_read.id = ''.join(['@', forward_adapter.seq, reverse_adapter.seq, ':', reverse_read.id[1:]])

        # Write entries to fastq outputs
        forward_output.next(forward_read)
        reverse_output.next(reverse_read)

    if count % 100000 != 0:

        sys.stdout.write(print_prefix + time.strftime('%X') + "    " + "Discard Rate:" + str(round(float(discard) / float(count), 3) * 100) + "%    Count:" + str(count) + "\n")
    
    # Close Fastq Files
    forward_input.close()
    reverse_input.close()
    forward_output.close()
    reverse_output.close()
