import argparse
import printer
import fastq
import nucleotide
import re
import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    desc = "Quality Control on the adapter sequences in a pair of fastq files"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-i", "--input",
        nargs=1,
        required=True,
        help="Takes a single collapsed fastq file"
        )
    parser.add_argument(
        "-as", "--adapter_sequence",
        type=str,
        default="NNNWSMRWSYWKMWWT",
        help="The adapter sequence following IUPAC DNA naming"
        )
    args = parser.parse_args()
    ref_adapter = ''.join([args.adapter_sequence, args.adapter_sequence])

    # Check If Input is Gzip and call appropariate FastqOpen
    read = 'r'
    is_input_one_gzipped = not not re.search('.*\.gz', args.input[0])
    if is_input_one_gzipped:
        read = ''.join([read, 'g']);

    forward_input = fastq.FastqOpen(args.input[0], read)
    positive = {}
    negative = {}
    positive_adapter = {}
    negative_adapter = {}
    for read in forward_input:
    	tmp = read.id.split(":");
    	item = ':'.join(tmp[1:5])
    	if tmp[6] == "+":
    		if item in positive:
    			positive[item] += 1
    		else:
    			positive[item] = 1
    		positive_adapter[item] = tmp[0][1:]
    	else :
     		if item in negative:
    			negative[item] += 1
    		else:
    			negative[item] = 1   		
    		negative_adapter[item] = tmp[0][1:]

    count = 0
    base_wise_mismatch = [0] * len(ref_adapter)
    mismatch_frequency = [0] * len(ref_adapter)
    base_wise_mismatch_matrix = [[0 for i in range(len(ref_adapter))] for j in range(len(ref_adapter))]
    mismatch_frequency_matrix = [[0 for i in range(len(ref_adapter))] for j in range(len(ref_adapter))]
    for item in positive:
    	if item in negative:
            if positive[item] == 1 and negative[item] == 1:
                count += 1
                for i in range((len(ref_adapter)/2)):
                    base_wise_mismatch_matrix[i][(i*2):32] = nucleotide.base_mismatch(positive_adapter[item][(0+i):16] + positive_adapter[item][(16+i):32], negative_adapter[item][(16+i):32] + negative_adapter[item][(0+i):16], base_wise_mismatch_matrix[i][(i*2):32])
                    mismatch_frequency_matrix[i][nucleotide.distance(positive_adapter[item][(0+i):16] + positive_adapter[item][(16+i):32], negative_adapter[item][(16+i):32] + negative_adapter[item][(0+i):16])] += 1

    fig = plt.figure()
    ax = fig.add_subplot(111)
    N = len(ref_adapter)
    ind = np.arange(N)
    y = [ float(item) / float(count) for item in base_wise_mismatch ]
    width = 0.35
    ax.bar(ind, y, width, color='blue')
    xTickMarks = list(ref_adapter)
    ax.set_xticks(ind+width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, fontsize=10)
    plt.savefig('base_wise_mismatch')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.set_yscale('log')
    N = len(ref_adapter)
    ind = np.arange(N)
    y = [ float(item) / float(count) for item in mismatch_frequency ]
    width = 0.35
    ax.bar(ind, y, width, color='blue')
    xTickMarks = [str(i) for i in range(N)]
    ax.set_xticks(ind+width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, fontsize=10)
    plt.savefig('mismatch_frequency')