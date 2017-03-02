import argparse
import printer
import fastq
import nucleotide
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == '__main__':

    desc = "Quality Control on the adapter sequences in a pair of fastq files"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-i", "--input",
        required=True,
        nargs='+',
        help="Takes many collapsed fastq file (one per pair)"
        )
    parser.add_argument(
        "-as", "--adapter_sequence",
        type=str,
        default="NNNWSMRWSYWKMWWT",
        help="The adapter sequence following IUPAC DNA naming"
        )
    parser.add_argument(
        "-ap", "--adapter_position",
        type=str,
        default="---++++++++++++-",
        help="The positions in the adapter sequence to include in distance calcualtions"
        )
    parser.add_argument(
        "-op", "--output_with_prefix",
        type=str,
        help="The output directory with a Prefix"
        )

    args = parser.parse_args()
    ref_adapter = ''.join([args.adapter_sequence, args.adapter_sequence])
    ref_indexes = list(''.join([args.adapter_position, args.adapter_position]))
    ref_indexes = [ i for i in range(len(ref_indexes)) if ref_indexes[i] == "+" ]

    # Check If Input Files are Gzipped and call appropariate FastqOpen
    read = ['r'] * len(args.input)
    for i in range(len(args.input)):
        is_input_one_gzipped = not not re.search('.*\.gz', args.input[i])
        if is_input_one_gzipped:
            read[i] = ''.join([read[i], 'g']);

    ### Obtain Metric Data For Plotting

    items = range(len(args.input))
    positive = [ {} for i in items]
    negative = [ {} for i in items]
    positive_set = [ {} for i in items]
    negative_set = [ {} for i in items]
    positive_count = [ {} for i in items]
    negative_count = [ {} for i in items]
    positive_adapter = [ {} for i in items]
    negative_adapter = [ {} for i in items]
    uid_counts = [[] for i in items]
    
    for i in items:

        forward_input = fastq.FastqOpen(args.input[i], read[i])
        
        for element in forward_input:
        
            tmp = element.id.split(":")
            item = ':'.join(tmp[1:5])
            uid_counts[i].append(int(tmp[5]))

            if tmp[6] == "+" :
                if not item in positive[i]:
                    positive[i][item] = 0
                    positive_count[i][item] = int(tmp[5])
                    positive_set[i][item] = []
                    positive_adapter[i][item] = []
                positive[i][item] += 1
                positive_count[i][item] += int(tmp[5])           
                positive_set[i][item].append(tmp[7:])
                positive_adapter[i][item].append(tmp[0][1:])

            else :
                if not item in negative[i]:
                    negative[i][item] = 0 
                    negative_count[i][item] = int(tmp[5])
                    negative_set[i][item] = []
                    negative_adapter[i][item] = []
                negative[i][item] += 1
                negative_count[i][item] += int(tmp[5])
                negative_set[i][item].append(tmp[7:])
                negative_adapter[i][item].append(tmp[0][1:])


    distribution = [[ [0] * 32 for i in range(11) ] for j in items]
    family_distribution = [ [ [0] * 32 for i in range(11) ] for j in items]
    for i in items:
        for j in negative[i]:
            if negative[i][j] >= 3:
                print(negative[i][j])
                print(negative_set[i][j])
                distribution[i][negative[i][j]][min(nucleotide.dist_pairwise(ref_indexes, *negative_set[i][j]))] += 1
                family_distribution[i][negative[i][j]][min(nucleotide.dist_list(negative_adapter[i][j], ref_indexes))] += 1

    data_frame_one = []
    data_frame_two = []
    for item in items:
        for family in range(2,6):
            for distance in range(len(distribution[item][family])):
                counts = distribution[item][family][distance]
                counts_two = family_distribution[item][family][distance]
                data_frame_one.append({"File" : item, "Family" : family, "Distance" : distance, "Counts" : counts })
                data_frame_two.append({"File" : item, "Family" : family, "Distance" : distance, "Counts" : counts_two})

    import pandas as pd 

    data_for_plot = pd.DataFrame(data_frame_one)
    data_for_plot_two = pd.DataFrame(data_frame_two)
    sns.tsplot(time="Distance", condition="Family",value="Counts", unit="File", data=pd.DataFrame(data_frame_one))
    sns.plt.savefig(".".join([args.output_with_prefix, "minFamilyDistance.png"]))
    sns.plt.close()
    sns.tsplot(time="Distance", condition="Family",value="Counts", unit="File", data=pd.DataFrame(data_frame_two))
    sns.plt.savefig(".".join([args.output_with_prefix, "minFamilyConsensusDistance.png"]))
    sns.plt.close()

    
    cur_min = float("inf")
    cur_max = float("-inf")
    for i in items:
        cur_min = min([cur_min, min(uid_counts[i])])
        cur_max = max([cur_max, max(uid_counts[i])])

    count_data = [ [0] * (cur_max+1) for i in items ]
    for i in items:
        for j in uid_counts[i]:
            count_data[i][j] += 1

    data_frame_tmp = []
    for i in items:
        for j in range(cur_min, cur_max+1):
            data_frame_tmp.append({"File" : i, "Family" : j, "Counts" : count_data[i][j]})


    g = sns.factorplot(data=pd.DataFrame(data_frame_tmp), x="Family", y="Counts", units="File", kind="bar")
    sns.plt.savefig(".".join([args.output_with_prefix, "uidDistribution.png"]))
    sns.plt.close()


#     read_count = [0] * 100

#     for i in positive_count:
#         if positive_count[i] >= 100:
#             read_count[99] += 1
#         else:
#             read_count[positive_count[i]-1] += 1

#     for i in negative_count:
#         if negative_count[i] >= 100:
#             read_count[99] += 1
#         else:
#             read_count[negative_count[i]-1] += 1

#     count = 0
#     base_wise_mismatch = [0] * len(ref_adapter)
#     mismatch_frequency = [0] * len(ref_adapter)
#     base_wise_mismatch_matrix = [[0 for i in range(len(ref_adapter))] for j in range(len(ref_adapter))]
#     mismatch_frequency_matrix = [[0 for i in range(len(ref_adapter))] for j in range(len(ref_adapter))]
#     for item in positive:
#     	if item in negative:
#             if positive[item] == 1 and negative[item] == 1:
#                 count += 1
#                 for i in range((len(ref_adapter)/2)):
#                     base_wise_mismatch_matrix[i][(i*2):32] = nucleotide.base_mismatch(positive_adapter[item][(0+i):16] + positive_adapter[item][(16+i):32], negative_adapter[item][(0+i):16] + negative_adapter[item][(16+i):32], base_wise_mismatch_matrix[i][(i*2):32])
#                     mismatch_frequency_matrix[i][nucleotide.distance(positive_adapter[item][(0+i):16] + positive_adapter[item][(16+i):32], negative_adapter[item][(0+i):16] + negative_adapter[item][(16+i):32])] += 1

#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     N = len(ref_adapter)
#     ind = np.arange(N)
#     y = [ float(item) / float(count) for item in base_wise_mismatch_matrix[0] ]
#     width = 0.35
#     ax.bar(ind, y, width, color='blue')
#     xTickMarks = list(ref_adapter)
#     ax.set_xticks(ind+width)
#     xtickNames = ax.set_xticklabels(xTickMarks)
#     plt.setp(xtickNames, fontsize=10)
#     plt.savefig('base_wise_mismatch')

#     for j in range(14):
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#         #ax.set_yscale('log')
#         N = len(ref_adapter)
#         ind = np.arange(N)
#         y = [ float(item) / float(count) for item in mismatch_frequency_matrix[j] ]
#         width = 0.35
#         ax.bar(ind, y, width, color='blue')
#         xTickMarks = [str(i) for i in range(N)]
#         ax.set_xticks(ind+width)
#         xtickNames = ax.set_xticklabels(xTickMarks)
#         plt.setp(xtickNames, fontsize=10)
#         plt.savefig(''.join(["mismatch_frequency", str(j)]))

# data = [ min(nucleotide.dist_list(nucleotide.random_unambiguous("NN",2))) for i in range(100000) ]
# import itertools
# counts = [0] * 24
# for i in data:
#     counts[i] += 1

# N = len(counts)
# x = range(N)
# width = 1/1.5
# import matplotlib.pyplot as plt
# plt.bar(x,counts,width)
# plt.show()
