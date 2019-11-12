from collections import defaultdict
import gzip
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import argparse
import box_viz
import pandas as pd



def linear_search(key, L):
    for i  in range(len(L)):
        curr =  L[i]
        if key == curr:
            return i
    return -1


def binary_search(key, D):
    lo = -1
    hi = len(D)
    while (hi - lo > 1):
        mid = (hi + lo) // 2

        if key == D[mid][0]:
            return D[mid][1]

        if (key < D[mid][0]):
            hi = mid
        else:
            lo = mid

    return -1



def main():
    parser = argparse.ArgumentParser(description=
                                     'boxplot gene expression distribution across \
                                      tissue groups or tissue types for a gene')

    parser.add_argument('--fnt', '-file_name_gene_tmp',
                        type=str,
                        help='RNA Seq gene tmp file (gziped)',
                        required = False)

    parser.add_argument('--fnb', '-file_name_annotations',
                        type=str,
                        help='Sample attributes (tissue types) file (txt)',
                        required=True)

    parser.add_argument('--fna', '-file_name_ages',
                        type=str,
                        help='Sample phenotypes (ages) file (txt)',
                        required=True)

    parser.add_argument('--t', '-tissue',
                        type=str,
                        help='Input either tissue group (SMTS) or tissue type (SMTSD)',
                        required=False)

    parser.add_argument('--tg', '-target_gene',
                        type=str,
                        help='Input target gene',
                        required =False)

    parser.add_argument('--bfn', '-boxplot_file_name',
                        type=str,
                        help='Name of output boxplot file',
                        required =False)

    args = parser.parse_args()


#    tpm_file_name = args.fnt
    brain_tissue_file_name = args.fnb
    age_file_name = args.fna
#    group_col_name = args.t
#    gene_name = args.tg
#    boxplot_name = args.bfn

#    sample_id_col_name = 'SAMPID'



    #Create brain tissue to id dict. result looks like....
    # = {'Cortex':[id1, id2, id3], 'Medulla':[id4, id5, id6]}
    brain_tissue_to_id_dict = {}
    for line in open(brain_tissue_file_name):
        strip = line.rstrip().split('\t')
        tissue = strip[1]
        Id = strip[0]
        if tissue not in brain_tissue_to_id_dict:
            brain_tissue_to_id_dict[tissue] = []
        brain_tissue_to_id_dict[tissue].append(Id)
    #print(brain_tissue_to_id_dict)

    #Create age to id dict. result looks like....
    # = {'20-29':[id1, id2, id3], '30-39':[id4, id5, id6]}
    age_to_id_dict = {}
    with open(age_file_name) as afn:
        next(afn)
        for line in afn:
            strip = line.rstrip().split('\t')
            age = strip[2]
            Id = strip[0]
            if age not in age_to_id_dict:
                age_to_id_dict[age] = []
            age_to_id_dict[age].append(Id)
        #print(age_to_id_dict)



######
#binary_serach
#    samples = []
#    for line in open(brain_tissue_file_name):
#        samples.append(line.rstrip().split('\t'))

#    print(samples)
#    group_col_idx = linear_search(group_col_name, sample_info_header)
#    sample_id_col_idx = linear_search(sample_id_col_name, sample_info_header)

#    groups = []
#    members = []

#    for row_idx in range(len(samples)):
#        sample = samples[row_idx]
#        sample_name = sample[sample_id_col_idx]
#        curr_group = sample[group_col_idx]

#        curr_group_idx = linear_search(curr_group, groups)

#        if curr_group_idx == -1:
#            curr_group_idx = len(groups)
#            groups.append(curr_group)
#            members.append([])

#        members[curr_group_idx].append(sample_name)

#    version = None
#    dim = None
#    data_header = None

#    gene_name_col = 1


#    group_counts = [ [] for i in range(len(groups)) ]

#    for l in gzip.open(tpm_file_name, 'rt'):
#        if version == None:
#            version = l
#            continue

#        if dim == None:
#            dim = [int(x) for x in l.rstrip().split()]
#            continue

#        if data_header == None:
#            data_header = []
#            i = 0
#            for field in l.rstrip().split('\t'):
#                data_header.append([field, i])
#                i += 1
#            data_header.sort(key=lambda tup: tup[0])

#            continue

#        A = l.rstrip().split('\t')

#        if A[gene_name_col] == gene_name:
#            for group_idx in range(len(groups)):
#                for member in members[group_idx]:
#                    member_idx = binary_search(member, data_header)
#                    if member_idx != -1:
#                        group_counts[group_idx].append(int(A[member_idx]))
#            break


    #box_viz.boxplot(group_counts, groups, group_col_name, gene_name, boxplot_name)

if __name__ == '__main__':
    main()
