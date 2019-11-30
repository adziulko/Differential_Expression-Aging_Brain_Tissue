from collections import defaultdict
import gzip
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import argparse
import box_viz
import pandas as pd
import itertools
from itertools import islice
import operator

# python gtex_main.py --fnt ID_MIR_WASH_tpm.txt --fnb ID_Brain_Nerve_Tissue.txt --fna GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt --tg 'MIR6859-1' 'WASH7P' --bfn '_testplot.png'


def main():
    parser = argparse.ArgumentParser(description=
                                     'boxplot gene expression distribution across \
                                      tissue groups or tissue types for a gene')

    parser.add_argument('--fnt', '-file_name_gene_tmp',
                        type=str,
                        help='RNA Seq gene tmp file (gziped)',
                        required=True)

    parser.add_argument('--fnb', '-file_name_annotations',
                        type=str,
                        help='Sample attributes (tissue types) file (txt)',
                        required=True)

    parser.add_argument('--fna', '-file_name_ages',
                        type=str,
                        help='Sample phenotypes (ages) file (txt)',
                        required=True)




    # test file uses genes 'MIR6859-1' and 'WASH7P'
    parser.add_argument('--tg', '-target_genes',
                        type=str, nargs='*',
                        help='Input list of target genes',
                        required =False)

    parser.add_argument('--bfn', '-boxplot_file_name',
                        type=str,
                        help='Name of output boxplot file',
                        required =False)

    args = parser.parse_args()


    tpm_file_name = args.fnt
    brain_tissue_file_name = args.fnb
    age_file_name = args.fna





    #Create brain tissue to id dictionary. result looks like....
    ## = {'Cortex':[id1, id2, id3], 'Medulla':[id4, id5, id6]}
    #I believe the original data was already sorted, so no need to sort
    #short_to_long_id_dict made to fix shortened ids in age file. This dict takes the first...
    ##elements of the id and maps to list of lengthened id. Example...
    ## = {'GTEX-1117F': ['GTEX-1117F-0011-R10a-SM-AHZ7F', 'GTEX-1117F-0011-R10b-SM-CYKQ8',
    ##                                                          'GTEX-1117F-3226-SM-5N9CT']}
    brain_tissue_to_id_dict = {}
    short_to_long_id_dict = {}
    for line in open(brain_tissue_file_name):
        strip = line.rstrip().split('\t')
        tissue = strip[1]
        sample_id = strip[0]
        if tissue not in brain_tissue_to_id_dict:
            brain_tissue_to_id_dict[tissue] = []
        brain_tissue_to_id_dict[tissue].append(sample_id)
        short = '-'.join(sample_id.split('-')[:2])
        if short not in short_to_long_id_dict:
            short_to_long_id_dict[short] = []
        short_to_long_id_dict[short].append(sample_id)
    #print(short_to_long_id_dict)



    #Create age to id dictionary. next(afn) skips header row. result looks like....
    # = {'20-29':[id1, id2, id3], '30-39':[id4, id5, id6]}
    #I believe the original data was already sorted, so no need to sort
    #The short_to_long_id_dict replaces the short id with the long id list made above
    age_to_id_dict = {}
    with open(age_file_name) as afn:
        next(afn)
        for line in afn:
            strip = line.rstrip().split('\t')
            age = strip[2]
            sample_id = strip[0]
            #print(sample_id)
            if age not in age_to_id_dict:
                age_to_id_dict[age] = []
            if sample_id in short_to_long_id_dict:
                #print(sample_id)
                age_to_id_dict[age] = age_to_id_dict[age] + short_to_long_id_dict[sample_id]
    #print(age_to_id_dict)


    #Create nested gene to sample_id to TPM dictionary. result looks like....
    # = {'gene1': {'id1':TPM, 'id2':TPM, 'id3':TPM},
    #    'gene2': {'id4':TPM, 'id5':TPM, 'id6':TPM}}
    gene_to_id_to_tpm_dict = {}
    header = None
    for line in open(tpm_file_name):
        if header == None:
            header = []
            for field in line.rstrip().split('\t'):
                header.append(field)
        strip = line.rstrip().split('\t')
        if strip[1] == 'Description':
            continue
        gene = strip[1]
        gene_to_id_to_tpm_dict[gene] = {}
        for i in range(2, len(strip)):
            gene_to_id_to_tpm_dict[gene][header[i]] = float(strip[i])
    #print(gene_to_tpm_dict['MIR6859-1'])
    #print(header)
    #for x in list(gene_to_tpm_dict)[0:1]:
    #    print(gene_to_tpm_dict[x])
        #print ("key {}, value {} ".format(x, gene_to_tpm_dict[x]))


    # ACL code - rearrange dictionary
    # create nested brain tissue to age to id dictionary:
    # = {'Brain - Frontal Cortex (BA9)': {'60-69': {'GTEX-13QIC-0011-R10a-SM-5O9C7', '30-39': {'GTEX-16YQH-0011-R10a-SM-AHZMD', ...}}}
    brain_tissue_to_age_to_id_dict = {}
    for tissue in brain_tissue_to_id_dict:
        if tissue not in brain_tissue_to_age_to_id_dict:
            brain_tissue_to_age_to_id_dict[tissue] = {}
        for age in age_to_id_dict:
            brain_tissue_to_age_to_id_dict[tissue][age] = set(brain_tissue_to_id_dict[tissue]) & \
                set(age_to_id_dict[age])
    print(brain_tissue_to_age_to_id_dict)

    #Mock age_to_brain_tissue_to_gene_to_tpm_dictionary:
    #{'Brain - Frontal Cortex (BA9)': {'MIR6859-1': {'60-69': [], '50-59': [], '40-49': [], '20-29': [], '30-39': #[], '70-79': []}, 'WASH7P': {'60-69': [], '50-59': [], '40-49': [], '20-29': [], '30-39': [], '70-79': []}} }
    brain_tissue_to_gene_to_age_to_tpm_dictionary = {}
    for tissue in brain_tissue_to_age_to_id_dict:
        if tissue not in brain_tissue_to_gene_to_age_to_tpm_dictionary:
            brain_tissue_to_gene_to_age_to_tpm_dictionary[tissue] = {}
        for gene in gene_to_id_to_tpm_dict:
            if gene not in brain_tissue_to_gene_to_age_to_tpm_dictionary:
                brain_tissue_to_gene_to_age_to_tpm_dictionary[tissue][gene] = {}
            for age in brain_tissue_to_age_to_id_dict[tissue]:
                if age not in brain_tissue_to_gene_to_age_to_tpm_dictionary:
                    brain_tissue_to_gene_to_age_to_tpm_dictionary[tissue][gene][age] = []

    for gene in gene_to_id_to_tpm_dict:
        for tissue in brain_tissue_to_age_to_id_dict:
            for age in brain_tissue_to_age_to_id_dict[tissue]:
                for sample_id in brain_tissue_to_age_to_id_dict[tissue][age]:
                    if sample_id in gene_to_id_to_tpm_dict[gene]:
                        brain_tissue_to_gene_to_age_to_tpm_dictionary[tissue][gene][age].append(gene_to_id_to_tpm_dict[gene][sample_id])
    #print(brain_tissue_to_gene_to_age_to_tpm_dictionary)


    # attempt to make list with all genes with difference between young and old samples

    # significant_genes = []
    # young = '20-29'
    # old = '70-79'
    # for age in age_to_brain_tissue_to_gene_to_tpm_dictionary:
    #     if str(age) == young:
    #
    #     if str(age) == old:
    #         print(age)




    # Alison's code to make boxplots

    # gene_name_list = args.tg  # a list strings
    # boxplot_name_base = args.bfn

    # use gene_name to filter data into appropriate lists for plotting

    # if dictionary goes tissue: gene: age: count can easily loop through all
    # tissue keys, looking for info from the desired gene, then plot each age bracket info

    # if dictionary goes tissue:gene:age:counts
    # for tissue in brain_tissue_to_gene_to_age_to_tpm_dictionary:
    #     for gene in gene_name_list:
    #             if gene in brain_tissue_to_gene_to_age_to_tpm_dictionary[tissue]:
    #             make plot for gene, with multiple boxes side by side for ages

    #         boxplot_name = str(tissue) + boxplot_name_base
    #         title = age
    #         x_axis = tissue
    #         y_axis = 'tpms'
    #         x_ticks = gene_name_list
    #         lists = []
    #         for gene in gene_name_list:
    #             if gene in age_to_brain_tissue_to_gene_to_tpm_dictionary[age][tissue]:
    #                 lists.append(age_to_brain_tissue_to_gene_to_tpm_dictionary[age][tissue][gene])
    #
    #         box_viz.boxplot(boxplot_name, title, x_axis, y_axis, lists, x_ticks)



if __name__ == '__main__':
    main()
