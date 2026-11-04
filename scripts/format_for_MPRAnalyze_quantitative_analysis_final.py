#!/usr/bin/env python
#

import os
import pysam
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import sys
import itertools 
import pickle
import random

cells="HOb"  # NOT the cell type. Only gates the rep1-vs-rep2 oligo-count sanity
             # check below (prints only, no effect on output). Left as "HOb" on
             # purpose so the check still runs for adipo — see PIPELINE_PLAN.md.


# local paths for Ry

# wynton paths
original_oligo_list = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/associations/MH_oligos.fa" 
data_df_path = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/UMI/UMI_std2_filter.txt"


# SF 29Jan2023: moved functions out of main function
def get_bcseq(f):
    oligo_split = f.split("_")
    bc_seq = oligo_split[-1]  # 15 bp bc seq is the last item always
    return bc_seq

def get_oligo(f):
    oligo_split = f.split("_")
    oligo_to_fuse = oligo_split[:-1]  # remove the bc_seq, always the last item in the name
    oligo = "_".join(oligo_to_fuse)
    return oligo

def main():    
    ### prep the data for looping ###

    df = pd.read_csv(data_df_path, sep='\t', index_col = 'oligo_bc')

    # add easily accessible info for oligo and bcsequence in new columns
    df['bc_seq'] = df.index.map(get_bcseq)
    df['oligo'] = df.index.map(get_oligo)

    #figure out the oligo with the most barcodes
    max_barcodes = 0
    df_oligo_grouped = df.groupby('oligo')
    for each_oligo, barcodes in df_oligo_grouped:
        number_barcodes = len(barcodes.index)
        if number_barcodes > max_barcodes:
            max_barcodes = number_barcodes
    #print('For the', each_allele, 'allele', max_oligo_id, 'had the most barcodes, with', max_barcodes, 'barcodes.')

    # make a table using the very first original oligos, create a list of 'loci'
    fasta = pysam.FastaFile(original_oligo_list)
    fasta_reference_names = fasta.references
    new_df_RNA = pd.DataFrame(index=fasta_reference_names)  # initialize two dfs, with each row being an 'oligo', no columns yet. base for MPRAnalyze dataframe
    new_df_DNA = pd.DataFrame(index=fasta_reference_names)

    meta_info = df.loc[:, ['bc_seq', 'oligo']]

    # loop over each rep, make new df and then concat together at the end

    cell_dfs = {}
    sanity_check_1 = {}
    for each_rep in ['rep1', 'rep2', 'rep3']:
        print(each_rep)
        rep_data = df.copy().filter(regex=each_rep)   # outputs dataframe for each rep/celltype, includes RNA and DNA
        merged_data = rep_data.merge(meta_info, on='oligo_bc')   # re-add the meta info
        merged_data = merged_data.set_index('oligo')
        merged_data_list = merged_data.groupby('oligo').aggregate(lambda a: a.tolist())
        merged_data_list = merged_data_list.reset_index()
        all_oligo = merged_data.index.tolist()  # get list of all the oligo to loop over that we have data for (there will be repeats since AH and MH oligos)
        identifier = each_rep
        col_names = [identifier+'_'+str(num) for num in list(range(1,max_barcodes+1))]   # add columns with number from 1-max col number, + identifier that includes cell type and rep
        RNA_df = new_df_RNA.copy().reindex(columns=col_names)  # make new dfs for RNA and DNA, add the number of columns found earlier of max number of ancestral or derived barcodes
        DNA_df = new_df_DNA.copy().reindex(columns=col_names)
        all_oligo_unique = list(dict.fromkeys(all_oligo))  # unique list of all the oligo we have data for
        print(len(all_oligo_unique))
        sanity_check = {}
        for i, each_oligo in enumerate(all_oligo_unique):
            if i % 10000 == 0: print(i)
            oligo_df = merged_data_list[merged_data_list['oligo'] == each_oligo]  # get the rows of the df that are from that locus
            sanity_check[each_oligo] = len(oligo_df.iloc[0]['bc_seq'])  # for each oligo we save the number of barcodes
            RNA_data = oligo_df['RNA' + '_' + identifier].values[0]  # all barcodes RNA values per locus is converted into a list    iloc[:,0] is all the rows of the first column (there's only one column from the filtering)
            DNA_data = oligo_df['DNA' + '_' + identifier].values[0] # same for DNA
            RNA_data.extend([0] * (max_barcodes - len(RNA_data)))  # pad the list with zeroes to the same length as the number of columns in the df
            RNA_df.loc[each_oligo,:] = RNA_data  # slot the list as a row in the new df
            DNA_data.extend([0] * (max_barcodes - len(DNA_data)))
            DNA_df.loc[each_oligo,:] = DNA_data
        sanity_check_1[identifier] = sanity_check

        cell_dfs[identifier] = [RNA_df, DNA_df]   # save the dfs in a dict under cell_type, rep# and allele

    # print(sanity_check_1)
    if cells == 'HOb':
        test1 = sanity_check_1['rep1'] # gives the sanity chekc dict belonging to this ID #KL
        test2 = sanity_check_1['rep2']
        len1 = len(test1) # gives us the number of oligos #KL
        len2 = len(test2)
        if len1 != len2:
            print('dicts not the same length')
        for num in range(1,10):
            locus = random.choice(list(test1)) #gives us random oligo #KL
            num1 = test1[locus] #gives us length of oligo_df belonging to this oligo #KL
            num2 = test2[locus]
            if num1 != num2: #not sure why this is used as sanity check #KL
                print('something went wrong')
            if num1 == num2:
                print('everything fine...for now')
    
    return cell_dfs

def cell_types_separate(cell_dfs):
    def return_rep(f):
        split = f.split("_") 
        rep = split[0]
        return rep
        
    def return_barcode(f):
        split = f.split("_")
        barcode = split[1]
        return barcode
    
    def make_annot_table(df):
        col_names = list(df.columns.values)   # cells + '_' + each_rep + '_' + each_allele + '_' num
        annot_table = pd.DataFrame(index = col_names, columns = ['replicate', 'barcode'])
        annot_table['replicate'] = annot_table.index.map(return_rep)
        annot_table['barcode'] = annot_table.index.map(return_barcode)
        #annot_table.to_csv(f'./{cells}/{library}{adaptor}/output/mpranalyze_quantitative/{cells}_{library}{adaptor}_col_annot_quantitative.txt', sep='\t', header=True)
        annot_table.to_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_quantitative/adipo_col_annot_quantitative.txt', sep='\t', header=True)

    to_fuse_RNA = []
    to_fuse_DNA = []
    for each_id, list_of_dfs in cell_dfs.items():
        to_fuse_RNA.append(list_of_dfs[0])
        to_fuse_DNA.append(list_of_dfs[1])
    final_df_RNA = pd.concat(to_fuse_RNA, axis=1)
    final_df_DNA = pd.concat(to_fuse_DNA, axis=1)
    final_df_RNA.fillna(0, inplace=True)
    final_df_DNA.fillna(0, inplace=True)
    
    shape_RNA = final_df_RNA.shape
    shape_DNA = final_df_DNA.shape
    
    if not shape_RNA == shape_DNA:   # the RNA and DNA matrices need to be the exact same size (same number of columns and rows), check here
        print('something went wrong and you need to fix that bug')
    
    make_annot_table(final_df_RNA)   # the same for both, only needs to be done once, just feed RNA in (could have been the DNA table too)
    
    final_df_RNA.to_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_quantitative/adipo_RNA_MPRAnalyze_quantitative.txt',
        sep='\t', header=True, index_label='locus')
    final_df_DNA.to_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_quantitative/adipo_DNA_MPRAnalyze_quantitative.txt',
        sep='\t', header=True, index_label='locus')

cell_dfs =  main()
print('onto the second step')
cell_types_separate(cell_dfs)

