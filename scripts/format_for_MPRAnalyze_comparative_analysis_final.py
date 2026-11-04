#!/usr/bin/env python

import os
import pysam
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import sys
import itertools 
import pickle
import random

#library = sys.argv[1]
#adaptor = sys.argv[2]
cells = "osteoblast"

# wynton paths
original_oligo_list = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/associations/MH_oligos.fa"
data_df_path = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/UMI/UMI_std2_filter.txt"
pickle_path = '/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/MPRAnalyze_dict_nocontrols_improved2.pickle'

# local paths
#original_oligo_list = "/Users/reasterlin/Desktop/ahituv_lab/AH_MPRA/final_library_files/archaic_derived_mpra_oligos.fa"
#data_df_path = "/Users/reasterlin/Desktop/ahituv_lab/AH_MPRA_experimental/HOB_AH_data/barcode_counts_UMI_small.txt"
#pickle_path = '/Users/reasterlin/Desktop/ahituv_lab/AH_MPRA_experimental/HOB_AH_data/mpranalyze_comparative/MPRAnalyze_dict_nocontrols_improved.pickle'


# SF: This script runs on barcode_counts_UMI.txt. In this file each row is an oligo-barcode pair, and each column is a count.
# Number of columns is DNA/RNA*cell type * rep so in Carly's it is 2*3*3=18 cols, in the human mpra 2*3 = 6 (dna/rna * rep)
# 1. removing controls
# 2. Calculate the oligo with the most barcodes, per allele source (separately for "Ancestral" or "Derived")
# 3. Make a table using the very first original oligos, create a list of 'loci' that has the ancestral or derived categorization removed.
# why do we have to use the original fasta here and not the same data as in step 1+2?
# 4. Initialize two dfs, with each row being a 'locus', no columns yet. base for MPRAnalyze dataframe
# 5. create the columns - each column is a unique set of cell-rep-allele type. Different barcodes for the same set are populated in different columns

def main():    
    ### prep the data for looping ###
    
    def turn_into_locus(f):
        oligo_split = f.split("_")
        locus = oligo_split[1] 
        return locus

    def turn_into_locus_on_fasta(f): # have to use different function because fasta does not contain barcode
        oligo_split = f.split("_")
        locus = oligo_split[1]
        return locus
        
    def get_allele_source(f):
        if "Derived" in f: allele = "modern"
        elif "Ancestral" in f: allele = "archaic"
        return allele
    
    def get_bcseq(f):
        oligo_split = f.split("_")    
        bc_seq = oligo_split[-1]   # 15 bp bc seq is the last item always
        return bc_seq
    
    def get_oligo(f):
        oligo_split = f.split("_")
        oligo_to_fuse = oligo_split[:-1]   # remove the bc_seq, always the last item in the name
        oligo = "_".join(oligo_to_fuse)
        return oligo
    
    data_df = pd.read_csv(data_df_path, sep='\t', index_col = 'oligo_bc')
    
    print(len(data_df.index), 'number of barcodes in original df')
    
    #remove ctrls from the original dataframe by creating list without them, and dropping those index rows
    indices = list(data_df.index.values)
    ctrl= [index for index in indices if "ctrl" in index]
    scrambled_control = [index for index in indices if 'scrambled' in index]

    all_ctrls =  ctrl+scrambled_control
    print(len(all_ctrls), 'number of rows to drop')
    data_df.drop(all_ctrls, inplace=True)
    #data_df.drop(pos_controls, inplace=True)
    print(len(data_df.index), 'number of barcodes after dropping controls')
    
    #add easily accessible info for locus, allele source and bcsequence in new columns
    data_df['locus'] = data_df.index.map(turn_into_locus)
    data_df['allele_source'] = data_df.index.map(get_allele_source)
    data_df['bc_seq'] = data_df.index.map(get_bcseq)
    data_df['oligo'] = data_df.index.map(get_oligo)
    
    #figure out the oligo with the most barcodes, per allele source
    max_bc_info_dict = {'archaic': [], 'modern': []}
    df_allele_grouped = data_df.groupby('allele_source')
    for each_allele, allele_group in df_allele_grouped:
        max_barcodes = 0
        max_oligo_id = 'test'
        df_oligo_grouped = allele_group.groupby('oligo')
        for each_oligo, barcodes in df_oligo_grouped:
            number_barcodes = len(barcodes.index)
            if number_barcodes > max_barcodes:
                max_barcodes = number_barcodes
                max_oligo_id = each_oligo
        max_bc_info_dict[each_allele] = max_barcodes
        #print('For the', each_allele, 'allele', max_oligo_id, 'had the most barcodes, with', max_barcodes, 'barcodes.')
    
    # make a table using the very first original oligos, create a list of 'loci' that has the ancestral or derived categorization removed
    fasta = pysam.FastaFile(original_oligo_list)
    fasta_reference_names = fasta.references
    loci = []   # generate a list of all the loci (oligos without ancestral or derived info)
    for each_oligo in fasta_reference_names:
        if 'ctrl' in each_oligo:    
            continue
        if 'scrambled' in each_oligo: 
            continue

        locus = turn_into_locus_on_fasta(each_oligo)
        loci.append(locus)
    no_dups_loci = list(dict.fromkeys(loci))  # there will be dups because of ancestral or derived for each
    print('Length of loci: ', len(loci))
    print('Length of no_dups_loci: ', len(no_dups_loci))

    new_df_RNA = pd.DataFrame(index=no_dups_loci)  # initialize two dfs, with each row being a 'locus', no columns yet. base for MPRAnalyze dataframe
    new_df_DNA = pd.DataFrame(index=no_dups_loci)
    
    meta_info = data_df.loc[:, ['locus', 'allele_source', 'bc_seq', 'oligo']]
    
    # loop over each cell type + rep combo, make new df and then concat together at the end
        
    cell_dfs = {}
    sanity_check_1 = {}
    for each_rep in ['rep1', 'rep2', 'rep3']:
        print('rep: ', each_rep)

        rep_data = data_df.copy().filter(regex=each_rep)   # outputs dataframe for each rep/celltype, includes RNA and DNA
        merged_data = rep_data.merge(meta_info, on='oligo_bc')
        
        for each_allele in ['modern', 'archaic']:
            print('each_allele: ', each_allele)
            allele_df = merged_data[merged_data['allele_source'] == each_allele]   # only the rows where the allele_source col matches the looped allele
            num_cols = max_bc_info_dict[each_allele]
            identifier = each_rep + '_' + each_allele
            identifier_wo_allele = each_rep
            # add columns with number from 1-max col number, + identifier that includes cell type, rep and allele
            print(identifier, num_cols)
            col_names = [ identifier+'_'+str(num) for num in list(range(1,num_cols+1)) ]

            # make new dfs for RNA and DNA, add the number of columns found earlier of max number of ancestral or derived barcodes
            RNA_df = new_df_RNA.copy().reindex(columns=col_names)
            DNA_df = new_df_DNA.copy().reindex(columns=col_names)
            # get list of all the loci to loop over that we have data for (there will be repeats since AH and MH oligos)
            allele_df = allele_df.set_index('locus')
            allele_df_list = allele_df.groupby('locus').aggregate(lambda a: a.tolist())
            allele_df_list = allele_df_list.reset_index()
            all_loci = allele_df.index.tolist()  # get list of all the oligo to loop over that we have data for (there will be repeats since AH and MH oligos)            # print('all loci: ',len(all_loci), all_loci)
            # unique list of all the loci we have data for
            all_loci_unique = list(dict.fromkeys(all_loci))
            # print('all unique loci: ',len(all_loci_unique), all_loci_unique)

            sanity_check = {}
            for i, each_locus in enumerate(all_loci_unique):

                if i % 10000 == 0: print(i)

                locus_df = allele_df_list[allele_df_list['locus'] == each_locus]   # get the rows of the df that are from that locus
                sanity_check[each_locus] = len(locus_df.iloc[0]['bc_seq'])
                # all barcodes RNA values per locus is converted into a list
                # iloc[:,0] is all the rows of the first column (there's only one column from the filtering)
                RNA_data = locus_df['RNA' + '_' + identifier_wo_allele].values[0]
                DNA_data = locus_df['DNA' + '_' + identifier_wo_allele].values[0]  # same for DNA

                # pad the list with zeroes to the same length as the number of columns in the df
                RNA_data.extend([0] * (num_cols - len(RNA_data)))

                # slot the list as a row in the new df
                RNA_df.loc[each_locus,:] = RNA_data

                DNA_data.extend([0] * (num_cols - len(DNA_data)))
                DNA_df.loc[each_locus,:] = DNA_data

            sanity_check_1[identifier] = sanity_check
            cell_dfs[identifier] = [RNA_df, DNA_df]   # save the dfs in a dict under cell_type, rep# and allele
    
    #print(sanity_check_1)
    if cells == 'osteoblast':
        test1 = sanity_check_1['rep1_modern'] # gives the sanity chekc dict belonging to this ID #KL
        test2 = sanity_check_1['rep2_modern']
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

    with open(pickle_path, 'wb') as handle:
        pickle.dump(cell_dfs, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    return cell_dfs

def cell_types_separate(cell_dfs):    
    
    def return_rep(f):
        split = f.split("_") 
        rep = split[0]
        return rep
    def return_allele(f):
        split = f.split("_")
        allele = split[1]
        return allele
    def return_barcode(f):
        split = f.split("_")
        barcode = split[2]
        return barcode
    
    def make_annot_table(df):
        col_names = list(df.columns.values)   # cells + '_' + each_rep + '_' + each_allele + '_' num
        annot_table = pd.DataFrame(index = col_names, columns = ['replicate', 'allele', 'barcode', 'barcode_allele'])
        annot_table['replicate'] = annot_table.index.map(return_rep)
        annot_table['allele'] = annot_table.index.map(return_allele)
        annot_table['barcode'] = annot_table.index.map(return_barcode)
        annot_table['barcode_allele'] = annot_table.allele.str.cat(annot_table.barcode, sep='.')
        #annot_table.to_csv('/Users/reasterlin/Desktop/ahituv_lab/AH_MPRA_experimental/HOB_AH_data/HOb_col_annot_comparative_improved.txt', sep='\t', header=True) local
        annot_table.to_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/adipo_col_annot_comparative_improved.txt', sep='\t', header=True) # wynton
    
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
    
    # local paths 
    #final_df_RNA.to_csv('/Users/reasterlin/Desktop/ahituv_lab/AH_MPRA_experimental/HOB_AH_data/HOb_RNA_MPRAnalyze_comparative_improved.txt', sep='\t', header=True, index_label='locus')
    #final_df_DNA.to_csv('/Users/reasterlin/Desktop/ahituv_lab/AH_MPRA_experimental/HOB_AH_data/HOb_DNA_MPRAnalyze_comparative_improved.txt', sep='\t', header=True, index_label='locus')
    
    # wynton paths
    final_df_RNA.to_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/adipo_RNA_MPRAnalyze_comparative_improved.txt', sep='\t', header=True, index_label='locus')
    final_df_DNA.to_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/adipo_DNA_MPRAnalyze_comparative_improved.txt', sep='\t', header=True, index_label='locus')
    

cell_dfs =  main()
print('onto the second step')
# run one of these  
cell_types_separate(cell_dfs)
#reps_divided(cell_dfs)
