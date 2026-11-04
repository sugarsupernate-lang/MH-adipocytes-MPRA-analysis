# the aim of this script is to apply an outlier filter (2 stds) on barcode_counts_UMI.txt and return a new file in the format of barcode_counts_UMI.txt that can then be used as input for the rest of the pipeline.
# in addition the script should also check the effect of the filter (e.g. how many barcodes are removed, correlation etc)

import pandas as pd
import numpy as np
import sys
#import matplotlib # added to prevent display error #Katharina 28.7.22
#matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
#import matplotlib.pyplot as plt
import os
#import seaborn as sns
from collections import Counter
from statistics import mean 


#todo
#add counts for other libraries x
#delete other filter adn barcode thresholds x
#save exploded df x
#statistics of filter: histogram of how many bcs per oligo are removed, oligo coverage, correlation between reps + average correlation between reps (for different dna count thresholds and activity thresholds)
#scatterpliot of correlation

def std_n_outliers(rna, dna, oligo_bc, n):
    mean = np.mean(rna)
    std = np.std(rna)
    lower_bound = mean - n * std
    upper_bound = mean + n * std
    indices = [i for i, v in enumerate(rna) if lower_bound < v < upper_bound]
    rna = [rna[i] for i in indices]
    dna = [dna[i] for i in indices]
    oligo_bc = [oligo_bc[i] for i in indices]
    return rna, dna, oligo_bc


# manually copied read counts form alignment output
def create_ratio_df():

    rep1 = pd.read_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/UMI/barcode_counts_UMI_rep1.txt', sep='\t', header=0, names=["oligo_bc_no_filter", "RNA_no_filter", "DNA_no_filter"])
    rep2 = pd.read_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/UMI/barcode_counts_UMI_rep2.txt', sep='\t', header=0, names=["oligo_bc_no_filter", "RNA_no_filter", "DNA_no_filter"])
    rep3 = pd.read_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/UMI/barcode_counts_UMI_rep3.txt', sep='\t', header=0, names=["oligo_bc_no_filter", "RNA_no_filter", "DNA_no_filter"])

    UMI_all_reps = pd.read_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/UMI/barcode_counts_UMI.txt', sep='\t')
    UMI_all_reps.set_index('oligo_bc', inplace=True)

    rep1DNA = rep1["DNA_no_filter"].sum()
    rep2DNA = rep2["DNA_no_filter"].sum()
    rep3DNA = rep3["DNA_no_filter"].sum()
    rep1RNA = rep1["RNA_no_filter"].sum()
    rep2RNA = rep2["RNA_no_filter"].sum()
    rep3RNA = rep3["RNA_no_filter"].sum()


    # get number of reads 
    DNA_reads = [rep1DNA, rep2DNA, rep3DNA]
    RNA_reads = [rep1RNA, rep2RNA, rep3RNA]
    read_dict = {"DNA_reads": DNA_reads, "RNA_reads": RNA_reads}
    print(DNA_reads)
    print(RNA_reads)

    # group oligos and sum RNA, DNA counts. Also counts how many barcodes we have per oligo. Calculate ratio and perform log2 on ratio
    replicates = []
    exploded_df = []
    
    for i, rep in enumerate([rep1, rep2, rep3]):
        rep[f'oligo_bc_filtered_std2'] = None
        for molecule in ["RNA", "DNA"]:
            rep[f'{molecule}_filtered_std2'] = None
        rep = rep[(rep['DNA_no_filter'] > 0) | (rep['RNA_no_filter'] > 0)]
        rep[['oligo', 'bc']] = rep.oligo_bc_no_filter.str.rsplit("_", 1, expand=True)
        rep = rep.set_index('oligo')
        rep_list = rep.groupby('oligo').aggregate(lambda a: a.tolist())
        rep_list = rep_list.reset_index()
        
        for index, row in rep_list.iterrows():
            rep_list.at[index, 'RNA_filtered_std2'], rep_list.at[index, 'DNA_filtered_std2'], rep_list.at[index, 'oligo_bc_filtered_std2'] = std_n_outliers(
                row['RNA_no_filter'], row['DNA_no_filter'], row['oligo_bc_no_filter'], 2)

        for outlier_filter in ["no_filter","filtered_std2"]:
            to_explode = rep_list[[col for col in rep_list.columns if f"{outlier_filter}" in col]].copy()
            to_explode = to_explode[to_explode[f'oligo_bc_{outlier_filter}'].astype(bool)]
            exploded = to_explode.explode([f"RNA_{outlier_filter}", f"DNA_{outlier_filter}", f'oligo_bc_{outlier_filter}'], ignore_index=True)
            exploded_df.append(exploded)
            rep_list = rep_list.drop([f'oligo_bc_{outlier_filter}'], axis=1)
        rep_list = rep_list.drop(['bc'], axis=1)
    
    
        replicates.append(rep_list)

    rep1_counts = replicates[0]
    rep2_counts = replicates[1]
    rep3_counts = replicates[2]

    # merge the different outputs
    # save the exploded df (2std) - so we get a file in the same format as barcode_counts_UMI.txt - which than can be used to create the inout for MPRAnalyze and the other steps of the analysis (get one file for all the replicates)
    
    rep1_exploded_filtered = exploded_df[1]
    rep2_exploded_filtered = exploded_df[3]
    rep3_exploded_filtered = exploded_df[5]
    exploded_df_filter = rep1_exploded_filtered.merge(rep2_exploded_filtered,on='oligo_bc_filtered_std2',how="outer",suffixes=("_rep1",None)).merge(rep3_exploded_filtered,on='oligo_bc_filtered_std2',how="outer",suffixes=("_rep2","_rep3"))
    # print(list(exploded_df_filter.columns))
    exploded_df_filter.fillna(0.0, inplace=True)
    exploded_df_filter = exploded_df_filter[["oligo_bc_filtered_std2", "RNA_filtered_std2_rep1", "RNA_filtered_std2_rep2","RNA_filtered_std2_rep3","DNA_filtered_std2_rep1","DNA_filtered_std2_rep2","DNA_filtered_std2_rep3"]]
    exploded_df_filter.columns = ["oligo_bc","RNA_rep1","RNA_rep2","RNA_rep3","DNA_rep1","DNA_rep2","DNA_rep3"]
    exploded_df_filter.set_index('oligo_bc', inplace=True)
    UMI_comb_filter = exploded_df_filter.combine_first(UMI_all_reps)
    UMI_comb_filter.reset_index(inplace=True)
    UMI_comb_filter.to_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/UMI/UMI_std2_filter.txt', index=False, sep='\t')

    rep1_exploded_no_filter = exploded_df[0]
    rep2_exploded_no_filter = exploded_df[2]
    rep3_exploded_no_filter = exploded_df[4]
    
    exploded_df_no_filter = rep1_exploded_no_filter.merge(rep2_exploded_no_filter,on='oligo_bc_no_filter',how="outer",suffixes=("_rep1",None)).merge(rep3_exploded_no_filter,on='oligo_bc_no_filter',how="outer",suffixes=("_rep2","_rep3"))
    exploded_df_no_filter.fillna(0.0, inplace=True)
    exploded_df_no_filter = exploded_df_no_filter[["oligo_bc_no_filter", "RNA_no_filter_rep1", "RNA_no_filter_rep2","RNA_no_filter_rep3","DNA_no_filter_rep1","DNA_no_filter_rep2","DNA_no_filter_rep3"]]
    exploded_df_no_filter.columns = ["oligo_bc","RNA_rep1","RNA_rep2","RNA_rep3","DNA_rep1","DNA_rep2","DNA_rep3"]
    exploded_df_no_filter.set_index('oligo_bc', inplace=True)
    UMI_comb_no_filter = exploded_df_no_filter.combine_first(UMI_all_reps)
    UMI_comb_no_filter.reset_index(inplace=True)
    UMI_comb_no_filter.to_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/UMI/UMI_no_filter.txt', index=False, sep='\t')

    return

ratio_df = create_ratio_df()
    
    
