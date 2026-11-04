import pandas as pd
import sys
from scipy import stats
import numpy as np
import re
import statsmodels.stats.multitest as smm
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

# create list of active sequences (wo archaic/modern)
# do unique on list
# filter input for mpranalyze comparative using this list on locus column
activity_folder_path = '/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_quantitative'
comb_df_path = f'{activity_folder_path}/adipo_quant_comb_df_adjusted_fdr.csv'
mpranalyze_output = '/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_quantitative/adipo_results_full_quantitative_nobc_fdr_annot_factor_new_design.txt'

output = pd.read_csv(mpranalyze_output, sep='\t', header=0)
output['oligo'] = output.index
#print("SHOWING MPRAnalyze-quant output", output.head().to_string())

# read in RNA-DNA file
RNA_DNA = pd.read_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/allreps_Filter.txt', usecols=["oligo_column", "RNA_rep1", "DNA_rep1", "RNA_rep2", "DNA_rep2", "RNA_rep3", "DNA_rep3", "RNA_rep_comb", "DNA_rep_comb", "ratio_log_rep1", "ratio_log_rep2", "ratio_log_rep3", "ratio_log_rep_comb", "count_rep1", "count_rep2", "count_rep3", "count_rep_comb"])
        
# merge rna/dna ratio to dataframe
merged_with_ratio = output.merge(RNA_DNA, left_on = "oligo", right_on = "oligo_column", how = "left")
# drop oligo column
merged_with_ratio = merged_with_ratio.drop(['oligo_column'], axis=1)

print("SHOWING FINAL MERGED DF", merged_with_ratio.head())

# add control annotations

merged_with_ratio.loc[merged_with_ratio['oligo'].str.contains('ctrl'), 'control_type'] = 'Ctrl'
merged_with_ratio.loc[merged_with_ratio['oligo'].str.contains('scrambled'), 'control_type'] = 'scrambled_control'
merged_with_ratio['control_type'] = merged_with_ratio['control_type'].fillna(value='No control')

# print number of controls for each type
control_types = ["Ctrl", "scrambled_control"]

for control in control_types:
    rows = merged_with_ratio['control_type'] == control
    row_count = rows.sum()
    print(control, ":    ", row_count)


# add active/non active
merged_with_ratio["activity"] = "non_active"
merged_with_ratio.loc[merged_with_ratio["fdr.mad"] < 0.05, "activity"] = "active"
print(merged_with_ratio[merged_with_ratio['activity'] =="active"])

# add fdr that does not consider controls and oligos with less than 5 DNA conts
#print("MERGED WITH RATIO BEFORE TAKING OUT COUTNS WITH LESS THAN 5", len(merged_with_ratio))
filtered_merged_with_ratio = merged_with_ratio[(merged_with_ratio['DNA_rep_comb'] >= 5) & (merged_with_ratio['control_type'] == 'No control')& (~merged_with_ratio['pval.mad'].isna())]
print(len(filtered_merged_with_ratio.index))
rej, pval_corr = smm.fdrcorrection(filtered_merged_with_ratio["pval.mad"])
print(len(pval_corr))
merged_with_ratio.loc[(merged_with_ratio['DNA_rep_comb'] >= 5) & (merged_with_ratio['control_type'] == 'No control')& (~merged_with_ratio['pval.mad'].isna()), 'fdr.mad_adjusted'] = pval_corr
merged_with_ratio['fdr.mad_adjusted'].fillna(1, inplace=True)

# add new activity columns
merged_with_ratio.loc[merged_with_ratio["fdr.mad_adjusted"] < 0.05, "activity_adjusted"] = "active"
merged_with_ratio.loc[(merged_with_ratio['DNA_rep_comb'] >= 5) & (merged_with_ratio['control_type'] == 'No control') & (~merged_with_ratio['pval.mad'].isna()) & (merged_with_ratio["fdr.mad_adjusted"] >= 0.05), 'activity_adjusted'] = "non_active"



# add column whether oligo is input or not for mpranalyze comp
merged_with_ratio["input_comparative"] = "no"
merged_with_ratio.loc[merged_with_ratio["fdr.mad_adjusted"] < 0.05, "input_comparative"] = "yes"
merged_with_ratio.loc[merged_with_ratio['control_type'] == "PosCtrl_diff", "input_comparative"] = "yes"

merged_with_ratio.to_csv(comb_df_path, header=True, index = False)

#print("quant_results done", "\n", merged_with_ratio)


# use quant output instead of direct mpra output
output_quantitative = pd.read_csv(comb_df_path)



output_quantitative['oligo'] = output_quantitative['oligo'].str.extract(r'(seq\d+)') # remove allele

active_df = output_quantitative[output_quantitative["input_comparative"]=="yes"] 
inactive_df = output_quantitative[output_quantitative["input_comparative"]=="no"] 

active_oligos = list(set(active_df["oligo"].tolist())) #turn into list using only unique oligos
inactive_oligos = list(set(inactive_df["oligo"].tolist())) 

#print(active_df.head().to_string())
#print(active_df.info())
print("Number active oligos",len(active_oligos))

#print(inactive_df.head().to_string())
#print(inactive_df.info())
print("Number inactive oligos",len(inactive_oligos))


# read in input for mpranalyze comparative
input_comparative_rna = pd.read_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/adipo_RNA_MPRAnalyze_comparative_improved.txt', sep='\t', header =0)
input_comparative_dna = pd.read_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/adipo_DNA_MPRAnalyze_comparative_improved.txt', sep='\t', header =0)

print(input_comparative_rna['locus'])
print(active_oligos)
input_comparative_rna_active = input_comparative_rna[input_comparative_rna['locus'].isin(active_oligos)]
input_comparative_dna_active = input_comparative_dna[input_comparative_dna['locus'].isin(active_oligos)]

print(input_comparative_rna_active.info())
print(input_comparative_dna_active.info())

input_comparative_rna_active.to_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/adipo_RNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt', sep='\t', header=True, index=False)
input_comparative_dna_active.to_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/adipo_DNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt', sep='\t', header=True, index=False)
