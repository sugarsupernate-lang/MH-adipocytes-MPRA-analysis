import pandas as pd
import sys
import statsmodels.stats.multitest as smm


# read in concataneted output of mpranalyze comparative for all libraries
output = pd.read_csv('/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/results/adipo_comp_results.tsv', sep='\t', header=0, names = ['oligo','statistic','pval','fdr','df.test','df.dna','df.rna.full','df.rna.red','logFC'])

#TAKE OUT ROWS WHERE WE DONT HAVE FOLDCHANGE
print("LENGTH OF OUTPUT FILE", len(output))
print("Number of rows containing a logFC of NA", output['logFC'].isna().sum())
output = output.dropna(subset=['logFC'])
print("Length after dropping logFC NAs", len(output))

# take out rows that were previously defined as active but now not anymore - basically only take rows that are now defined as active
# create set with sequences defined as active according to new calculation

## RYDER: Don't need to do this step, because this is done for massive MPRA where each library was processed separately
#output_quant_combined = pd.read_csv(f'./{cells}/quantitative_analysis_combined/comb_df_combined_fdr.csv')
#output_quant_combined['oligo'] = output_quant_combined['oligo'].str.replace("_archaic|_modern", '')
#active_df = output_quant_combined[(output_quant_combined['activity_adjusted_combined'] == 'active')] 
#active_oligos = list(set(active_df["oligo"].tolist())) #turn into list using only unique oligos

#output = output[output['oligo'].isin(active_oligos)]
#print(len(output.index))

# do correction for all
rej, pval_corr = smm.fdrcorrection(output["pval"])
output["differential.full"] = rej
output["fdr.full"] = pval_corr
print(len(pval_corr))

# do correlation for all except controls
output_wo_controls = output[~output['oligo'].str.contains('Ctrl')]
print(len(output_wo_controls.index))
rej, pval_corr = smm.fdrcorrection(output_wo_controls["pval"])
print(len(pval_corr))
output.loc[~output['oligo'].str.contains('Ctrl'), 'fdr.wo_controls'] = pval_corr
output.loc[~output['oligo'].str.contains('Ctrl'), 'differential.wo_controls'] = rej

# count differntial/non differential
num_differential = (output["differential.wo_controls"] == True).sum()
num_not_differential = (output["differential.wo_controls"] == False).sum()
print(f'# differential: {num_differential}')
print(f'# not differential: {num_not_differential}')
print(f'# total: {num_differential+num_not_differential}')
print(f'% differential: {num_differential/(num_differential+num_not_differential)}')

# save new df
output.to_csv("/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/results/adipo_MH_comp_fdr_final.csv", header=True, index=False)
