#!/usr/bin/env python
"""
adipo MH — sort_comp_input_barcodes.py (canonical, RNA/DNA swap already fixed).

BUG in the original: the two to_csv() calls at the end had filenames swapped:
  input_comparative_rna_to_save.to_csv('.../HOb_DNA_..._sorted.txt', ...)  # WRONG
  input_comparative_dna_to_save.to_csv('.../HOb_RNA_..._sorted.txt', ...)  # WRONG

This version writes each dataframe to the correctly-labeled file. All other
logic is unchanged.

Reads (already filtered) from:
  output/mpranalyze_comparative_rerun/{adipo_RNA,adipo_DNA}_MPRAnalyze_comparative_filter_adjusted_fdr.txt

Writes sorted files to the same directory:
  adipo_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt
  adipo_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt
"""
import pandas as pd
import os

BASE = '/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative_rerun'

input_comparative_dna = pd.read_csv(f'{BASE}/adipo_DNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt', sep='\t', header=0)
input_comparative_rna = pd.read_csv(f'{BASE}/adipo_RNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt', sep='\t', header=0)

# Count non-zero barcodes per oligo
input_comparative_rna["count_bcs_rna"] = (input_comparative_rna != 0).sum(axis=1)
input_comparative_dna["count_bcs_dna"] = (input_comparative_dna != 0).sum(axis=1)
input_comparative_dna["count_bcs_rna"] = input_comparative_rna["count_bcs_rna"]

print("SORTING")
input_comparative_rna.sort_values(by='count_bcs_rna', ascending=False, inplace=True)
input_comparative_dna.sort_values(by='count_bcs_rna', ascending=False, inplace=True)

print(input_comparative_rna[['locus', 'count_bcs_rna']].head(10))
print(input_comparative_dna[['locus', 'count_bcs_dna', 'count_bcs_rna']].head(10))

input_comparative_rna_to_save = input_comparative_rna.drop(columns=['count_bcs_rna'])
input_comparative_dna_to_save = input_comparative_dna.drop(columns=['count_bcs_dna', 'count_bcs_rna'])

# >>> FIX: RNA data goes to RNA-named file, DNA data goes to DNA-named file <<<
input_comparative_rna_to_save.to_csv(
    f'{BASE}/adipo_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt',
    sep='\t', header=False, index=False)
input_comparative_dna_to_save.to_csv(
    f'{BASE}/adipo_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt',
    sep='\t', header=False, index=False)

# Compute split positions for csplit (unchanged from original)
count_bcs_rna_list = input_comparative_rna['count_bcs_rna'].tolist()
chunk_sum = 0
split_indices = []
max_bc = 8000

for idx, count_bcs in enumerate(count_bcs_rna_list):
    if chunk_sum + count_bcs > max_bc:
        if idx != 0:
            split_indices.append(idx)
        chunk_sum = 0
    chunk_sum += count_bcs

split_indices_plus_one = [x + 1 for x in split_indices]
print(f"n chunks: {len(split_indices_plus_one) + 1}")

with open(f'{BASE}/split_positions.txt', 'w') as f:
    for line in split_indices_plus_one:
        f.write(f"{line} ")
