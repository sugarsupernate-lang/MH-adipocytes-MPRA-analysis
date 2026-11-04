import os
import pysam
import pandas as pd
import numpy as np
import matplotlib # added to prevent display error #Katharina 4.8.22
matplotlib.use('Agg') # added to prevent display error #Katharina 4.8.22
import matplotlib.pyplot as plt
import sys
import itertools 

#library = sys.argv[1]
#adaptor = sys.argv[2]
#cells = sys.argv[3]
#min_associations = sys.argv[4]

# the path with the alignments of R1 and R3 to the barcodes (oligo to barcode pairs)
#alignment_dir = f'./{cells}/{library}{adaptor}/output/alignments'
alignment_dir = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/aligned_reads"

#project_dir = f'./{cells}/{library}{adaptor}/output/UMI'
project_dir = '/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/UMI'

# the path with the R2s of each experiment
#UMI_fastq_dir = f'./{cells}/{library}{adaptor}/input/DNA_RNA'
UMI_fastq_dir = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/fastq_data/UMI_reads"

# oligo to barcode pairs fasta from the DNA only part
#barcode_mapping_fasta = f'./neurons/{library}{adaptor}/output/DNA_barcode_associations_{min_associations}/oligos_to_barcodes_comb_{library}{adaptor}_30_{min_associations}.fasta' #katharina 6.3.23 - association output is saved in neuron folder
barcode_mapping_fasta = "/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/associations/oligos_barcodes.fa"
# data_df_path = f'{project_dir}/barcode_counts_comb_UMI.txt'
data_df_path = f'{project_dir}/barcode_counts_UMI.txt' # Katharina 4.8.22 

minimum_mapping_quality = 20   
cigar_req = '15M'
MD_req = '15'

def loop_over_bam(bam, fastq):
    
    not_proper_pair = 0
    cigar_pairs = 0
    MD_pairs = 0
    number_pairs_not_passing_mapq = 0
    total_passing_filter = 0
    barcode_UMI = {}
    barcode_counts = {}
    mol_dict = {}

    for each_mol in fastq:
        mol = each_mol.name
        UMI = each_mol.sequence
        mol_dict[mol] = UMI
    
    for i, read1 in enumerate(bam):
        
        if not read1.is_read1:  # go through reads as pairs only, skip every other read (don't look at read 2 if read 1 not present)
            continue
        
        if not read1.is_proper_pair:  # filter out cases that aren't 'proper pairs' (make sure paired-endedness worked ok?)
            not_proper_pair += 1
            continue
        
        read2 = next(bam)    
        if not read2.is_read2:
            continue # make sure that we're operating on pairs of reads
        
        read1_cigar = read1.cigarstring 
        read1_MD = read1.get_tag('MD')
        read2_cigar = read2.cigarstring
        read2_MD = read2.get_tag('MD')
            
        read1_mapq = read1.mapping_quality
        read2_mapq = read2.mapping_quality
        
        if not read1_cigar == cigar_req or not read2_cigar == cigar_req:
            cigar_pairs += 1
            continue
        if not read1_MD == MD_req and not read2_MD == MD_req:
            MD_pairs += 1
            continue
        if not read1_mapq >= minimum_mapping_quality or not read2_mapq >= minimum_mapping_quality:   # need at least one of the reads to have mapq greater than threshold  ### 2020-3-11 that's not what this does. it requires both reads to have mapq of 20 or more
            number_pairs_not_passing_mapq += 1
            continue
        
        total_passing_filter += 1
        
        oligo_bcseq = read1.reference_name
        molecule = read1.query_name
        UMI_seq = mol_dict[molecule]
        
        if oligo_bcseq in barcode_UMI:
            barcode_UMI[oligo_bcseq].append(UMI_seq)
        if oligo_bcseq not in barcode_UMI:
            barcode_UMI[oligo_bcseq] = [UMI_seq]
    
    how_diff_UMI = []
    
    for each_bc, UMIs in barcode_UMI.items():
        number_UMI_total = len(UMIs)
        number_unique_UMI = len(set(UMIs))
        percent_diff = float(number_unique_UMI / number_UMI_total)
        how_diff_UMI.append(percent_diff)
        barcode_counts[each_bc] = number_unique_UMI
    
    print(not_proper_pair, 'pairs filtered from not proper pair')
    print(cigar_pairs, 'pairs filtered for cigar')
    print(MD_pairs, 'pairs filtered for MD')
    print(number_pairs_not_passing_mapq, 'pairs filtered for mapq')
    print(total_passing_filter, 'total read pairs passing filter')
    
    return barcode_counts, how_diff_UMI
        
def loop_through_bams():
    
    fasta = pysam.FastaFile(barcode_mapping_fasta)
    barcode_coverage_per_file = {}
    fasta_reference_names = fasta.references
    total_number_bc = len(fasta_reference_names)
    print('total_number_bc: ', total_number_bc)
    all_data_df = pd.DataFrame(index=fasta_reference_names)
    UMI_fastqs = os.listdir(UMI_fastq_dir)
    UMI_diffs = []
    print(os.listdir(alignment_dir)) 
    for each_file in os.listdir(alignment_dir):
        if '.bam' in each_file:  # skip the fasta ref file
            print(each_file)
            with open(f'{alignment_dir}/{each_file}', 'r') as f:
                # for a given alignment file - find the corresponding UMI (R2) file
                if each_file == "DNA_aligned_rep1.bam": UMI_fastq_filename = "DNA1_R2.fastq.gz"; rep_info = "DNA_rep1"
                elif each_file == "DNA_aligned_rep2.bam": UMI_fastq_filename = "DNA2_R2.fastq.gz"; rep_info = "DNA_rep2"
                elif each_file == "DNA_aligned_rep3.bam": UMI_fastq_filename = "DNA3_R2.fastq.gz"; rep_info = "DNA_rep3"
                elif each_file == "RNA_aligned_rep1.bam": UMI_fastq_filename = "RNA1_R2.fastq.gz"; rep_info = "RNA_rep1"
                elif each_file == "RNA_aligned_rep2.bam": UMI_fastq_filename = "RNA2_R2.fastq.gz"; rep_info = "RNA_rep2"
                elif each_file == "RNA_aligned_rep3.bam": UMI_fastq_filename = "RNA3_R2.fastq.gz"; rep_info = "RNA_rep3"
                else: "ERROR: could not identify UMI file from alignment filename"

                #rep_info = each_file[8:-4]
                #rep_info_parsed = rep_info.split("_")   #  0 - RNA or DNA, 1 - cell, 2 - rep#
                #print(rep_info)
                #print(rep_info_parsed[0])
                #print(rep_info_parsed[1])
                #print(rep_info_parsed[2])
                #UMI_fastq_filename_list = [name for name in UMI_fastqs if rep_info_parsed[0] in name and rep_info_parsed[2] in name and "R2" in name] # should only be 1 item
                #UMI_fastq_filename = UMI_fastq_filename_list[0]
                
                print(UMI_fastq_filename)
                UMI_fastq = pysam.FastxFile(f'{UMI_fastq_dir}/{UMI_fastq_filename}')
                
                bam_file = pysam.AlignmentFile(f)
                barcode_count_dict, how_diff_UMI = loop_over_bam(bam_file, UMI_fastq)
                UMI_diffs.extend(how_diff_UMI)
                
                all_data_df[rep_info] = all_data_df.index.map(barcode_count_dict)
                all_data_df.fillna(0, inplace=True)
                number_bc_found = float((all_data_df!=0)[rep_info].astype(int).sum(axis=0))
                barcode_coverage_per_file[rep_info] = [number_bc_found, total_number_bc]  
    
    plt.clf()
    plt.hist(UMI_diffs, bins=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    plt.savefig(f'{project_dir}/how_diff_UMI_hist.pdf')
    
    reordered_df = reorder_df(all_data_df)
#    reordered_df=all_data_df # SF for debug of one sample
    reordered_df.to_csv(f'{project_dir}/barcode_counts_UMI.txt', sep='\t', header=True, index_label='oligo_bc')
    
    df_rep1 = reordered_df.filter(regex='rep1')
    df_rep1.to_csv(f'{project_dir}/barcode_counts_UMI_rep1.txt', sep='\t', header=True, index_label='oligo_bc')

    df_rep2 = reordered_df.filter(regex='rep2')
    df_rep2.to_csv(f'{project_dir}/barcode_counts_UMI_rep2.txt', sep='\t', header=True, index_label='oligo_bc')

    df_rep3 = reordered_df.filter(regex='rep3')
    df_rep3.to_csv(f'{project_dir}/barcode_counts_UMI_rep3.txt', sep='\t', header=True, index_label='oligo_bc')

    return barcode_coverage_per_file, reordered_df
    
def reorder_df(df):
    
    new_order = [
    f'RNA_rep1', 
    f'RNA_rep2',
    f'RNA_rep3',
    f'DNA_rep1',
    f'DNA_rep2',
    f'DNA_rep3'
    ]
    
    df = df[new_order]
    
    return df

def main():
#    if not os.path.exists(data_df_path):
        barcode_coverage, data_df = loop_through_bams()
        print(barcode_coverage, 'barcode coverage per file')
#    else:
#        data_df = pd.read_csv(f'{project_dir}/barcode_counts_UMI.txt', sep='\t', index_col = 'oligo_bc')
main()


