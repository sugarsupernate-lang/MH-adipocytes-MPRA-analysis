

# three paths are output directory, MPRAnalyze array, and split index file
csplit -f /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/rna_chunks_filter_sorted/rna -s -n 4 /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/adipo_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt $(cat /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/split_positions.txt)

csplit -f /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/dna_chunks_filter_sorted/dna -s -n 4 /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/adipo_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt $(cat /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/split_positions.txt)

