#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r n
#$ -N adipo_mh_count_chain
#$ -hold_jid adipo_mh_align
#$ -q long.q
#$ -l mem_free=64G
#$ -l h_rt=24:00:00
#$ -l scratch=10G
#$ -o /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/count_chain.log
#$ -wd /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo
# =============================================================================
# adipo_mh_count_chain — front half of the MH adipocyte MPRA pipeline.
# Held on adipo_mh_align; fires once all 6 BAMs land in aligned_reads/.
#
#   1. RNA_DNA_processing.py   count   -> output/UMI/barcode_counts_UMI*.txt
#   2. filter_UMI_outliers.py  std-2   -> output/UMI/UMI_std2_filter.txt
#   3. calculate_ratio_df.py   ratios  -> output/allreps_Filter.txt
#
# NOTE there is no explode step here. That exists only in AH_counts/adipo, where
# calculate_ratio_df.py had been retargeted to a divergent filename. The MH
# lineage reads UMI/UMI_std2_filter.txt directly — filter_UMI_outliers.py already
# emits one row per barcode.
#
# If step 1 proves too heavy for a single job (AH hit this at ~436M pairs/rep),
# split it per replicate the way AH_counts/adipo did with count_one_rep.py +
# merge_reps.py. The 0715 run is lighter (DNA 41-46M, RNA 94-135M reads/rep),
# so the monolithic job is expected to hold.
# =============================================================================
set -euo pipefail

BASE=/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo
SCRIPTS=$BASE/scripts
BAMDIR=$BASE/aligned_reads

echo "============================================"
echo "  adipo_mh_count_chain START $(date)  host=$(hostname)"
echo "============================================"

for b in DNA_aligned_rep1 DNA_aligned_rep2 DNA_aligned_rep3 \
         RNA_aligned_rep1 RNA_aligned_rep2 RNA_aligned_rep3; do
    if [ ! -s "$BAMDIR/$b.bam" ]; then
        echo "ERROR: missing or empty BAM $BAMDIR/$b.bam — aborting." >&2
        exit 1
    fi
done
echo "All 6 BAMs present."
ls -la "$BAMDIR"/*.bam

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mpra_association_seq
export PYTHONUNBUFFERED=1
cd "$SCRIPTS"

echo ""; echo ">>> STEP 1: RNA_DNA_processing.py (count) $(date)"
python3 -u RNA_DNA_processing.py
test -s "$BASE/output/UMI/barcode_counts_UMI.txt"      || { echo "ERROR: barcode_counts_UMI.txt not produced" >&2; exit 1; }
test -s "$BASE/output/UMI/barcode_counts_UMI_rep1.txt" || { echo "ERROR: per-rep counts not produced" >&2; exit 1; }

echo ""; echo ">>> STEP 2: filter_UMI_outliers.py (std-2) $(date)"
python3 -u filter_UMI_outliers.py
test -s "$BASE/output/UMI/UMI_std2_filter.txt"         || { echo "ERROR: UMI_std2_filter.txt not produced" >&2; exit 1; }

echo ""; echo ">>> STEP 3: calculate_ratio_df.py (ratios) $(date)"
python3 -u calculate_ratio_df.py
test -s "$BASE/output/allreps_Filter.txt"              || { echo "ERROR: allreps_Filter.txt not produced" >&2; exit 1; }

echo ""
echo "============================================"
echo "  adipo_mh_count_chain COMPLETE $(date)"
echo "    $BASE/output/UMI/barcode_counts_UMI.txt (+ per-rep)"
echo "    $BASE/output/UMI/UMI_std2_filter.txt"
echo "    $BASE/output/allreps_Filter.txt   <-- EYEBALL THIS before MPRAnalyze"
echo "============================================"
