#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r n
#$ -N adipo_mh_quant
#$ -q long.q
#$ -l mem_free=48G
#$ -l h_rt=24:00:00
#$ -l scratch=4G
#$ -o /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/quant.log
#$ -wd /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo
# =============================================================================
# adipo_mh_quant — QUANTITATIVE (activity) stage.  *** STAGED, DO NOT AUTO-FIRE ***
# Submit by hand only after allreps_Filter.txt and UMI_std2_filter.txt have been
# eyeballed.
#
#   1. format_for_MPRAnalyze_quantitative_analysis_final.py
#        -> output/mpranalyze_quantitative/adipo_{RNA,DNA}_MPRAnalyze_quantitative.txt
#        -> output/mpranalyze_quantitative/adipo_col_annot_quantitative.txt
#   2. MPRAnalyze_quantitative.R
#        -> output/mpranalyze_quantitative/adipo_results_full_quantitative_nobc_fdr_annot_factor_new_design.txt
#        -> output/mpranalyze_quantitative/adipo_alpha_full_quantitative_nobc_annot_factor_new_design.txt
#
# The R script's FDR now excludes scrambled controls from the BH correction
# inline, so scripts_rerun/fix_quant_fdr_posthoc.R is NOT needed on this branch.
# Its devtools/install_github calls are commented out — required, because Wynton
# compute nodes have no internet.
# =============================================================================
set -euo pipefail

BASE=/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo
SCRIPTS=$BASE/scripts
QUANT=$BASE/output/mpranalyze_quantitative
mkdir -p "$QUANT"

echo "=== adipo_mh_quant START $(date) host=$(hostname) ==="
test -s "$BASE/output/UMI/UMI_std2_filter.txt" || { echo "ERROR: UMI_std2_filter.txt missing — run the count chain first." >&2; exit 1; }

source ~/miniconda3/etc/profile.d/conda.sh
export PYTHONUNBUFFERED=1
cd "$SCRIPTS"

echo ">>> STEP 1: format_for_MPRAnalyze_quantitative_analysis_final.py $(date)"
conda activate mpra_association_seq
python3 -u format_for_MPRAnalyze_quantitative_analysis_final.py
test -s "$QUANT/adipo_RNA_MPRAnalyze_quantitative.txt" || { echo "ERROR: quant RNA table not produced" >&2; exit 1; }

echo ">>> STEP 2: MPRAnalyze_quantitative.R $(date)"
conda deactivate 2>/dev/null || true
module purge 2>/dev/null || true
module load CBI r/4.3.2
Rscript MPRAnalyze_quantitative.R adipo NA adipo
test -s "$QUANT/adipo_results_full_quantitative_nobc_fdr_annot_factor_new_design.txt" || { echo "ERROR: quant results not produced" >&2; exit 1; }

echo "=== adipo_mh_quant COMPLETE $(date) ==="
echo "  Results: $QUANT/adipo_results_full_quantitative_nobc_fdr_annot_factor_new_design.txt"
echo "  Next: bash scripts_rerun/launch_pipeline.sh  (comparative orchestrator)"
