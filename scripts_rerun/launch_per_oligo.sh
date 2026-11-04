#!/bin/bash
# =============================================================================
# MH HOb per-oligo rerun — one SGE array task per oligo.
#
# Prerequisites: launch_pipeline.sh has already run through step 1d (sort),
# producing:
#   $COMP_RERUN/${CELL_PREFIX}_{RNA,DNA}_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt
#   $COMP_RERUN/${CELL_PREFIX}_col_annot_comparative_improved.txt
#
# This launcher:
#   1. Pre-computes depth factors once (fast, reads the full matrix)
#   2. Generates an array job with one task per oligo
#   3. Submits finalize held on array completion
# =============================================================================

set -e
source /etc/profile 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true

CELL=adipo
CELL_PREFIX=adipo
BASE=/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo
OUTPUT_SUB=output
SCRIPTS_RERUN=$BASE/scripts_rerun
COMP_RERUN=$BASE/$OUTPUT_SUB/mpranalyze_comparative_rerun

# For MH HOb the comp matrix prefix is HOb; for BJ it's also HOb (see BJ patch)
# Override below if running NPC (uses adipo_)
COMP_FILE_PREFIX=adipo

ANNOT=$COMP_RERUN/${CELL_PREFIX}_col_annot_comparative_improved.txt
DNA_FULL=$COMP_RERUN/${COMP_FILE_PREFIX}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt
RNA_FULL=$COMP_RERUN/${COMP_FILE_PREFIX}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt
DEPTH_RDS=$COMP_RERUN/depth_factors.rds

PER_OLIGO_DIR=$COMP_RERUN/per_oligo
RESULTS_DIR=$PER_OLIGO_DIR/results
LOG_DIR=$PER_OLIGO_DIR/logs
JOB_DIR=$PER_OLIGO_DIR/jobs

mkdir -p $RESULTS_DIR $LOG_DIR $JOB_DIR

N_OLIGOS=$(wc -l < $DNA_FULL)
echo "=== MH $CELL_PREFIX per-oligo rerun ==="
echo "N oligos: $N_OLIGOS"

# ─────────────────────────────────────────────────────────────────────────────
# 1) Pre-compute depth factors (one job)
# ─────────────────────────────────────────────────────────────────────────────
cat > $JOB_DIR/job_precompute.sh << HEREDOC_END
#!/bin/bash
#\$ -S /bin/bash
#\$ -cwd
#\$ -j y
#\$ -N MH_${CELL_PREFIX}_precompute
#\$ -l h_rt=02:00:00
#\$ -l mem_free=32G
#\$ -o $LOG_DIR/precompute.log

source /etc/profile 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda deactivate 2>/dev/null || true
conda deactivate 2>/dev/null || true
module purge
module load CBI r/4.3.2

if [ -f $DEPTH_RDS ]; then echo "Depth factors already exist, skipping"; else
Rscript $SCRIPTS_RERUN/precompute_depth_factors.R $ANNOT $DNA_FULL $RNA_FULL $DEPTH_RDS
fi
HEREDOC_END

PRECOMP_ID=$(qsub -terse $JOB_DIR/job_precompute.sh)
echo "Precompute: $PRECOMP_ID"

# ─────────────────────────────────────────────────────────────────────────────
# 2) Per-oligo array job (held on precompute)
# ─────────────────────────────────────────────────────────────────────────────
cat > $JOB_DIR/array_per_oligo.sh << HEREDOC_END
#!/bin/bash
#\$ -S /bin/bash
#\$ -N MH_${CELL_PREFIX}_oligo
#\$ -o $LOG_DIR
#\$ -e $LOG_DIR
#\$ -j y
#\$ -l mem_free=4G
#\$ -l h_rt=04:00:00
#\$ -r y
#\$ -t 1-$N_OLIGOS

# Skip if result already exists
OUT=$RESULTS_DIR/oligo_\${SGE_TASK_ID}.txt
if [ -f "\$OUT" ]; then
  echo "Already done: \$OUT"
  exit 0
fi

source /etc/profile 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda deactivate 2>/dev/null || true
conda deactivate 2>/dev/null || true
module purge
module load CBI r/4.3.2

Rscript $SCRIPTS_RERUN/MPRAnalyze_one_oligo.R $DEPTH_RDS $DNA_FULL $RNA_FULL \$SGE_TASK_ID \$OUT
HEREDOC_END

ARRAY_ID=$(qsub -terse -hold_jid $PRECOMP_ID $JOB_DIR/array_per_oligo.sh | head -1 | cut -d. -f1)
echo "Per-oligo array: $ARRAY_ID ($N_OLIGOS tasks)"

# ─────────────────────────────────────────────────────────────────────────────
# 3) Finalize (concat + BH + sign)
# ─────────────────────────────────────────────────────────────────────────────
cat > $JOB_DIR/job_finalize.sh << HEREDOC_END
#!/bin/bash
#\$ -S /bin/bash
#\$ -N MH_${CELL_PREFIX}_finalize
#\$ -o $LOG_DIR/finalize.log
#\$ -j y
#\$ -l h_rt=04:00:00
#\$ -l mem_free=8G

source /etc/profile 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda deactivate 2>/dev/null || true
module purge
module load CBI r/4.3.2

cd $RESULTS_DIR
# Header from first existing result file (not hardcoded oligo_1)
FIRST_RES=\$(ls oligo_*.txt 2>/dev/null | head -1)
if [ -z "\$FIRST_RES" ]; then echo "ERROR: no result files"; exit 1; fi
head -1 \$FIRST_RES > all_results_combined.txt
MISSING=0
for i in \$(seq 1 $N_OLIGOS); do
  f=oligo_\$i.txt
  if [ -f "\$f" ]; then tail -n +2 "\$f" >> all_results_combined.txt; else MISSING=\$((MISSING+1)); fi
done
echo "Missing results: \$MISSING / $N_OLIGOS"
echo "Concat rows: \$(wc -l < all_results_combined.txt)"

Rscript -e '
res <- read.table("all_results_combined.txt", header=TRUE, sep="\t", row.names=1)
cat("Total:", nrow(res), "\n")
n_degen <- sum(is.na(res\$logFC))
cat("NA logFC dropped:", n_degen, "\n")
res <- res[!is.na(res\$logFC), ]
res\$fdr <- p.adjust(res\$pval, method="BH")
cat("Sig @ FDR<0.05:", sum(res\$fdr < 0.05), "\n")
write.table(res, "${CELL}_comparative_rerun_final.txt",
            sep="\t", quote=FALSE, col.names=NA)
'
HEREDOC_END

FINAL_ID=$(qsub -terse -hold_jid $ARRAY_ID $JOB_DIR/job_finalize.sh | head -1)
echo "Finalize: $FINAL_ID"
echo ""
echo "Output: $RESULTS_DIR/${CELL}_comparative_rerun_final.txt"
