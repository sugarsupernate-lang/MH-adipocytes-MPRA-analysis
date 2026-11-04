#!/bin/bash
# =============================================================================
# Chunked MPRAnalyze comparative rerun with LIBRARY-WIDE depth factors.
# Replaces the per-oligo architecture.
#
# - Depth factors estimated ONCE on the full filtered matrix (Weiss 2021 design)
# - csplit into chunks (~50 oligos each)
# - 48h per-task wall time
# - Resumable (skips chunks with existing comp_result_*.txt)
#
# Assumes prior pipeline has already run through sort step, producing:
#   $COMP_RERUN/${CELL_PREFIX}_{RNA,DNA}_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt
#   $COMP_RERUN/${CELL_PREFIX}_col_annot_comparative_improved.txt
#   $COMP_RERUN/split_positions.txt
# =============================================================================

set -e
source /etc/profile 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true

# ── cell-type config (EDIT PER PIPELINE) ──
DATASET=MH
CELL=adipo
CELL_PREFIX=adipo
COMP_FILE_PREFIX=adipo
BASE=/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo
OUTPUT_SUB=output
# ──────────────────────────────────────────

SCRIPTS_RERUN=$BASE/scripts_rerun
COMP_RERUN=$BASE/$OUTPUT_SUB/mpranalyze_comparative_rerun

ANNOT=$COMP_RERUN/${CELL_PREFIX}_col_annot_comparative_improved.txt
DNA_FULL=$COMP_RERUN/${COMP_FILE_PREFIX}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt
RNA_FULL=$COMP_RERUN/${COMP_FILE_PREFIX}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt
DEPTH_RDS=$COMP_RERUN/depth_factors.rds

RNA_CHUNK_DIR=$COMP_RERUN/rna_chunks
DNA_CHUNK_DIR=$COMP_RERUN/dna_chunks
JOB_DIR=$COMP_RERUN/chunk_jobs
LOG_DIR=$COMP_RERUN/chunk_logs
RESULTS=$COMP_RERUN/results

mkdir -p $RNA_CHUNK_DIR $DNA_CHUNK_DIR $JOB_DIR $LOG_DIR $RESULTS

echo "=== Chunked libdepth rerun — $CELL_PREFIX ==="

# ─────────────────────────────────────────────────────────────────────────────
# JOB 1: Prep (precompute depth + csplit + generate array)
# ─────────────────────────────────────────────────────────────────────────────
cat > $JOB_DIR/job_prep.sh << HEREDOC_END
#!/bin/bash
#\$ -S /bin/bash
#\$ -cwd
#\$ -j y
#\$ -N ${DATASET}_${CELL}_libdepth_prep
#\$ -l h_rt=02:00:00
#\$ -l mem_free=32G
#\$ -o $LOG_DIR/prep.log

set -e
source /etc/profile 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda deactivate 2>/dev/null || true
conda deactivate 2>/dev/null || true
module purge
module load CBI r/4.3.2

echo "=== Precompute depth factors (library-wide) ==="
if [ -f $DEPTH_RDS ]; then
  echo "Depth factors already exist at $DEPTH_RDS — skipping"
else
  Rscript $SCRIPTS_RERUN/precompute_depth_factors.R $ANNOT $DNA_FULL $RNA_FULL $DEPTH_RDS
fi

echo "=== csplit into chunks ==="
# Clear any stale chunk files
rm -f $RNA_CHUNK_DIR/rna* $DNA_CHUNK_DIR/dna*
csplit -f $RNA_CHUNK_DIR/rna -s -n 4 $RNA_FULL \$(cat $COMP_RERUN/split_positions.txt)
csplit -f $DNA_CHUNK_DIR/dna -s -n 4 $DNA_FULL \$(cat $COMP_RERUN/split_positions.txt)
N=\$(ls $RNA_CHUNK_DIR/rna* | wc -l)
echo "Created \$N chunks"

conda activate mpra_association_seq
python3 -u $SCRIPTS_RERUN/generate_chunk_jobs.py
HEREDOC_END

PREP_ID=$(qsub -terse $JOB_DIR/job_prep.sh)
echo "Prep: $PREP_ID"

# ─────────────────────────────────────────────────────────────────────────────
# JOB 2: Launch array (held on prep)
# ─────────────────────────────────────────────────────────────────────────────
cat > $JOB_DIR/job_launch_array.sh << HEREDOC_END
#!/bin/bash
#\$ -S /bin/bash
#\$ -cwd
#\$ -j y
#\$ -N ${DATASET}_${CELL}_libdepth_launch
#\$ -l h_rt=00:30:00
#\$ -l mem_free=4G
#\$ -o $LOG_DIR/launch.log

cd $JOB_DIR
ARRAY_ID=\$(qsub -terse array_job.sh | cut -d. -f1)
echo "chunked array submitted: \$ARRAY_ID"
FIN_ID=\$(cat $JOB_DIR/finalize.jobid 2>/dev/null || true)
if [ -n "\$FIN_ID" ]; then
  qalter -hold_jid "\$ARRAY_ID" "\$FIN_ID" && \\
    echo "repointed finalize \$FIN_ID onto array \$ARRAY_ID"
else
  echo "WARNING: no finalize.jobid -- finalize may run before the array finishes"
fi
HEREDOC_END

LAUNCH_ID=$(qsub -terse -hold_jid $PREP_ID $JOB_DIR/job_launch_array.sh)
echo "Array launcher: $LAUNCH_ID"

# ─────────────────────────────────────────────────────────────────────────────
# JOB 3: Finalize (held on array via name)
# ─────────────────────────────────────────────────────────────────────────────
cat > $JOB_DIR/job_finalize.sh << HEREDOC_END
#!/bin/bash
#\$ -S /bin/bash
#\$ -cwd
#\$ -j y
#\$ -N ${DATASET}_${CELL}_libdepth_finalize
#\$ -l h_rt=02:00:00
#\$ -l mem_free=16G
#\$ -o $LOG_DIR/finalize.log

set -e
source /etc/profile 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda deactivate 2>/dev/null || true
conda deactivate 2>/dev/null || true
module purge
module load CBI r/4.3.2

cd $RESULTS
FIRST=\$(ls comp_result_*.txt 2>/dev/null | head -1)
if [ -z "\$FIRST" ]; then echo "ERROR: no comp_result files"; exit 1; fi
head -1 \$FIRST > all_results_combined.txt
N_MISSING=0
for f in comp_result_*.txt; do
  tail -n +2 "\$f" >> all_results_combined.txt
done
N_DONE=\$(ls comp_result_*.txt 2>/dev/null | wc -l)
echo "Concat done: \$(wc -l < all_results_combined.txt) rows from \$N_DONE chunks"

Rscript -e '
res <- read.table("all_results_combined.txt", header=TRUE, sep="\t", row.names=1)
cat("Total:", nrow(res), "\n")
cat("NA logFC dropped:", sum(is.na(res\$logFC)), "\n")
res <- res[!is.na(res\$logFC), ]
res\$fdr <- p.adjust(res\$pval, method="BH")
cat("Sig @ FDR<0.05:", sum(res\$fdr < 0.05), "\n")
write.table(res, "${CELL}_comparative_rerun_libdepth_final.txt",
            sep="\t", quote=FALSE, col.names=NA)
'
HEREDOC_END

# SGE resolves a hold_jid NAME at SUBMIT time. The chunked array does not exist
# yet here -- it is created later by the launcher job -- so a hold on its NAME is
# silently treated as already satisfied and finalize can run against an empty
# results/ directory. This already fired once in bj_count (0-byte combined table
# dated a day BEFORE its own chunk results). Hold on the launcher only, record the
# finalize job id, and let the launcher qalter it onto the array's real numeric id.
FINAL_ID=$(qsub -terse -hold_jid "$LAUNCH_ID" $JOB_DIR/job_finalize.sh)
echo "$FINAL_ID" > $JOB_DIR/finalize.jobid
echo "Finalize: $FINAL_ID"

echo ""
echo "Monitor: qstat -u \$USER"
echo "Output: $RESULTS/${CELL}_comparative_rerun_libdepth_final.txt"
