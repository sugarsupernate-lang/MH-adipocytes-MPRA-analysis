#!/bin/bash
# =============================================================================
# MH adipo MPRAnalyze rerun — fire and forget
#
# Fixes applied:
#   1. sort_comp_input_barcodes.py RNA/DNA filename swap (was writing RNA to
#      DNA-named file and vice versa → MPRAnalyze received inverted channels)
#   2. MPRAnalyze_quantitative.R FDR included scrambled controls (inflated FDR;
#      now uses non-scrambled subset for BH)
#   3. NA-logFC dropped before BH in comp post-processing
#
# All output → mpranalyze_*_rerun/ dirs. Originals untouched.
#
# Invoke: bash scripts_rerun/launch_pipeline.sh
# =============================================================================

set -e

# Ensure SGE + modules are available for ssh-command shells
source /etc/profile 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true

CELL=adipo
CELL_PREFIX=adipo
BASE=/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo
OUTPUT_SUB=output
SCRIPTS=$BASE/scripts
SCRIPTS_RERUN=$BASE/scripts_rerun

QUANT_ORIG=$BASE/$OUTPUT_SUB/mpranalyze_quantitative
QUANT_RERUN=$BASE/$OUTPUT_SUB/mpranalyze_quantitative_rerun
COMP_ORIG=$BASE/$OUTPUT_SUB/mpranalyze_comparative
COMP_RERUN=$BASE/$OUTPUT_SUB/mpranalyze_comparative_rerun

# Library-wide depth factors — same contract as launch_chunked_libdepth.sh.
# generate_chunk_jobs.py points every chunk at this RDS, and
# MPRAnalyze_comp_chunked.R injects it via setDepthFactors(), so depth must be
# estimated ONCE on the full filtered matrix rather than per chunk.
ANNOT=$COMP_RERUN/${CELL_PREFIX}_col_annot_comparative_improved.txt
DNA_FULL=$COMP_RERUN/${CELL_PREFIX}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt
RNA_FULL=$COMP_RERUN/${CELL_PREFIX}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt
DEPTH_RDS=$COMP_RERUN/depth_factors.rds

mkdir -p $QUANT_RERUN $COMP_RERUN/{rna_chunks,dna_chunks,chunk_jobs,chunk_logs,results}

echo "=== MH adipo Pipeline Rerun ==="
echo "Quant: $QUANT_RERUN"
echo "Comp:  $COMP_RERUN"

# ─────────────────────────────────────────────────────────────────────────────
# JOB 1: Prep (quant FDR fix + format + filter + sort + csplit + generate array)
# ─────────────────────────────────────────────────────────────────────────────
cat > $COMP_RERUN/job_prep.sh << HEREDOC_END
#!/bin/bash
#\$ -S /bin/bash
#\$ -cwd
#\$ -j y
#\$ -N MH_adipo_prep
#\$ -q long.q
#\$ -l h_rt=04:00:00
#\$ -l mem_free=32G
#\$ -o $COMP_RERUN/prep.log
#\$ -wd $BASE

set -e
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
echo "=== STEP 1a: Quant FDR post-hoc fix (scrambled excluded) ==="
conda deactivate 2>/dev/null || true
module load CBI r/4.3.2

ORIG=$QUANT_ORIG/${CELL}_results_full_quantitative_nobc_fdr_annot_factor_new_design.txt
CORRECTED=$QUANT_RERUN/${CELL}_quant_corrected.txt

Rscript $SCRIPTS_RERUN/fix_quant_fdr_posthoc.R \$ORIG \$CORRECTED
module unload r 2>/dev/null || true
echo "=== STEP 1a done ==="

echo "=== STEP 1b: Re-format comparative input ==="
conda activate mpra_association_seq
cd $SCRIPTS
python3 -u format_for_MPRAnalyze_comparative_analysis_final.py
cp $COMP_ORIG/${CELL_PREFIX}_RNA_MPRAnalyze_comparative_improved.txt $COMP_RERUN/
cp $COMP_ORIG/${CELL_PREFIX}_DNA_MPRAnalyze_comparative_improved.txt $COMP_RERUN/
cp $COMP_ORIG/${CELL_PREFIX}_col_annot_comparative_improved.txt $COMP_RERUN/
echo "=== STEP 1b done ==="

echo "=== STEP 1c: Filter for active (uses corrected quant FDR) ==="
QUANT_ORIG_FILE=$QUANT_ORIG/${CELL}_results_full_quantitative_nobc_fdr_annot_factor_new_design.txt
if [ ! -f \$QUANT_ORIG_FILE.bak_rerun ]; then
    cp \$QUANT_ORIG_FILE \$QUANT_ORIG_FILE.bak_rerun
fi
cp \$CORRECTED \$QUANT_ORIG_FILE
python3 -u filter_mpra_comp_input_for_active.py
cp \$QUANT_ORIG_FILE.bak_rerun \$QUANT_ORIG_FILE
cp $COMP_ORIG/${CELL_PREFIX}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt $COMP_RERUN/
cp $COMP_ORIG/${CELL_PREFIX}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt $COMP_RERUN/
echo "=== STEP 1c done ==="

echo "=== STEP 1d: Sort with FIXED filenames ==="
python3 -u $SCRIPTS/sort_comp_input_barcodes.py
echo "=== STEP 1d done ==="

echo "=== STEP 1e: csplit into chunks ==="
csplit -f $COMP_RERUN/rna_chunks/rna -s -n 4 \\
    $COMP_RERUN/${CELL_PREFIX}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt \\
    \$(cat $COMP_RERUN/split_positions.txt)
csplit -f $COMP_RERUN/dna_chunks/dna -s -n 4 \\
    $COMP_RERUN/${CELL_PREFIX}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt \\
    \$(cat $COMP_RERUN/split_positions.txt)
N_CHUNKS=\$(ls $COMP_RERUN/rna_chunks/rna* | wc -l)
echo "Created \$N_CHUNKS chunks"

echo "=== STEP 1f: Precompute library-wide depth factors ==="
conda deactivate 2>/dev/null || true
module purge 2>/dev/null || true
module load CBI r/4.3.2
if [ -f $DEPTH_RDS ]; then
  echo "Depth factors already exist at $DEPTH_RDS — skipping"
else
  Rscript $SCRIPTS_RERUN/precompute_depth_factors.R $ANNOT $DNA_FULL $RNA_FULL $DEPTH_RDS
fi
conda activate mpra_association_seq

echo "=== STEP 1g: Generate array job ==="
python3 -u $SCRIPTS_RERUN/generate_chunk_jobs.py
echo "=== PREP DONE ==="
HEREDOC_END

PREP_ID=$(qsub -terse $COMP_RERUN/job_prep.sh)
echo "Prep job: $PREP_ID"

# ─────────────────────────────────────────────────────────────────────────────
# JOB 2: Comparative array (held until prep finishes)
# ─────────────────────────────────────────────────────────────────────────────
cat > $COMP_RERUN/job_launch_array.sh << HEREDOC_END
#!/bin/bash
#\$ -S /bin/bash
#\$ -cwd
#\$ -j y
#\$ -N MH_adipo_launch
#\$ -l h_rt=00:30:00
#\$ -l mem_free=4G
#\$ -o $COMP_RERUN/launch.log
#\$ -wd $BASE

cd $COMP_RERUN/chunk_jobs
# Capture the array's REAL job id. See the hold-ordering note below -- the
# finalize job cannot depend on this array by NAME, because it was submitted
# before this array existed.
ARRAY_ID=\$(qsub -terse array_job.sh | cut -d. -f1)
echo "chunked array submitted: \$ARRAY_ID"
FIN_ID=\$(cat $COMP_RERUN/finalize.jobid 2>/dev/null || true)
if [ -n "\$FIN_ID" ]; then
  qalter -hold_jid "\$ARRAY_ID" "\$FIN_ID" && \
    echo "repointed finalize \$FIN_ID to depend on array \$ARRAY_ID"
else
  echo "WARNING: no finalize.jobid found -- finalize may run before the array finishes"
fi
HEREDOC_END

LAUNCH_ID=$(qsub -terse -hold_jid $PREP_ID $COMP_RERUN/job_launch_array.sh)
echo "Array launcher: $LAUNCH_ID (held on $PREP_ID)"

# ─────────────────────────────────────────────────────────────────────────────
# JOB 3: Finalize (concat chunks, drop NA-logFC, BH, sign-flip handled by R)
# ─────────────────────────────────────────────────────────────────────────────
cat > $COMP_RERUN/job_finalize.sh << HEREDOC_END
#!/bin/bash
#\$ -S /bin/bash
#\$ -cwd
#\$ -j y
#\$ -N MH_adipo_finalize
#\$ -q long.q
#\$ -l h_rt=02:00:00
#\$ -l mem_free=16G
#\$ -o $COMP_RERUN/finalize.log
#\$ -wd $BASE

set -e
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
cd $COMP_RERUN/results

# Concat with header preserved
head -1 comp_result_0000.txt > all_results_combined.txt
for f in comp_result_*.txt; do tail -n +2 "\$f" >> all_results_combined.txt; done
echo "Concat done: \$(wc -l < all_results_combined.txt) rows"

conda deactivate 2>/dev/null || true
module load CBI r/4.3.2

Rscript -e '
res <- read.table("all_results_combined.txt", header=TRUE, sep="\t", row.names=1)
cat("Total:", nrow(res), "\n")
n_degen <- sum(is.na(res\$logFC))
cat("Degenerate (NA logFC) dropped:", n_degen, "\n")
res <- res[!is.na(res\$logFC), ]
cat("After drop:", nrow(res), "\n")
# NOTE: unlike AH, do NOT negate logFC here — MH factor default (archaic<modern
# alphabetical) gives logFC = log2(modern/archaic), which is the MH convention.
res\$fdr <- p.adjust(res\$pval, method="BH")
cat("Sig @ FDR<0.05:", sum(res\$fdr < 0.05), "\n")
write.table(res, "${CELL}_comparative_rerun_final.txt",
            sep="\t", quote=FALSE, col.names=NA)
'
echo "=== FINALIZE DONE ==="
HEREDOC_END

# ── HOLD ORDERING BUG, FIXED 2026-08-23 (RE) ────────────────────────────────
# Was: -hold_jid "$LAUNCH_ID,MH_adipo_libdepth_comp"
#
# SGE resolves a hold_jid NAME at SUBMIT time. MH_adipo_libdepth_comp is the
# chunked array, which does not exist yet at this point -- it is created later,
# by the launcher job. A hold on a nonexistent job name is not an error: SGE
# treats it as ALREADY SATISFIED. So finalize was held only on the launcher, and
# went eligible the moment the launcher exited -- while all 191 array tasks were
# still running. It would then have concatenated an EMPTY results/ directory and
# written a final table that looked finished.
#
# Observed live on 2026-08-23: finalize sat in `qw`, not `hqw`, with the array
# still queued. This is the same failure that produced BJ's 115-byte corrupt
# "canonical" output.
#
# Fix: submit finalize held only on the launcher, write its job id to disk, and
# have the launcher qalter it onto the array's REAL numeric job id once that
# array exists. Numeric ids cannot be silently unresolvable.
FINAL_ID=$(qsub -terse -hold_jid "$LAUNCH_ID" $COMP_RERUN/job_finalize.sh)
echo "$FINAL_ID" > $COMP_RERUN/finalize.jobid
echo "Finalize: $FINAL_ID (launcher will repoint it onto the array job id)"

echo ""
echo "Pipeline submitted. Monitor with: qstat -u $USER"
echo "  Prep:     $PREP_ID"
echo "  Launch:   $LAUNCH_ID"
echo "  Finalize: $FINAL_ID"
echo ""
echo "Output files (after completion):"
echo "  Quant corrected: $QUANT_RERUN/${CELL}_quant_corrected.txt"
echo "  Comp final:      $COMP_RERUN/results/${CELL}_comparative_rerun_final.txt"
