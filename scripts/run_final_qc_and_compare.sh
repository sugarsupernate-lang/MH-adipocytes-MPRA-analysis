#!/bin/bash
#$ -S /bin/bash
#$ -N adipo_final_qc
#$ -cwd
#$ -j y
#$ -o /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/qc/final_qc.log
#$ -q long.q
#$ -l mem_free=64G
#$ -l h_rt=12:00:00
# Full QC suite (all 26 activity sections incl. 2.5/2.7) plus the independent-pass
# comparison, run once the stage-4 comparative has finalised.
set -uo pipefail
BASE=/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo
set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate mpra_association_seq; set -u
export PYTHONUNBUFFERED=1
SUM=$BASE/output/qc/SUMMARY.txt
: > "$SUM"
FINAL=$BASE/output/mpranalyze_comparative_rerun/results/adipo_comparative_rerun_final.txt
NCHUNK=$(ls $BASE/output/mpranalyze_comparative_rerun/results/comp_result_*.txt 2>/dev/null | wc -l)
{
  echo "=== stage 4 integrity ==="
  echo "chunk results: $NCHUNK/191"
  if [ "$NCHUNK" -lt 191 ]; then echo "WARNING: incomplete array -- final table is PARTIAL"; fi
  if [ -s "$FINAL" ]; then
    newest=$(ls -t $BASE/output/mpranalyze_comparative_rerun/results/comp_result_*.txt | head -1)
    if [ "$FINAL" -nt "$newest" ]; then echo "timestamp OK: final postdates every chunk"; else echo "WARNING: final does NOT postdate all chunks"; fi
    echo "final rows: $(($(wc -l < "$FINAL") - 1))"
  else
    echo "ERROR: no final table"; exit 1
  fi
} >> "$SUM"

python3 -u $BASE/scripts/mpra_qc_suite.py && echo "QC_OK" >> "$SUM" || { echo "QC_FAILED" >> "$SUM"; }
{
  echo; echo "=== QC verdicts ==="
  grep -E "^  -> " $BASE/output/qc/qc_report.txt | sed 's/^  -> //'
  echo; echo "=== not evaluated ==="
  grep "NOT EVALUATED" $BASE/output/qc/qc_report.txt | sed 's/^  //'
  echo; echo "=== sections 2.5 / 2.7 ==="
  awk '/2\.5  Differential activity distribution/,/-> /' $BASE/output/qc/qc_report.txt
  awk '/2\.7  Volcano plot/,/-> /' $BASE/output/qc/qc_report.txt
} >> "$SUM" 2>&1

echo >> "$SUM"; echo "=== comparison to independent MPRAsnakeflow pass ===" >> "$SUM"
python3 -u $BASE/scripts/compare_to_josei.py >> "$SUM" 2>&1 || echo "COMPARE_FAILED" >> "$SUM"
echo "=== ALL DONE $(date) ===" >> "$SUM"
