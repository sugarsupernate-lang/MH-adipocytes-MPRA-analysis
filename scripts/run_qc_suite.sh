#!/bin/bash
#$ -S /bin/bash
#$ -N adipo_qc
#$ -cwd
#$ -j y
#$ -o /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/qc/qc_suite.log
#$ -q long.q
#$ -l mem_free=64G
#$ -l h_rt=12:00:00
#$ -r y
# Gokhman Lab MPRA QC book (chapters 1-2), evaluated on the MH adipocyte MPRA.
# https://gokhmanlaborganization.github.io/MPRA_QC_book/
set -uo pipefail
mkdir -p /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/qc
set +u; source ~/miniconda3/etc/profile.d/conda.sh; conda activate mpra_association_seq; set -u
export PYTHONUNBUFFERED=1
if python3 -u /wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/scripts/mpra_qc_suite.py; then
  echo "=== QC SUITE DONE $(date) ==="
else
  echo "=== QC SUITE FAILED (exit $?) $(date) ==="; exit 1
fi
