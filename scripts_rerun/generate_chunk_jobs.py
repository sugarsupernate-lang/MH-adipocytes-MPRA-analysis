#!/usr/bin/env python
"""
Generate SGE array job for MPRAnalyze comparative chunks.
Uses LIBRARY-WIDE pre-computed depth factors (depth_factors.rds) injected per chunk.
48h per-task wall time.
"""
import os, glob, sys

COMP_DIR  = '/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative_rerun'
SCRIPTS   = '/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/scripts_rerun'
PREFIX    = 'adipo'
DATASET   = 'MH'
CELL      = 'adipo'  # lowercase for job-name matching with finalize hold_jid

rna_dir   = os.path.join(COMP_DIR, 'rna_chunks')
dna_dir   = os.path.join(COMP_DIR, 'dna_chunks')
jobs_dir  = os.path.join(COMP_DIR, 'chunk_jobs')
logs_dir  = os.path.join(COMP_DIR, 'chunk_logs')
results   = os.path.join(COMP_DIR, 'results')
depth_rds = os.path.join(COMP_DIR, 'depth_factors.rds')

for d in (jobs_dir, logs_dir, results):
    os.makedirs(d, exist_ok=True)

rna_chunks = sorted(glob.glob(os.path.join(rna_dir, 'rna*')))
if not rna_chunks:
    sys.exit(f"ERROR: no rna chunks in {rna_dir}")
n = len(rna_chunks)
print(f"{n} chunks")

with open(os.path.join(jobs_dir, 'chunk_manifest.tsv'), 'w') as f:
    f.write("task_id\tdna_chunk\trna_chunk\tresult_file\n")
    for i, rna_path in enumerate(rna_chunks):
        chunk_id = os.path.basename(rna_path).replace('rna', '')
        dna_path = os.path.join(dna_dir, f'dna{chunk_id}')
        result   = os.path.join(results, f'comp_result_{chunk_id}.txt')
        f.write(f"{i+1}\t{dna_path}\t{rna_path}\t{result}\n")

r_script = os.path.join(SCRIPTS, 'MPRAnalyze_comp_chunked.R')

array_script = os.path.join(jobs_dir, 'array_job.sh')
with open(array_script, 'w') as f:
    f.write(f"""#!/bin/bash
#$ -S /bin/bash
#$ -N {DATASET}_{CELL}_libdepth_comp
#$ -o {logs_dir}
#$ -e {logs_dir}
#$ -j y
#$ -l mem_free=8G
#$ -l h_rt=48:00:00
#$ -r y
#$ -t 1-{n}

MANIFEST={jobs_dir}/chunk_manifest.tsv
LINE=$(awk -v id=$SGE_TASK_ID 'NR==id+1' $MANIFEST)
DNA_CHUNK=$(echo "$LINE" | cut -f2)
RNA_CHUNK=$(echo "$LINE" | cut -f3)
RESULT=$(echo "$LINE"    | cut -f4)

if [ -f "$RESULT" ]; then echo "Already done: $RESULT"; exit 0; fi

source /etc/profile 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda deactivate 2>/dev/null || true
conda deactivate 2>/dev/null || true
module purge
module load CBI r/4.3.2
conda activate mpra_association_seq

echo "Task $SGE_TASK_ID: DNA=$DNA_CHUNK RNA=$RNA_CHUNK"
Rscript {r_script} {depth_rds} "$DNA_CHUNK" "$RNA_CHUNK" "$RESULT"
""")

print(f"Wrote {array_script} ({n} tasks)")
