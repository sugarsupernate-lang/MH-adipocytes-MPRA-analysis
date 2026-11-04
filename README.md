# MH adipocytes MPRA — analysis scripts

This repository is an archive of the scripts provided from the Wynton analysis tree for the MH adipocyte MPRA. It documents the activity and differential-activity methods for QC review. Scripts are included **without modifications to their original source**, and retain hard-coded Wynton paths. This is not a turnkey portable workflow.

## Source

The scripts were retrieved from `.../ryder/archaic_mpra/MH_counts/adipo/`. Some scripts document MH-specific corrections and reruns.

## Pipeline overview

1. **Counting and filtering** — `job_count_chain.sh` invokes `scripts/RNA_DNA_processing.py`, `scripts/filter_UMI_outliers.py`, and `scripts/calculate_ratio_df.py` in sequence. This starts from pre-aligned BAMs and produces barcode counts, an outlier-filtered count table, and per-oligo RNA/DNA ratios. The SGE script has a hold on an upstream alignment job that is **not included** in this archive.
2. **Activity / quantitative analysis** — `job_quant.sh` invokes `scripts/format_for_MPRAnalyze_quantitative_analysis_final.py` and `scripts/MPRAnalyze_quantitative.R`. The R script uses `analyzeQuantification(rnaDesign = ~1, dnaDesign = ~replicate)` and `testEmpirical(..., useControls = TRUE)`. The uploaded R script recalculates `fdr.mad` with BH correction over non-scrambled oligos; scrambled controls are excluded from the BH correction.
3. **Comparative / differential activity rerun** — `scripts_rerun/launch_pipeline.sh` orchestrates a post-hoc correction of quantitative FDR, comparative input formatting, activity filtering, barcode sorting, splitting, precomputation of library-wide depth factors, chunked testing, and final global BH adjustment. The R analysis is `scripts_rerun/MPRAnalyze_comp_chunked.R`, with `dnaDesign = ~ replicate + barcode_allele`, `rnaDesign = ~ allele`, `correctControls = FALSE`, `reducedDesign = ~1`. It uses `BiocParallel::SerialParam()` and `setDepthFactors()` with factors estimated from the full filtered library. Before final BH correction, records with missing `logFC` are removed. The script documents `logFC = log2(Derived / Ancestral)` and does not invert the sign at finalization.
4. **QC** — `scripts/mpra_qc_suite.py` and related scripts implement supplemental QC from the pipeline outputs. These are QC tools, not the core differential-activity model.

## Differential result checked on HPC

The user verified the existence of `output/mpranalyze_comparative_rerun/results/adipo_comparative_rerun_final.txt` (August 24, 2026), with **815 lines including header (814 result rows)**, and **191 `comp_result_*.txt` files**. The final table has `pval`, `fdr`, and `logFC` columns. These checks establish file existence and shape, **not** a complete provenance audit, execution-log validation, or equivalence to a separate analysis. Data files and logs are not included here.

## Important caveats for interpretation

- The original quantitative script and `scripts_rerun/fix_quant_fdr_posthoc.R` both describe excluding scrambled controls when adjusting activity FDR. The rerun launcher still performs a post-hoc correction, so the exact input/output provenance should be checked before describing which correction was operative for an individual reported table.
- `scripts/sort_comp_input_barcodes.py` explicitly documents fixing swapped RNA/DNA output filenames relative to an earlier version.
- `scripts_rerun/launch_pipeline.sh` documents changes to handling missing `logFC`, library-wide depth normalization and SGE scheduling.
- Alternate experimental/helper implementations are retained (e.g., `MPRAnalyze_comp.R`, `MPRAnalyze_comp_chunked_serial.R`, `MPRAnalyze_one_oligo.R`, and alternative launchers) but should **not** all be interpreted as scripts executed for the checked final file.
- Environment and dependency versions were not fully captured in this archive. The quantitative launch script references R 4.3.2 and an `mpra_association_seq` Conda environment, but package versions and upstream commits were not verified.

## Directory map

- `job_count_chain.sh`, `job_quant.sh`: entry points for count and quantitative stages.
- `scripts/`: original counting, input formatting, quantitative and comparative scripts, QC utilities.
- `scripts_rerun/`: rerun launcher, full-library depth-factor estimation, chunked comparative model, FDR correction and helper/alternative scripts.

## Reproduction

To rerun the workflow, obtain the original reference/association files, aligned sequencing data, intermediate output tables, original software environments and Wynton/SGE configuration. Replace the project-specific absolute paths as appropriate. No raw sequencing or result matrices are distributed in this code-only archive.
