#!/usr/bin/env Rscript
# Post-hoc FDR fix for existing MH HOb quant results.
# No need to re-run MPRAnalyze — just re-apply BH with scrambled excluded.
#
# Reads the original results file, preserves every column except fdr.mad,
# re-computes fdr.mad from pval.mad with scrambled excluded, and writes
# a new _rerun.txt file. Original is not modified.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript fix_quant_fdr_posthoc.R <orig_results.txt> <out_corrected.txt>")
}
orig_file <- args[1]
out_file  <- args[2]

res <- read.table(orig_file, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "", row.names = 1)
cat("Rows:", nrow(res), "\n")

old_fdr_sig  <- sum(res$fdr.mad < 0.05, na.rm = TRUE)
non_neg_ctrls <- !grepl("scrambled", rownames(res))
cat("Scrambled:", sum(!non_neg_ctrls), " Non-scrambled:", sum(non_neg_ctrls), "\n")

res$fdr.mad <- NA_real_
res$fdr.mad[non_neg_ctrls] <- p.adjust(res$pval.mad[non_neg_ctrls], method = "BH")

new_fdr_sig <- sum(res$fdr.mad < 0.05, na.rm = TRUE)
cat("Active (old FDR < 0.05):", old_fdr_sig, "\n")
cat("Active (new FDR < 0.05):", new_fdr_sig, "\n")

# NOTE: col.names=TRUE (not NA) to match the original MPRAnalyze_quantitative.R
# output format, which the downstream filter_mpra_comp_input_for_active.py
# relies on (pandas auto-detects the misaligned header and uses oligo names as index).
write.table(res, file = out_file, sep = "\t", quote = FALSE,
            na = "", row.names = TRUE, col.names = TRUE)
cat("Wrote", out_file, "\n")
