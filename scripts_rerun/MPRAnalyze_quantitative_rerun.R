#!/usr/bin/env Rscript
# MH NPC — MPRAnalyze_quantitative.R FIXED VERSION.
#
# Bug in original (line 83):
#   res$fdr.mad <- p.adjust(res$pval.mad, method='BH')
# included scrambled controls in the BH correction, inflating FDRs for real
# oligos. The correct line was left commented out:
#   # res$fdr.mad[non_neg_ctrls] <- p.adjust(res$pval.mad[non_neg_ctrls], method='BH')
#
# This rerun script: same model, but FDR computed on non-scrambled oligos only.
# Output file has _rerun suffix so the original is preserved.

suppressPackageStartupMessages({
  library(MPRAnalyze)
  library(glue)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript MPRAnalyze_quantitative_rerun.R <annot> <dna> <rna> <alpha_out> <results_out>")
}
annot_file  <- args[1]
dna_file    <- args[2]
rna_file    <- args[3]
alpha_out   <- args[4]
results_out <- args[5]

annot <- read.csv(annot_file, sep = "\t", row.names = 1)
annot$replicate <- as.factor(annot$replicate)
annot$barcode   <- as.factor(annot$barcode)

dcounts_table <- read.table(dna_file, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, quote = "", row.names = 1)
rcounts_table <- read.table(rna_file, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, quote = "", row.names = 1)

dcounts <- as.matrix(dcounts_table)
rcounts <- as.matrix(rcounts_table)

neg_ctrls    <- grepl("scrambled", rownames(dcounts))
cat("Scrambled oligos:", sum(neg_ctrls), "\n")
cat("Non-scrambled:",    sum(!neg_ctrls), "\n")

obj <- MpraObject(dnaCounts = dcounts, rnaCounts = rcounts,
                  colAnnot = annot, controls = neg_ctrls)
obj <- estimateDepthFactors(obj, lib.factor = "replicate", which.lib = "both")
obj <- analyzeQuantification(obj, rnaDesign = ~ 1, dnaDesign = ~ replicate)

alpha_res <- getAlpha(obj)
res       <- testEmpirical(obj, useControls = TRUE)
stopifnot(identical(rownames(alpha_res), rownames(res)))
res$alpha <- alpha_res$alpha

# >>> FIX: exclude scrambled from FDR calculation <<<
non_neg_ctrls <- !grepl("scrambled", rownames(res))
res$fdr.mad <- NA_real_
res$fdr.mad[non_neg_ctrls] <- p.adjust(res$pval.mad[non_neg_ctrls], method = "BH")

cat("Old (all-oligo) BH result would have been:\n")
cat("  Active @ FDR<0.05:", sum(p.adjust(res$pval.mad, method = "BH") < 0.05, na.rm = TRUE), "\n")
cat("NEW (non-scrambled) BH result:\n")
cat("  Active @ FDR<0.05:", sum(res$fdr.mad < 0.05, na.rm = TRUE), "\n")

write.table(alpha_res, file = alpha_out, sep = "\t", quote = FALSE,
            na = "", row.names = TRUE, col.names = TRUE)
write.table(res,       file = results_out, sep = "\t", quote = FALSE,
            na = "", row.names = TRUE, col.names = TRUE)
