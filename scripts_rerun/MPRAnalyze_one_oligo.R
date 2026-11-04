#!/usr/bin/env Rscript
# Per-oligo MPRAnalyze comparative: extracts one oligo from full matrices,
# injects pre-computed depth factors, runs testLrt, writes one-row result.
#
# Usage:
#   Rscript MPRAnalyze_one_oligo.R <depth_rds> <dna_full> <rna_full> <row_index> <out>
#
# row_index is 1-based (matches the order of oligos in dna_full / rna_full).

suppressPackageStartupMessages({ library(MPRAnalyze) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript MPRAnalyze_one_oligo.R <depth_rds> <dna_full> <rna_full> <row_index> <out>")
}
depth_rds <- args[1]
dna_file  <- args[2]
rna_file  <- args[3]
row_idx   <- as.integer(args[4])
out_file  <- args[5]

dep <- readRDS(depth_rds)

# Read full DNA/RNA; extract just the one oligo's row.
# Files are headerless; row 1 is oligo 1 etc.
dna_df <- read.table(dna_file, header = FALSE, sep = "\t",
                     stringsAsFactors = FALSE, quote = "", row.names = 1)
rna_df <- read.table(rna_file, header = FALSE, sep = "\t",
                     stringsAsFactors = FALSE, quote = "", row.names = 1)

oligo_name <- rownames(dna_df)[row_idx]
cat("Processing row", row_idx, "oligo:", oligo_name, "\n")
stopifnot(oligo_name == rownames(rna_df)[row_idx])

dna1 <- as.matrix(dna_df[row_idx, , drop = FALSE])
rna1 <- as.matrix(rna_df[row_idx, , drop = FALSE])
colnames(dna1) <- dep$col_names
colnames(rna1) <- dep$col_names

obj <- MpraObject(dnaCounts = dna1, rnaCounts = rna1, colAnnot = dep$colAnnot)
obj <- setDepthFactors(obj, dnaDepth = dep$dnaDepth, rnaDepth = dep$rnaDepth)

obj <- analyzeComparative(obj,
                          dnaDesign     = ~ replicate + barcode_allele,
                          rnaDesign     = ~ allele,
                          correctControls = FALSE,
                          reducedDesign = ~ 1)

res <- testLrt(obj)
res$oligo <- oligo_name

write.table(res, file = out_file, sep = "\t", quote = FALSE,
            na = "", row.names = TRUE, col.names = TRUE)
cat("Wrote", out_file, "\n")
