#!/usr/bin/env Rscript
# Pre-compute MPRAnalyze depth factors from the FULL filtered library.
# Saves as RDS so per-oligo jobs can inject them via setDepthFactors().
#
# Usage:
#   Rscript precompute_depth_factors.R <annot> <dna_full> <rna_full> <out_rds>

suppressPackageStartupMessages({ library(MPRAnalyze) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript precompute_depth_factors.R <annot> <dna_full> <rna_full> <out_rds>")
}
annot_file <- args[1]
dna_file   <- args[2]
rna_file   <- args[3]
out_rds    <- args[4]

annot <- read.table(annot_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "", row.names = 1)
dna   <- as.matrix(read.table(dna_file, header = FALSE, sep = "\t",
                              stringsAsFactors = FALSE, quote = "", row.names = 1))
rna   <- as.matrix(read.table(rna_file, header = FALSE, sep = "\t",
                              stringsAsFactors = FALSE, quote = "", row.names = 1))

cat("Full library dims: DNA", dim(dna), " RNA", dim(rna), "\n")

# Column names on dna/rna come back as V1, V2, ... — need to align with annot row names
colnames(dna) <- rownames(annot)
colnames(rna) <- rownames(annot)

obj <- MpraObject(dnaCounts = dna, rnaCounts = rna, colAnnot = annot)
obj <- estimateDepthFactors(obj, lib.factor = "replicate")

saveRDS(list(dnaDepth = obj@dnaDepth, rnaDepth = obj@rnaDepth,
             colAnnot = annot, col_names = rownames(annot)),
        file = out_rds)
cat("Saved depth factors to", out_rds, "\n")
cat("  n DNA depth:", length(obj@dnaDepth), "\n")
cat("  n RNA depth:", length(obj@rnaDepth), "\n")
