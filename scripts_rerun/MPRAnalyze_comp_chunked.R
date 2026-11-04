#!/usr/bin/env Rscript
# Chunked MPRAnalyze comparative with LIBRARY-WIDE pre-computed depth factors.
# Matches the Weiss 2021 design: depth estimated once on the full filtered matrix,
# then injected into each chunk via setDepthFactors().
#
# Usage:
#   Rscript MPRAnalyze_comp_chunked.R <depth_rds> <dna_chunk> <rna_chunk> <result_out>

suppressPackageStartupMessages({ library(MPRAnalyze); library(BiocParallel) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript MPRAnalyze_comp_chunked.R <depth_rds> <dna_chunk> <rna_chunk> <result_out>")
}
depth_rds  <- args[1]
dna_file   <- args[2]
rna_file   <- args[3]
result_out <- args[4]

cat("Depth RDS:", depth_rds, "\nDNA chunk:", dna_file, "\nRNA chunk:", rna_file, "\n")

dep <- readRDS(depth_rds)

dna <- as.matrix(read.table(dna_file, header = FALSE, sep = "\t",
                            row.names = 1, stringsAsFactors = FALSE, quote = ""))
rna <- as.matrix(read.table(rna_file, header = FALSE, sep = "\t",
                            row.names = 1, stringsAsFactors = FALSE, quote = ""))
cat("Chunk dims: DNA", dim(dna), " RNA", dim(rna), "\n")

colnames(dna) <- dep$col_names
colnames(rna) <- dep$col_names

obj <- MpraObject(dnaCounts = dna, rnaCounts = rna, colAnnot = dep$colAnnot)
# Use pre-computed library-wide depth factors (not per-chunk)
obj <- setDepthFactors(obj, dnaDepth = dep$dnaDepth, rnaDepth = dep$rnaDepth)

# BiocParallel's registered default is MulticoreParam with as many workers as the
# node has cores (125 observed on Wynton). analyzeComparative defaults to
# BPPARAM = NULL, so without this line it forks that many workers inside a
# one-slot cgroup, which deadlocks the task at zero CPU. Observed 2026-08-24:
# 8 of 191 array tasks stalled this way; all 8 completed once SerialParam was set.
obj <- analyzeComparative(obj,
                          BPPARAM       = BiocParallel::SerialParam(),
                          dnaDesign     = ~ replicate + barcode_allele,
                          rnaDesign     = ~ allele,
                          correctControls = FALSE,
                          reducedDesign = ~ 1)

res <- testLrt(obj)
write.table(res, file = result_out, sep = "\t", quote = FALSE,
            na = "", row.names = TRUE, col.names = TRUE)
cat("Wrote", result_out, "(", nrow(res), "rows )\n")
