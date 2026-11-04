#install.packages('devtools')

# Load packages

#library("devtools")

# install the correct version of MPRAnalyze
# devtools::install_github("YosefLab/MPRAnalyze", force=TRUE)

library("MPRAnalyze")
library("glue")

sessionInfo()
#packageVersion("MPRAnalyze")

args = commandArgs(trailingOnly=TRUE)
print(args)
lib <- args[1]
adaptor <- args[2]
cells <- args[3]



# Read in data

annot = read.csv(
	glue("/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_quantitative/adipo_col_annot_quantitative.txt"), 
	sep = '\t', 
	row.names = 1
	)

annot$replicate<-as.factor(annot$replicate)
annot$barcode<-as.factor(annot$barcode)

dcounts_table = read.table(
	glue("/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_quantitative/adipo_DNA_MPRAnalyze_quantitative.txt"), 
	header=T, 
	sep='\t', 
	stringsAsFactors=F, 
	quote="", 
	row.names=1
	)
rcounts_table = read.table(
	glue("/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_quantitative/adipo_RNA_MPRAnalyze_quantitative.txt"), 
	header=T, 
	sep='\t', 
	stringsAsFactors=F, 
	quote="", 
	row.names=1
	)

# Coerce to matrices for later

dcounts = as.matrix(as.data.frame(dcounts_table))
rcounts = as.matrix(as.data.frame(rcounts_table))
neg_ctrls = grepl("scrambled", row.names(dcounts))   # could be done on the dcounts or the rcounts, they're the same
non_neg_ctrls = !grepl("scrambled", row.names(dcounts))  # all the rows that aren't negative controls

# Estimate and extract depth factors on entire dataset

obj <- MpraObject(dnaCounts = dcounts, rnaCounts = rcounts, colAnnot = annot, controls=neg_ctrls)
obj <- estimateDepthFactors(obj, lib.factor = "replicate", which.lib="both")
obj <- analyzeQuantification(
	obj,
	rnaDesign = ~1,
	dnaDesign = ~ replicate
	)
alpha_res <- getAlpha(obj)
res <- testEmpirical(obj, useControls=TRUE)

# Add alphas to res
identical(rownames(alpha_res), rownames(res)) # should be TRUE
res$alpha <- cbind(alpha_res$alpha)

# save(list = c("res", "alpha_res","obj"), file = "data") #KL 30.1.23 - not sure why we this is done. not really needed

# added by SF: this line was previously run before the mpranalyze model (on dcounts), before running mpranalyze, which ommits several enahncers. Now moved here and run on res
non_neg_ctrls = !grepl("scrambled", row.names(res))  # all the rows that aren't negative controls

# Fix FDRs to correct values excluding the negative controls from that calculation (we don't need to include them in the tests)
# RE 2026-08-19 (adipo/MH port): scrambled controls are now EXCLUDED from the BH
# correction, applied inline here rather than as a post-hoc patch. This is the same
# correction that MH npc applies after the fact via scripts_rerun/fix_quant_fdr_posthoc.R
# and that AH applies inside its rerun prep job. Including scrambled inflates the
# null and deflates the number of oligos called active.
res$fdr.mad <- NA_real_
res$fdr.mad[non_neg_ctrls] <- p.adjust(res$pval.mad[non_neg_ctrls], method='BH')

# Write results to out

write.table(
			alpha_res, 
            file=glue("/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_quantitative/adipo_alpha_full_quantitative_nobc_annot_factor_new_design.txt"), 
            sep="\t", 
            quote = F, 
            na = "", 
            row.names = T, 
            col.names = T
            )

write.table(
			res, 
            file=glue("/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_quantitative/adipo_results_full_quantitative_nobc_fdr_annot_factor_new_design.txt"), 
            sep="\t", 
            quote = F, 
            na = "", 
            row.names = T, 
            col.names = T
            )

