#library(readxl)
library(MPRAnalyze)
library(glue)

args = commandArgs(trailingOnly=TRUE)

dna_file <- args[1]
rna_file <- args[2]
output_file <- args[3]

mpra_colannot <- read.table(
  glue("/wynton/group/ahituv/ryder/archaic_mpra/MH_counts/adipo/output/mpranalyze_comparative/adipo_col_annot_comparative_improved.txt"),
  header=T, sep='\t', stringsAsFactors=F, quote="", row.names=1)
mpra_data_DNA <- read.table(
  dna_file,
  header=F, sep='\t', stringsAsFactors=F, quote="", row.names=1)
mpra_data_RNA <- read.table(
  rna_file,
  header=F, sep='\t', stringsAsFactors=F, quote="", row.names=1)

print('done reading in data for full model')

# convert to matrix from data.frame
ce.dnaCounts <- as.matrix(as.data.frame(mpra_data_DNA))
ce.rnaCounts <- as.matrix(as.data.frame(mpra_data_RNA))

# assign column annotations, can stay as data.frame
ce.colAnnot <- mpra_colannot


### WITHOUT CONTROLS ###
# initialize MpraObject without controls
obj <- MpraObject(dnaCounts = ce.dnaCounts, rnaCounts = ce.rnaCounts, 
                 colAnnot = ce.colAnnot)
                  
print('made big MPRAobject')
# get depth factors
obj <- estimateDepthFactors(obj, lib.factor = "replicate")

d_depth <- obj@dnaDepth
r_depth <- obj@rnaDepth
print('calculated depth factors')

d = as.matrix(as.data.frame(mpra_data_DNA))
r = as.matrix(as.data.frame(mpra_data_RNA))

obj2 <- MpraObject(dnaCounts = d, rnaCounts = r, colAnnot = ce.colAnnot)
print('made MPRAobject')
obj2 <- setDepthFactors(obj2, dnaDepth = d_depth, rnaDepth = r_depth)

# 14.09.23 - changed the design so that replicates are not modelled in the RNA design
obj2 <- analyzeComparative(obj2, 
                           dnaDesign = ~ replicate+barcode_allele, 
                           rnaDesign = ~ allele,
                           correctControls = FALSE,
                           reducedDesign = ~ 1)
#(rnaDesign = ~ replicate + allele, dnaDesign = ~ replicate + barcode + allele, reducedDesign = ~ replicate)

### SF DESIGN ### 
#analyzeComparative(obj2,
#                  dnaDesign = ~ barcode_allele + replicate,
#                  rnaDesign = ~ allele,
#                  correctControls = FALSE,
#                  reducedDesign = ~1)



### CARLY'S DESIGN ###
#obj2 <- analyzeComparative(
#	obj2,
#	dnaDesign = ~ replicate + barcode_allele,
#	rnaDesign = ~ replicate + allele,
#	correctControls = FALSE,
#	reducedDesign = ~ replicate
#	)
print('analyze comparative done')

# test for differential activity
res <- testLrt(obj2)

# save results
write.table(res, file=output_file, 
            sep="\t", quote = F,na = "", row.names = T, col.names = T)
