library("edgeR")
setwd("/Users/jing/Desktop/bcbc/Projects/Basil_HiFi/trichomes_RNAseq/Normalization")
tomct <- read.csv("gene_count_matrix.csv", row.names=1, header=TRUE)
my_group <- (c("Sweet_basil_root","Sweet_basil_trichome",	"Sweet_basil_leaf",	"Sweet_basil_leaf_trichome"))
tom_dgl <- DGEList(counts=tomct, group=my_group)
tom_dgl
head(tom_dgl$counts)
dim(tom_dgl$counts)
tom_dgl$samples
tom_dgl.rawdata <- tom_dgl
keep <- rowSums(cpm(tom_dgl)>1) >= 2
table(keep)
tom_dgl <- tom_dgl[keep, , keep.lib.sizes=FALSE]
tom_dgl <- calcNormFactors(tom_dgl)
cpm <- cpm(tom_dgl, log = FALSE, normalized.lib.sizes=TRUE)
write.csv(cpm,file="trichome_RNAseq_CPM_Adrian.csv" )
logcounts <- cpm(tom_dgl,log=TRUE)
write.csv(logcounts, file="log_normalized_counts_trichome_Adrian.csv")


