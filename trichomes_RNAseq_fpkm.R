#change this to correct path
library("edgeR")
setwd("/Users/jing/Desktop/bcbc/Projects/Basil_HiFi/trichomes_RNAseq/FPKM/")

#Get counts with featureCount subread package, use as input
##Read in files, group, and format
y <- read.delim("gene_count_matrix.csv", row.names=1, sep = ",")
keep <- rowSums(cpm(y) > 1) >= 1

expressioncount <- read.delim("Perrie_cds_length.txt", row.names=1)
result <- rpkm(y, normalized.lib.sizes=FALSE, gene.length=expressioncount$length)
write.csv(result, file="trichomes_RNAseq_fpkm.csv")
