####Script for Lena for DE
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install("statmod")
BiocManager::install("corrplot")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("ShortRead")
library("edgeR")
library("corrplot")

  setwd("/Users/jing/Desktop/bcbc/Projects/Basil_HiFi/3primer_RNAseq/Normalization_comparison_500slop")

#Get counts with featureCount subread package, use as input
##Read in files, group, and format
y <- read.delim("ht_seq_500_count.csv", row.names=1, stringsAsFactors=FALSE, sep = ",")
raw_counts = y

#make object
samples <- read.delim("samples.csv", header=TRUE, sep =",")
group <- factor(samples$group)
cbind(samples, Group=group)
y <- DGEList(counts=y, group=samples$group)

#Filter
keep <- rowSums(cpm(y) > 1) >= 2  #keep genes with at least 5 cpm in at least 3 samples
y <- y[keep,]
table(keep)
dim(keep)

##Re-compute the library sizes:
y$samples$lib.size <- colSums(y$counts)

#check library sizes
barplot(y$samples$lib.size, main="Raw Counts", xlab = "sample", ylab = "counts", ylim = c(0,8000000))

# Check distributions of samples using boxplots
# Let's add a blue horizontal line that corresponds to the median logCPM
#color by group
#par(mfrow=c(1,1),oma=c(200,0,0,0))

#par(cex.lab=3, cex.axis=3, mar=c(30, 4, 5, 2))

par(cex.lab=0.7, cex.axis=0.7, mar = c(2, 1, 0.4, 0.5))

group.col <- c("goldenrod1","deepskyblue", "firebrick", "darkolivegreen2", "firebrick2", "green3", "orchid1", "wheat1", "mistyrose", "gold1", "lightsalmon2", "darkgreen", "darkseagreen2", "darkgoldenrod", "lightgoldenrod1", "magenta2", "brown2", "coral", "royalblue4", "lightsteelblue1", "yellow", "dodgerblue3")[group]
png(filename = 'log2_cpm_boxplot_500slop.png', width = 2200, height = 1800, units = 'px')
boxplot(cpm(y,log = TRUE), xlab="", ylab="Log2 counts per million",las=2,col=group.col,
        )
abline(h=median(cpm(y,log = TRUE),col="blue"))
title("Boxplots of logCPMs\n(coloured by groups)",cex.main=1.6)
dev.off()

#Estimate normalization factors, By default, calcNormFactors uses the TMM method and the sample whose 75%-ile (of library-scale-scaled counts) is closest to the mean of 75%-iles as the reference.
y = calcNormFactors(y)
y$samples

#calc log counts
logcounts <- cpm(y,log=TRUE)
write.csv(logcounts, file="log_normalized_counts_500slop.csv")

#examine log expression of a sample versus all others
par(mfrow=c(1,2))
plotMD(logcounts,column=1)
abline(h=0,col="grey")
plotMD(y,column = 2)
abline(h=0,col="grey")
par(mfrow=c(1,1))


#estimate dipsersion
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
plotBCV(y)

dgListGroups <- c(rep("Control",3),rep("Treat",3))

#y<- estimateDisp(y, design, robust=TRUE)
summary(y$prior.df)
sqrt(y$common.disp)  #The square root of the common dispersion gives the coefficient of variation of biological variation (BCV).
png(filename = 'bcv_500slop.png', width = 800, height = 800, units = 'px')
plotBCV(y)
dev.off()

########QC############
#multidimensional scaling plot
png(filename = 'mds_500slop.png', width = 2400, height = 2400, units = 'px')
plotMDS(y, method="bcv")
dev.off()

png(filename = 'dendrogram_500slop.png', width = 1600, height = 1600, units = 'px')
tcounts <- t(as.table(as.matrix(y)))
counts.dist = hclust(dist(tcounts))
plot(counts.dist)
dev.off()

#scatterplot of reps
#png(filename = 'pnp1_pnp2_scatter.png', width = 800, height = 800, units = 'px')
y.cpm <-cpm(y, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25)
#plot(y.cpm[,"MtGv1"], y.cpm[,"MtGv2"], log="xy")
#plot(y.cpm[,"E08_WTs_Gv_1"], y.cpm[,"E09_WTs_Gv_2"], log="xy")
#dev.off()

write.csv(y.cpm, file="norm_cpm_all_500slop.csv")

#Correlation Matrix
png(filename = 'corr_matrix_500slop.png', width = 2400, height = 2400, units = 'px')
corrplot(cor(y.cpm), method="square", cl.lim=c(0,1), tl.col="black", addgrid.col="black", is.corr=FALSE, addCoef.col="white", number.cex = 0.5)
dev.off()

#############DE#######################
#Set up design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- glmFit(y, design)
logFC <- predFC(y,design,prior.count=1,dispersion=0.05)
cor(logFC[,1:7])  #correlation matrix of the pooled samples

#Estimate parameters
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y,design)


#Set up comparisons
# Get all genes that are dr in cfm in comparisons and pnp comparisons and upreg in cfm and pnp comparisons. 
my.contrasts <- makeContrasts(
  p2_perrie - p2_cardinal,
  p2_PB - p2_CB,
  p2_PA- p2_CA,
  p2_PB - p2_PA,
  p2_CB - p2_CA,
  levels=design
)

contrast.names <- colnames(my.contrasts)
d <- data.frame(matrix(nrow=nrow(y$counts), ncol=0))


#loop through contrasts
for (i in c(1:length(contrast.names))){
  
  #####To find genes diff:
  prefix.lrt  <- glmLRT(fit, contrast=my.contrasts[,contrast.names[i]])
  
  # total number of genes significantly up-regulated or down-regulated at 5% FDR
  summary(prefix.dt <- decideTestsDGE(prefix.lrt, adjust.method="fdr", p.value=0.05))
  
  #We can pick out which genes are DE:
  prefix.isDE <- as.logical(prefix.dt)
  prefix.DEnames <- rownames(y)[prefix.isDE]
  length(prefix.DEnames) <- nrow(y$counts) #make prefix.DEnames the same number of rows as d
  d <- cbind(d, prefix.DEnames)
  
  # plot all the logFCs against average count size, highlighting the DE genes
  plotSmear(prefix.lrt, de.tags=prefix.DEnames)
  abline(h=c(-1,1), col="blue") #The blue lines indicate 2-fold up or down.
  
  #export results
  outfile <- paste(contrast.names[i], ".csv", sep="")
  write.csv(prefix.lrt$table, file=outfile)
  #outfile2 <- paste(contrast.names[i], "sig.csv", sep="")
  #write.csv(prefix.DEnames, file=outfile2)
}

