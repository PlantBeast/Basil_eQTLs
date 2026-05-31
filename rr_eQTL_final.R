#read parameters from command line

args <- commandArgs(trailingOnly=TRUE)

#load libraries 

library(qtl)
library(rlecuyer)
library(snow)
library(foreach)  
library(doSNOW)

#set command line paramters to named variables

input_file <- args[1]
output_file <- args[2]
number_cpus_cores <- as.numeric(args[3])

#verify named variables
print(input_file)
print(output_file)
print(typeof(number_cpus_cores))

# --- NEW: CREATE OUTPUT DIRECTORIES ---
dir.create("Phenotype exploration", showWarnings = FALSE)
dir.create("Genotype exploration", showWarnings = FALSE)
dir.create("13_LOD scores_scan1", showWarnings = FALSE)
dir.create("18_LOD scores_model part", showWarnings = FALSE)
dir.create("21_Marker effect plots", showWarnings = FALSE)
dir.create("QTL files", showWarnings = FALSE)
# --------------------------------------

# --- OPTIMIZATION: Setup Cluster Once Globally ---
cl <- makeCluster(number_cpus_cores) 
print(cl)
registerDoSNOW(cl)
clusterEvalQ(cl, library(qtl)) # Load library on all nodes


#(2) Read in the genotype and phenotype data.
cross<-read.cross(format ="csvr", file="Your_Cross_File.csv", estimate.map=F)

#(3) Change population class from bc to RIL if you are working with a RIL population.
cross<-convert2riself(cross)

# --- OPTIMIZATION: Use all cores for Map Estimation ---
newmap=est.map(cross, map.function="kosambi", n.cluster=number_cpus_cores, tol=0.01, maxit=1000)
cross= replace.map(cross, newmap)

#(4) Jittermap
cross<-jittermap(cross)

#(5) Calc Genoprob
cross<-calc.genoprob(cross, map.function="kosambi")

############################
####CHECK PHENOTYPE DATA####
############################

#(6) Generate Histogram plots
pdf(file="Phenotype exploration/1_Phenotype Histograms_transformed.pdf", width=11, height=8.5)
for(i in 1:nphe(cross)) {
	plotPheno(cross, pheno.col=i)
}
dev.off()	

#(7) Shapiro-wilk normality test
norm<-{}
norm2<-{}
for(i in 1:nphe(cross)){
	x<-cross$pheno[i]
	x<-x[,1]
	y<-shapiro.test(x)
	ShapiroWilk.pvalue<-y$p.value
	Phenotype<-colnames(cross$pheno[i])
	norm<-cbind(Phenotype, ShapiroWilk.pvalue)
	norm2<-rbind(norm2, norm)
}
write.csv(file="Phenotype exploration/2_Normality_test_p-values_post_transformation.csv", norm2)

#(8) Batch effects
pdf(file="Phenotype exploration/3_Batch effects.pdf", width=11, height=8.5)
for(i in 1:nphe(cross)){
	par(mfrow=c(1,2), las=1, cex=0.8)
	means<-apply(cross$pheno[i], 1, mean)
	plot(means)
	plot(sample(means), xlab="Random Index", ylab="means", main=colnames(cross$pheno[i]))
}
dev.off()

###########################
####CHECK GENOTYPE DATA####
###########################

#(9) Graphical genotype
pdf(file="Genotype exploration/4_Graphical Genotype.pdf", width=11, height=8.5)
geno.image(cross)
dev.off()

#(10) Genetic map
pdf(file="Genotype exploration/5_Estimated Genetic Map.pdf", width=11, height=8.5)
par(cex=0.6)
plot.map(cross, show.marker.names=FALSE)
dev.off()

#(11.1) Recombination frequencies
png(file="Genotype exploration/6_recombination frequncies.png", width=2400, height=1900, pointsize = 48)
cross <- est.rf(cross)
plotRF(cross, col.scheme=c("redblue"))
dev.off()

#(11.2 - 11.4) RF plots and CSVs
for(i in 1:nchr(cross)){
    file<-paste("Genotype exploration/6_recombination frequncies_ch_", i, ".png", sep="")
    png(filename=file, width=2400, height=1900, pointsize = 48)
    plotRF(cross, chr=c(i), col.scheme=c("redblue"))
    dev.off()
    
    file_rf<-paste("Genotype exploration/6_rf_ch_", i, ".csv", sep="")
    rf <- pull.rf(subset(cross, chr=i))
    write.csv(file=file_rf, rf)
    
    file_lod<-paste("Genotype exploration/6_rf_lod_ch_", i, ".csv", sep="")
    rf_lod <- pull.rf(subset(cross, chr=i), what="lod")
    write.csv(file=file_lod, rf_lod)
}

#(12) Segregation distortion
sd<-geno.table(cross)
sd<-sd[ sd$P.value < 1e-5, ]
write.csv(file="Genotype exploration/7_Chi-square for segregation distortion.csv", sd)

#(13) Genotype comparisons histogram
pdf(file="Genotype exploration/8_Histogram of genotype comparisons.pdf")
genotype.comparisons<-comparegeno(cross)
hist(genotype.comparisons, breaks=200, xlab="Proportion of identical genotypes")
rug(genotype.comparisons)
dev.off()

#(14) Similar genotypes CSV
cg.high<-which(genotype.comparisons>0.95, arr.ind=TRUE)
proportion<-{}
for(i in 1:nrow(cg.high)){
	x<-genotype.comparisons[cg.high[i,1], cg.high[i,2]]
	proportion<-c(proportion, x)
}
Genotype1<-cg.high[,1]
Genotype2<-cg.high[,2]
cg.high<-cbind(Genotype1, Genotype2, proportion)
write.csv(file="Genotype exploration/9_Unusually similiar genotypes.csv", cg.high)

#(15) Poorly typed markers
x<-(nind(cross))
x<-x*0.15
missing.ind<-nmissing(cross, what="mar")
missing.ind<-missing.ind[missing.ind > x]
nmi<-"Number of Individuals without marker data"
missing.ind<-c(nmi, missing.ind)
write.csv(file="Genotype exploration/10_Poorly typed markers.csv", missing.ind)

#(16) Poorly typed individuals
x<-(totmar(cross))
x<-x*0.15
number.missing.mar<-nmissing(cross, what="ind")
Line.ID<-1:nind(cross)
missing.mar<-cbind(Line.ID, number.missing.mar)
missing.mar<-missing.mar[missing.mar[,2]>x,]
write.csv(file="Genotype exploration/11_Poorly typed individuals.csv", missing.mar)

############################################
####INITIAL GENOME SCAN FOR ADDITIVE QTL####
############################################

#(17) Single QTL Scan
cross.sc1<-scanone(cross, pheno.col=1:nphe(cross), method="hk", use="all.obs", model="np")

#(18) Permutations (Heavy computation)
# --- OPTIMIZATION: Use all cores for Permutations ---
cross.sc1.perms<-scanone(cross, pheno.col=1:nphe(cross), method="hk", n.perm=1000, verbose=TRUE, model="np", n.cluster=number_cpus_cores, tol=0.001)

#(19) Summary text files
sum<-summary(cross.sc1, threshold=0, perms=cross.sc1.perms, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=1.5, expandtomarkers=TRUE)
space<-" "
for(i in 1:nphe(cross)){
	x<-capture.output(sum[[i]])
	cat(colnames(cross$pheno[i]), file="QTL files/12_Initial QTL hits by phenotype.txt", sep="\n", append=TRUE)
	cat(x, file="QTL files/12_Summary of top hits by phenotype.txt", sep="\n", append=TRUE)
	cat(space, file="QTL files/12_Summary of top hits by phenotype.txt", sep="\n", append=TRUE)
}

#(20) LOD scores at every marker
data <- data.frame(cross.sc1)
for (i in 1:nphe(cross)){
	lodsi <- data[,c(1,2,i+2)]
	phenotype<-i
	pheno<-colnames(cross$pheno[phenotype])
	file<-paste("13_LOD scores_scan1/13_LOD scores for every marker_", pheno,".txt", sep="")
 	write.table(lodsi, file=file, col.names=NA, sep="\t")
}

#(21) QTL Plots
z<-summary(cross.sc1.perms, alpha=.05)
pdf(file="QTL files/14_QTL Plots.pdf", width=11, height=8.5)
for(i in 1:nphe(cross)){
    plot(cross.sc1, lodcolumn=i, lwd=1.5, gap=0, bandcol="gray70", incl.markers=TRUE, main=colnames(cross$pheno[i]), xlab=c("Threshold for alpha=.05 using 1000 permutations", z[i]))
    add.threshold(cross.sc1, perms=cross.sc1.perms, alpha=0.05, lodcolumn=i, gap=0)
}
dev.off()

#(22) Genetic Map Positions
newmap<-pull.map(cross)
for(i in 1:length(names(newmap))){
	snps<-names(newmap[[i]])
	gm<-c(snps, newmap[[i]])
	gm2<-matrix(gm, ncol=2)
	write.table(file="Genotype exploration/15_Genetic Map Positions.csv", sep=",", append=TRUE, gm2)
}

#(23) sim.geno
cross2<-sim.geno(cross)

######################################################################################
# MAIN PARALLEL LOOP (Optimization for 3000+ phenotypes)
######################################################################################

# Prepare Headers for Master Files in "QTL files" directory
cat("Phenotype,QTL,QTL_Name,QTL_LOD,QTL_Percent_Var", file="QTL files/20_Structured_QTL_Results.csv", sep="\n")
cat("trait,QTL,marker,AA_mean,BB_mean,AA_SE,BB_SE", file="QTL files/22_Structured_Means_SE.csv", sep="\n")
cat("Phenotype,LOD_Threshold_0.05", file="QTL files/23_threshold_for_every_trait.csv", sep="\n")


# Execute Loop in Parallel
clusterExport(cl, c("cross", "cross2", "cross.sc1.perms", "sum"))

results_list <- foreach(k = 1:nphe(cross), .packages="qtl") %dopar% {
    
    phenotype <- k
    pheno <- colnames(cross$pheno[phenotype])
    space <- " "
    
    # Strings to hold output for shared files
    out_16 <- c()
    out_17 <- c()
    out_19 <- c()
    out_20_txt <- c()
    out_20_csv <- c()
    out_22_txt <- c()
    out_22_csv <- c()
    
    # --- FILE 23: CALCULATE THRESHOLD FOR EVERY TRAIT ---
    sub.perms <- subset(cross.sc1.perms, lodcolumn=phenotype)
    perm_thresh <- summary(sub.perms, alpha=0.05)[1]
    out_23_csv <- paste(pheno, perm_thresh, sep=",")
    
    # (25) Make QTL object
    chromo<-sum[[pheno]]
    chr<-{}
    pos<-{}
    if(nrow(chromo) > 0){
        for(i in 1:nrow(chromo)){
            chr1<-chromo[i,1]
            chr2<-as.numeric(as.character(chr1))
            chr<-c(chr, chr2)
            pos1<-chromo[i,2]
            pos<-c(pos, pos1)
        }
    }
    
    if(!is.null(chr) && length(chr) > 0 && !is.na(chr[1])) {
        qtl<-makeqtl(cross, chr=chr, pos=pos, what="prob")
        createqtl<- paste("Q", 1:qtl$n.qtl, sep="")
        formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "+")))
        
        # (26) Scan for additional QTL
        cross.aq<-addqtl(cross, pheno.col=phenotype, qtl=qtl, formula=formula, method="hk", model="normal")
        
        xx<-capture.output(summary(cross.aq, perms=sub.perms, alpha=.05, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=1.5, expandtomarkers=TRUE))
        
        out_16 <- c(pheno, xx, space)
        
        # (27) Add additional QTL
        sum.aq<-summary(cross.aq, perms=sub.perms, alpha=.05, pvalues=TRUE, format="tabByCol", ci.function="lodint", drop=0.95, expandtomarkers=TRUE)
        
        chr.aq <- NULL
        pos.aq <- NULL
        
        if(length(sum.aq) > 0 && "lod" %in% names(sum.aq)) {
             new_chr <- sum.aq$lod[, 1]
             new_pos <- sum.aq$lod[, 2]
             chr.aq <- c(chr, new_chr[!is.na(new_chr)])
             pos.aq <- c(pos, new_pos[!is.na(new_pos)])
        } else {
             chr.aq <- chr
             pos.aq <- pos
        }

        qtl<-makeqtl(cross, chr=chr.aq, pos=pos.aq, what="prob")
        createqtl<- paste("Q", 1:qtl$n.qtl, sep="")
        formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "+")))
        
        # (28) Stepwise Selection
        pen<-summary(sub.perms)
        cross.sw<-stepwiseqtl(cross, pheno.col=phenotype, qtl=qtl, formula=formula, method="hk", penalties=pen, model="normal", additive.only=TRUE)
        swQTL<-capture.output(print(cross.sw))
        
        out_17 <- c(pheno, swQTL, space)
        
        # (29) Refine QTL
        sum.sw<-summary(cross.sw)
        if (length(sum.sw)==0){
             # No QTLs found after stepwise
             null_msg <- "	There were no LOD peaks above the threshold"
             out_16 <- c(out_16, null_msg, space)
             out_17 <- c(out_17, null_msg, space)
             # out_19 removed (empty for nulls)
             out_20_txt <- c(pheno, null_msg, space)
             out_22_txt <- c(pheno, null_msg, space)
        } else {
            chr.sw<-{}
            pos.sw<-{}
            for(i in 1:nrow(sum.sw)){
                chr.sw1<-sum.sw[i,2]
                chr.sw2<-as.numeric(as.character(chr.sw1))
                chr.sw<-c(chr.sw, chr.sw2)
                pos1.sw<-sum.sw[i,3]
                pos.sw<-c(pos.sw, pos1.sw)
            }

            qtl2<-makeqtl(cross, chr=chr.sw, pos=pos.sw, what="prob")
            createqtl<- paste("Q", 1:qtl2$n.qtl, sep="")
            formula<-as.formula(paste("y ~ ", paste(createqtl, collapse= "+")))
            rqtl<-refineqtl(cross, pheno.col=phenotype, qtl=qtl2, method="hk", model="normal")
            
            # (30) Write file 18 (Unique per pheno -> Direct Write to Folder)
            model_lods <- list()
            lodprof <- attr(rqtl,"lodprofile")
            if(!is.null(lodprof)){
                for (i in 1:length(lodprof)){
                    model_lods[[i]] <- lodprof[[i]]
                }
                model_lods.df <- as.data.frame(do.call(rbind, model_lods))
                file18 <- paste("18_LOD scores_model part/18_LOD scores at every marker of model part ", pheno, ".txt", sep="")
                write.table(model_lods.df, file=file18, col.names=NA, sep="\t")
            }

            # (31) File 19 Summary
            for (i in 1:rqtl$n.qtl){
                interval<-capture.output(lodint(rqtl, qtl.index=i, drop=1.5, expandtomarkers=TRUE))
                q_name<-paste("Q", i, sep="")
                out_19 <- c(out_19, pheno, q_name, interval, space)
            }
            
            # (32) FitQTL and File 20 (ANOVA)
            cross.ests<-fitqtl(cross, pheno.col=phenotype, qtl=rqtl, formula=formula, method="hk", dropone=TRUE, get.ests=TRUE, model="normal")
            ests<-capture.output(summary(cross.ests))
            out_20_txt <- c(pheno, ests, space)
            
            # --- STRUCTURED CSV 20 (Long Format + Smart Model Row) ---
            sum.fit <- summary(cross.ests)
            full.stats <- sum.fit$result.full
            model.lod <- full.stats["Model", "LOD"]
            model.var <- full.stats["Model", "%var"]

            if("result.drop" %in% names(sum.fit)){
                drop.table <- sum.fit$result.drop
                for(r in 1:nrow(drop.table)){
                    qtl.label <- paste("Q", r, sep="")
                    qtl.name <- rownames(drop.table)[r]
                    qtl.lod  <- drop.table[r, "LOD"]
                    qtl.var  <- drop.table[r, "%var"]
                    out_20_csv <- c(out_20_csv, paste(pheno, qtl.label, qtl.name, qtl.lod, qtl.var, sep=","))
                }
                out_20_csv <- c(out_20_csv, paste(pheno, "model", "model", model.lod, model.var, sep=","))
            } else {
                qtl.label <- "Q1"
                qtl.name <- rqtl$name
                out_20_csv <- c(out_20_csv, paste(pheno, qtl.label, qtl.name, model.lod, model.var, sep=","))
            }
            
            # (33) Effects Plots and File 22
            file21 <- paste("21_Marker effect plots/21_Marker effect plots ", pheno, ".pdf", sep="")
            pdf(file=file21, width=11, height=8.5)
            
            for(i in 1:length(chr.sw)) {
                b<-paste("Q",i, sep="")
                mar<-find.marker(cross2, chr=chr.sw[i], pos=pos.sw[i])
                plotPXG(cross, marker=mar, pheno.col=phenotype)
                
                phenoqtl<-paste(pheno, b)
                means<-effectplot(cross2, pheno.col=phenotype, mname1=mar, draw=FALSE)
                
                out_22_txt <- c(out_22_txt, phenoqtl, space)
                out_22_txt <- c(out_22_txt, capture.output(print(means$Means)), space, capture.output(print(means$SEs)), space)

                # --- STRUCTURED CSV 22 (Wide Format) ---
                mean_aa <- means$Means[1]
                mean_bb <- means$Means[2]
                se_aa   <- means$SEs[1]
                se_bb   <- means$SEs[2]
                out_22_csv <- c(out_22_csv, paste(pheno, b, mar, mean_aa, mean_bb, se_aa, se_bb, sep=","))
            }
            dev.off() 
        }
    } else {
        null_msg <- "	There were no LOD peaks above the threshold"
        out_16 <- c(pheno, null_msg, space)
        out_17 <- c(pheno, null_msg, space)
        out_20_txt <- c(pheno, null_msg, space)
        out_22_txt <- c(pheno, null_msg, space)
    }
    
    list(out16=out_16, out17=out_17, out19=out_19, 
         out20txt=out_20_txt, out20csv=out_20_csv, 
         out22txt=out_22_txt, out22csv=out_22_csv,
         out23csv=out_23_csv)
}

# Shut down cluster
stopCluster(cl)

# Post-Loop: Write Shared Files sequentially to "QTL files" folder
print("Writing shared files...")

for(res in results_list) {
    if(!is.null(res$out16)) cat(res$out16, file="QTL files/16_Additional QTL hits by phenotype.txt", sep="\n", append=TRUE)
    if(!is.null(res$out17)) cat(res$out17, file="QTL files/17_Additional QTL hits from stepwise analysis.txt", sep="\n", append=TRUE)
    if(!is.null(res$out19)) cat(res$out19, file="QTL files/19_Summary of Final QTL Intervals.txt", sep="\n", append=TRUE)
    if(!is.null(res$out20txt)) cat(res$out20txt, file="QTL files/20_ANOVA results and QTL effect estimates.txt", sep="\n", append=TRUE)
    if(!is.null(res$out20csv)) cat(res$out20csv, file="QTL files/20_Structured_QTL_Results.csv", sep="\n", append=TRUE)
    if(!is.null(res$out22txt)) cat(res$out22txt, file="QTL files/22_means and SE.txt", sep="\n", append=TRUE)
    if(!is.null(res$out22csv)) cat(res$out22csv, file="QTL files/22_Structured_Means_SE.csv", sep="\n", append=TRUE)
    if(!is.null(res$out23csv)) cat(res$out23csv, file="QTL files/23_threshold_for_every_trait.csv", sep="\n", append=TRUE)
}

print("Analysis Complete.")
##########################################################################################