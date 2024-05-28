
###--------------------------------------------------------------------------###
### environment settings
###--------------------------------------------------------------------------###

rm(list=ls())
options(stringsAsFactors = FALSE)

library(GWASTools)
# library(gdata)

###--------------------------------------------------------------------------###
### This script is loading previouly gds files and starting to do QC including
###
###   (1) Batch Quality Check (Chapter 3)
###   (2) Sample Quality Check (Chapter 4)
###   (3) Sample Identity Check (Chapter 5)
###
###   The first chunck is the way to load the gds files into R 
###       we will use this chunck of the code for all the QC in the future
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
### (1) Loading gds files
###--------------------------------------------------------------------------###

setwd("/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS5")
# Load the annotation files
load("GWAS_5_Annotation.RData")

# assign file name and the path of the gds files to an object 
geno.file <- "./GWAS5_geno.gds"
qxy.file <- "./GWAS5_qxy.gds"
bl.file <- "./GWAS5_bl.gds"

# get the genotype from the gds file created previously
gds1 <- GdsGenotypeReader(geno.file)
# annotated with scan and SNP annotation
genoData <- GenotypeData(gds1, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
genoData
# close(genoData)

# get the quality score, X and X intensity from the gds file created previously
gds2 <- GdsIntensityReader(qxy.file)
# annotated with scan and SNP annotation
qxyData <- IntensityData(gds2, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
qxyData
#close(qxyData)

# get the BAF and LRR from the gds file created previously
gds3 <- GdsIntensityReader(bl.file)
# annotated with scan and SNP annotation
blData <- IntensityData(gds3, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
blData
# close(blData)

###--------------------------------------------------------------------------###
### (2) Batch QC -- Chapter 3
###--------------------------------------------------------------------------###

##########################################################
### Chapter 3.1 Missing call rate for Samples and SNP  ###
##########################################################

##########################################################
### Calculate the number of missing calls for each snp ###  
###   over all samples for each sex separately         ###
##########################################################
miss <- missingGenotypeBySnpSex(genoData)

### Examine the results ###
# names(miss)
# head(miss$missing.counts)
# miss$scans.per.sex
# head(miss$missing.fraction)

### Add missing call rates into snp annotation file ###
# Make sure ordering matches snp annotation
allequal(snpAnnot$snpID, as.numeric(names(miss$missing.fraction)))
# TRUE

# Add a column for missing fraction  
snpAnnot$missing.n1 <- miss$missing.fraction
# Add metadata to describe the column
varMetadata(snpAnnot)["missing.n1", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all samples",
  "except that females are excluded for Y chr SNPs")

# Calculate summary of the missing call rate
# summary(snpAnnot$missing.n1)

# Plot the Missing call rate for all probess
pdf("./QC_check/Ch3.1-1 SNP missing call rate (n1).pdf")
hist(snpAnnot$missing.n1, probability = TRUE,
     xlab="SNP missing call rate", breaks = seq(0,1,0.01),
     main="Missing Call Rate for All Probes", las=1)
dev.off()

###########################
### Examine the results ###
###########################
# Find the number of SNPs with every call missing
length(snpAnnot$missing.n1[snpAnnot$missing.n1 == 1])
# 9935 SNPs
# ind <- which(snpAnnot$missing.n1 == 1)
# snpAnnot@data[ind,]
# write.table(snpAnnot@data[ind,], "./QC_check/number_SNPs_missing.n1_1.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Fraction of autosomal SNPs with missing call rate < 0.05
x <- snpAnnot$missing.n1[snpAnnot$chromosome < 23]
length(x[x < 0.05]) / length(x)
# 0.9920817

# Fraction of X chromosome SNPs with missing call rate < 0.05
x <- snpAnnot$missing.n1[snpAnnot$chromosome == 23]
length(x[x < 0.05]) / length(x)
# 0.9985893

# Fraction of Y chromosome SNPs with missing call rate < 0.05
x <- snpAnnot$missing.n1[snpAnnot$chromosome == 25]
length(x[x < 0.05]) / length(x)
# 0.9458958

##################################################################
### Calculate missing call rates for each sample by chromosome ###  
##################################################################

# Want to exclude all SNP probes with 100% missing call rate
# Check on how many SNPs to exclude
sum(snpAnnot$missing.n1 == 1)
# 9935

# Create a variable that contains the IDs of these SNPs to exclude
snpexcl <- snpAnnot$snpID[snpAnnot$missing.n1 == 1]
length(snpexcl)
# 9935

# Calculate the missing call rate per sample by chromosme 
miss <- missingGenotypeByScanChrom(genoData, snp.exclude=snpexcl)

### Examine the results ###
# names(miss)
# head(miss$missing.counts)
# head(miss$snps.per.chr)
# head(miss$missing.fraction)

# Check to make sure that the correct number of SNPs were excluded
# sum(miss$snps.per.chr)
# nrow(snpAnnot) - sum(miss$snps.per.chr)

### Add missing call rates into scan annotation file ###
# Check the ordering matches the sample annotation file
allequal(names(miss$missing.fraction), scanAnnot$scanID)

# Add the missing call rates vector to the sample annotation file
scanAnnot$missing.e1 <- miss$missing.fraction
varMetadata(scanAnnot)["missing.e1", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all snps with missing.n1<1",
  "except that Y chr SNPs are excluded for females")

# Calculate summary of the missing call rate
# summary(scanAnnot$missing.e1)

# Plot the Missing call rate for all samples
pdf("./QC_check/Ch3.1-2 Sample missing call rate (e1).pdf")
hist(scanAnnot$missing.e1,
     xlab="Fraction of missing calls over all probes",
     main="Histogram of Sample Missing Call Rate for all Samples",
     las=1,  probability = TRUE,  breaks = seq(0,1,0.01))
dev.off()

###########################
### Examine the results ###
###########################

# Look at missing.e1 for males
summary(scanAnnot$missing.e1[scanAnnot$sex == "M"])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002449 0.0004782 0.0007141 0.0010552 0.0009817 0.0058427 
# Look at missing.e1 for females
summary(scanAnnot$missing.e1[scanAnnot$sex == "F"])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002383 0.0005307 0.0008267 0.0013743 0.0011959 0.0404767

# Number of samples with missing call rate > 5%
sum(scanAnnot$missing.e1 > 0.05)
# 0

#######################################################################
### Calcualate the missing call rate for autosomes and              ###
###   the X chromosome separately                                   ###
### These values were calculated and added into the scan annotation ###
###   also, a logical "duplicated" variable was created             ###
#######################################################################

# Autosome 
auto <- colnames(miss$missing.counts) %in% 1:22
missa <- rowSums(miss$missing.counts[,auto]) / sum(miss$snps.per.chr[auto])
summary(missa)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002337 0.0005180 0.0007807 0.0012463 0.0011030 0.0407766 

# X chromosome
missx <- miss$missing.counts[,"X"] / miss$snps.per.chr["X"]
summary(missx)
#       Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#  0.0001603 0.0004969 0.0007133 0.0011898 0.0010500 0.0308102

# check they match sample annotation file
allequal(names(missa), scanAnnot$scanID)
# TRUE
allequal(names(missx), scanAnnot$scanID)
# TRUE

# Add these separate sample missing call rates to the scan annotation
scanAnnot$miss.e1.auto <- missa
scanAnnot$miss.e1.xchr <- missx

# Order scanAnnot by missing.e1 so duplicate subjectIDs
# with a higher missing rate are marked as duplicates

scanAnnot <- scanAnnot[order(scanAnnot$scanName, scanAnnot$missing.e1),]
scanAnnot$duplicated <- duplicated(scanAnnot$scanName)

table(scanAnnot$duplicated, useNA="ifany")
# FALSE  
#   248   
# no duplicates found

# Put scanAnnot back in scanID order; this is very important!!
scanAnnot <- scanAnnot[order(scanAnnot$scanID),]
allequal(scanAnnot$scanID, sort(scanAnnot$scanID))
varMetadata(scanAnnot)["duplicated", "labelDescription"] <- "TRUE for duplicate scan with higher missing.e1"

###########################################################
### Calculate missing call rate per SNP with missing.e1 ###
###   less than 0.05 over all samples.                  ###
###   (Smaples with missing call rate less than 0.05)   ###
###########################################################
# Find the samples with missing.e1 > .05 and make a vector of
# scanID to exclude from the calculation
scan.exclude <- scanAnnot$scanID[scanAnnot$missing.e1 > 0.05]

length(scan.exclude)
# 0

#######################################################################
### Extract the high missing call rate samples
#######################################################################

# no samples with high missing call rate
# d <- pData(scanAnnot)

# himis <- d[scan.exclude,c("corrected_scanName","sex","race","age","edu","center",
#                           "plate","well_position","missing.e1")]

# write.table(himis, "High_Missing_Call_Rate_Samples.txt", quote=FALSE, row.names = FALSE, sep="\t")

# Call missingGenotypeBySnpSex and save the output
miss <- missingGenotypeBySnpSex(genoData, scan.exclude=scan.exclude)

### Add missing call rates into snp annotation file ###
# Add the missing call rates vector to the sample annotation file
snpAnnot$missing.n2 <- miss$missing.fraction
varMetadata(snpAnnot)["missing.n2", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all samples with missing.e1<0.05",
  "except that females are excluded for Y chr SNPs")

# Calculate summary of the missing call rate
summary(snpAnnot$missing.n2)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#  0.000000 0.000000 0.000000 0.006564 0.000000 1.000000  

##############################################################
### Calculate missing call rate per sample with missing.n2 ###
###   less than 0.05 over all snps.                        ###
###   (SNP with missing call rate less than 0.05)          ###
##############################################################
# Create a vector of the SNPs to exclude.
snpexcl <- snpAnnot$snpID[snpAnnot$missing.n2 >= 0.05]
# length(snpexcl)
# 14962
miss <- missingGenotypeByScanChrom(genoData, snp.exclude=snpexcl)
### Add missing call rates into snp annotation file ###
# Add the missing call rates vector to the sample annotation file
scanAnnot$missing.e2 <- miss$missing.fraction
varMetadata(scanAnnot)["missing.e2", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all snps with missing.n2<0.05",
  "except that Y chr SNPs are excluded for females")

# Calculate summary of the missing call rate
summary(scanAnnot$missing.e2)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002108 0.0004058 0.0005909 0.0010693 0.0009294 0.0396922 

# Plot the Missing call rate for all samples (with SNP missing call rate < 0.05)

pdf("./QC_check/Ch3.1-3 Sample missing call rates with SNP missing rate < 0.05 (e2).pdf")
hist(scanAnnot$missing.e2, xlab="Fraction of missing calls over all probes
     with missing call rate < 0.05",
     main="Histogram of Sample Missing Call Rate for all Samples",
     las=1,  probability = TRUE,  breaks = seq(0,1,0.01))
dev.off()

##################################################
### Chapter 3.2 Missing call rate for by batch ###
##################################################

# do not incorporate batch information, skip this part for now

# varLabels(scanAnnot)
# Check how many batches exist and how many samples are in each batch
# length(unique(scanAnnot$plate))
table(scanAnnot$plate, useNA="ifany")
#  1  2  3 
# 96 96 56 

### Plot the distribution of the number of samples per batch.
pdf("./QC_check/Ch3.2-1 Sample per plate.pdf")
barplot(table(as.integer(scanAnnot$plate)),
        ylab="Number of Samples", xlab="Plate",
        main="Distribution of Samples per Batch")
dev.off()

### Examine the mean missing call rate per batch for all SNPs
batches <- unique(scanAnnot$plate)
bmiss <- rep(NA,length(batches))
names(bmiss) <- batches
bn <- rep(NA,length(batches)) 
names(bn) <- batches

for(i in 1:length(batches)) {
  x <- scanAnnot$missing.e1[is.element(scanAnnot$plate, batches[i])]
  bmiss[i] <- mean(x)
  bn[i] <- length(x)
}

# check relationship between sample size and missign call rate by batches
y <- lm(bmiss ~ bn)
anova(y)

# Analysis of Variance Table
# Response: bmiss
# Df     Sum Sq    Mean Sq F value Pr(>F)
# bn         1 1.7500e-09 1.7500e-09  0.0045 0.9572
# Residuals  1 3.8666e-07 3.8666e-07   

# Plot missing call rate ~ number of samples by plate
pdf("./QC_check/Ch3.2-2 Batch Effect.pdf")
plot(bn, bmiss,
     xlab="Number of samples per batch", ylab="Mean missing call rate",
     main="Mean Missing Call Rate vs\nSamples per Batch")
abline(y$coefficients)
dev.off()

####################################################
### Check Chi-square statistics for batch effect ###
####################################################
res <- batchChisqTest(genoData, batchVar="plate", return.by.snp=TRUE)

# chi-square values for each SNP
dim(res$chisq)
# 1828791       3
# genomic inflation factor
res$lambda
#         1         2         3 
# 0.7497860 0.8385354 0.3160044
# average chi-square test statistic for each of the batches
res$mean.chisq
#         1         2         3 
# 1.1627117 1.3347863 0.6068725

# Find out which plate having highest chi-square test statistics 
which.max(res$mean.chisq) 
# 2
# 2
max(res$mean.chisq)
# 1.334786

# Find out which plate having highest lambda 
which.max(res$lambda)
# 2
# 2
max(res$lambda)
# 0.8385354

#######################################################################
### Chapter 3.3 Chisquare test of Allele Freq Difference in Batches ### 
#######################################################################

# GWAS5 cohort is all White

# check races for each sample
x <- table(scanAnnot$race, useNA="ifany")

# Calculate the proportion of each race
x / sum(x)

# Count the races for each plate/bathes
x <- table(scanAnnot$race, scanAnnot$plate)
x

# Run an approximate chi-square test to see if there are ethnic effects between batches
chisq <- chisq.test(x)
chisq$p.value

# Calculate the fraction of samples in each batch that are White
batches <- unique(scanAnnot$plate)
eth <- rep(NA,length(batches))
names(eth) <- sort(batches)

for(i in 1:length(batches)){
  x <- scanAnnot$race[is.element(scanAnnot$plate, batches[i])]
  xl <- length(x[x == "White"])
  eth[i] <- xl / length(x)
}

# allequal(names(eth), names(res$mean.chisq))
# TRUE

# Plot the average Chi-Square test statistic against the
#     fraction of samples that are CEU
pdf("./QC_check/Ch3.3-1 Races vs Batches.pdf")
plot(eth, res$mean.chisq, xlab="Fraction of White Samples per Batch",
     ylab="Average Chi-square Test Statistic",
     main="Fraction of White Samples per Batch
     vs Average Chi-square Test Statistic")
abline(v=mean(eth), lty=2, col="red")
dev.off()

###--------------------------------------------------------------------------###
### (3) Sample QC -- Chapter 4
###--------------------------------------------------------------------------###

# genoData (gds1) :  genotype data
# qxyData (gds2) : quality data
# blData (gds3) : BL data

###################################################
### Chapter 4.1 Median quality score            ###
###################################################

# Calculate the median quality score for each sample
qual.results <- qualityScoreByScan(qxyData, genoData, snp.exclude = snpexcl)

# no high missing call rate samples, no need to exclude
# exc.qual.results <- qual.results[!row.names(qual.results)%in%scan.exclude,]

# Plot the histogram of the median quality score for each sample
# pdf("./QC_check/Ch4.1-1 Median Genotype Quality Scores of Samples.pdf")
# hist(exc.qual.results[,"median.quality"], main="Median Genotype Quality Scores
#      of Samples", xlab="Median Quality", las=1, ylim=c(0,2500), xlim=c(0.82, 0.84),
#      breaks = seq(0.82,0.84,0.001), axes=FALSE)
# axis(1, at=seq(0.82,0.84,0.001), cex.axis=0.7, las=2)
# axis(2, at=seq(0,2500,500), cex.axis=0.7, las=1)
# text(0.823, 1500, paste0("Min = ", round(min(exc.qual.results[,"median.quality"]),4))) 
# dev.off()

# Plot the histogram of the median quality score for all samples
pdf("./QC_check/Ch4.1-1 Median Genotype Quality Scores of samples.pdf")
hist(qual.results[,"median.quality"], 
     main="Median Genotype Quality Scores + of Samples", 
     xlab="Median Quality")
dev.off()

########################################################
### Chapter 4.2 B Allele Frequency variance analysis ###
########################################################

# nbins: a vector with integers corresponding to the number of bins for each chromosome. The values all must be even integers
nbins <- rep(12,23) 
# The length of nbins must be equal to one more than the number of autosomes in the genotype netcdf file

# This process identifes chromosome-sample pairs that have windows with very high BAF standard
# deviation, with "very high" defned as more than 4 standard deviations from the window's mean
# BAF standard deviation over all samples.
slidingBAF12 <- sdByScanChromWindow(blData, genoData, nbins=nbins, snp.exclude = snpexcl)

# names(slidingBAF12)
# dim(slidingBAF12[["2"]])

# calculates the mean and standard deviation of the BAF standard deviations in each window in each chromosome 
# over all samples. For the X chromosome, males and females are calculated separately, and we save the results split by sex.
sds.chr <- meanSdByChromWindow(slidingBAF12, scanAnnot$sex)

# sds.chr[["21"]]
# sds.chr[["X"]]

# identify windows within sample-chromosome pairs that have very high BAF standard
# deviations compared to the same window in other samples
res12bin4sd <- findBAFvariance(sds.chr, slidingBAF12, scanAnnot$sex, sd.threshold=4)

# head(res12bin4sd)
# table(res12bin4sd[, "chromosome"])

# Call chromIntensityPlot to plot the BAF of all SNPs on the indicated chromosome-sample pairs against position.
scanID <- res12bin4sd[, "scanID"]
chrom <- res12bin4sd[, "chromosome"]
chrom[res12bin4sd[, "chromosome"] == "X"] <- 23
bincode <- paste("Bin", res12bin4sd[, "bin"], sep = " ")

# Plot the results 
# Don't RUN !!!

# pdf("./QC_check/Ch4.2-1 High BAF sd regions.pdf")
# chromIntensityPlot(blData, scanID, chrom, info=bincode, ideogram=TRUE)
# dev.off()
# close(blData)

# Find out numbers of bad bins per sample

x <- sort(table(scanID), decreasing=TRUE)

y <- x[x>1]
col <- rep("black",length(y))
col[names(y)%in%scan.exclude] <- "red"

pdf("./QC_check/Ch4.2-2 Number of bad bins per sample.pdf")
barplot(y, axisnames = FALSE, axes = FALSE)
axis(2, at=seq(0,21,3), las=1, cex.axis=0.8)
Map(axis, side=1, at=(1:length(y))*1.2, col.axis=col, labels=names(y), lwd=0, las=2, cex.axis=0.7)
axis(1,at=(1:length(y))*1.2,labels=FALSE)
abline(h=seq(3,21,3), col="gray", lty="3434")
dev.off()

# Find the scanID and scanNames for those having lots of bad bins (>=15)

exclude.badbin <- as.integer(names(x)[x>=15])
length(exclude.badbin)
# 2
scan.exclude.2 <- sort(union(scan.exclude, exclude.badbin))
length(scan.exclude.2)
# 2

#########################################################
### 4.3 Missingness and heterozygosity within samples ###
#########################################################

# This step calculates the percent of missing and heterozygous genotypes in each chromosome of each sample.
miss <- missingGenotypeByScanChrom(genoData, snp.exclude = snpexcl)
# Calcualte the missing call rate for all samples
miss.rate <- t(apply(miss$missing.counts, 1, function(x) {
  x / miss$snps.per.chr}))

# no high missing rate samples
# ind <- which(row.names(miss$missing.counts)%in%scan.exclude)
# miss.rate <- as.data.frame(miss.rate)[-ind,]
miss.rate <- as.data.frame(miss.rate)

# Plot missingness in the autosomes
# Select autosomes 
cols <- names(miss.rate) %in% c(1:22,"X","XY")
pdf("./QC_check/Ch4.3-1 Missingness by Chromosome.pdf")
boxplot(miss.rate[,cols], main="Missingness by Chromosome",
        ylab="Proportion Missing", xlab="Chromosome", las=1,
        pch=16, cex.axis=0.7, cex=0.7)
dev.off()

# Plot X chromosome missingness for each sex
pdf("./QC_check/Ch4.3-2 X Chromosome Missingness.pdf")
boxplot(miss.rate$X ~ scanAnnot$sex,
        main="X Chromosome Missingness by Sex",
        ylab="Proportion Missing", las=1, pch=16, cex=0.7)
dev.off()

# Heterozygosity !!!
# Calculate heterozygosity by scan by chromosome
het.results <- hetByScanChrom(genoData, snp.exclude = snpexcl)

# Ensure heterozygosity results are ordered correctly
het.results2 <- het.results[order(as.integer(rownames(het.results))),]
allequal(scanAnnot$scanID, rownames(het.results2))

# Write autosomal and X chr heterozygosity to sample annot
scanAnnot$het.A <- het.results2[,"A"]
scanAnnot$het.X <- het.results2[,"X"]
varMetadata(scanAnnot)["het.A", "labelDescription"] <-
  "fraction of heterozygotes for autosomal SNPs"
varMetadata(scanAnnot)["het.X", "labelDescription"] <-
  "fraction of heterozygotes for X chromosome SNPs"

"%nin%" <- function(x,y){
  !x%in%y
}

# Plot heterozygosity over the autosomes
pdf("./QC_check/Ch4.3-3 Autosomal Heterozygosity.pdf")
boxplot(scanAnnot$het.A ~ scanAnnot$race,
        main="Autosomal Heterozygosity", las=1, ylim=c(0,0.151),
        ylab="% of heterozygosity")
dev.off()

# Plot female heterozygosity on the X chromosome
pdf("./QC_check/Ch4.3-4 X Chromosome Heterozygosity in Females.pdf")
ind <- which(scanAnnot$sex == "F")
boxplot(scanAnnot$het.X[ind] ~ scanAnnot$race[ind],
        main="X Chromosome Heterozygosity in Females",  ylim=c(0,0.151),
        las=1, ylab="% of heterozygosity")
dev.off()

### Examine results ###
# range(scanAnnot$het.A[ind])
# range(scanAnnot$het.X[ind])

###--------------------------------------------------------------------------###
### (4) Identity QC -- Chapter 5
###--------------------------------------------------------------------------###

inten.by.chrom <- meanIntensityByScanChrom(qxyData, snp.exclude = snpexcl)
# names(inten.by.chrom)

mninten <- inten.by.chrom[[1]]  # mean intensities
# dim(mninten)

### Check to be sure sample ordering is consistent
mninten2 <- as.data.frame(mninten[order(as.integer(rownames(mninten))),])
allequal(scanAnnot$scanID, rownames(mninten2))
# TRUE
mninten2$scanID <- row.names(mninten2)

# Merge the annotation file and the mean intensity file
mean_inten <- merge(pData(scanAnnot), mninten2, by="scanID", all=TRUE)

# if the scan.exclude is not NULL, remove those scans
# Use this as a primary data set
# mean_inten_2 <- mean_inten[-scan.exclude,]

#####################################################################
### example mis-labeled sex plate by plate, race by race
#####################################################################

nx <- sum(snpAnnot$chromosome == 23)
ny <- sum(snpAnnot$chromosome == 25)

pdf("./QC_check/Ch5.1-2 Mislabeled-Sex check by plate.pdf")
par(mfrow=c(2,2))

for(i in 1:max(as.integer(mean_inten$plate))){
  
  # Assign each sex a color for each plate
  temp <- subset(mean_inten, plate==i)
  xcol <- rep(NA, nrow(temp))
  xcol[temp$sex == "M"] <- "blue"
  xcol[temp$sex == "F"] <- "red"
  
  #All intensities
  x1 <- temp$X
  y1 <- temp$Y
  main1 <- paste0("Mean X vs \nMean Y Chromosome Intensity \n on plate ", i)
  #Het on X vs X intensity
  x2 <- temp$X
  y2 <- temp$het.X
  main2 <- paste0("Mean X Chromosome Intensity vs \nMean X Chromosome Heterozygosity \n on plate ", i)
  # Het on X vs Y intensity
  y3 <- temp$Y
  x3 <- temp$het.X
  main3 <- paste0("Mean X Chromosome Heterozygosity vs \nMean Y Chromosome Intensity \n on plate ", i)
  # X vs A het
  x4 <- temp$het.A[temp$sex == "F"]
  y4 <- temp$het.X[temp$sex == "F"]
  main4 <- paste0("Mean Autosomal Heterozygosity vs \nMean X Chromosome Heterozygosity \n on plate ", i)
  cols <- c("blue","red")
  mf <- c("male", "female")
  xintenlab <- paste("X intensity (n=", nx, ")", sep="")
  yintenlab <- paste("Y intensity (n=", ny, ")", sep="")
  
  plot(x1, y1, xlab=xintenlab, ylab=yintenlab,
       main=main1, col=xcol, cex.main=0.8)
  legend("topright",mf,col=cols,pch=c(1,1))
  plot(x2, y2, col=xcol, xlab=xintenlab,
       ylab="X heterozygosity", main=main2, cex.main=0.8)
  plot(x3, y3, col=xcol, ylab=yintenlab,
       xlab="X heterozygosity", main=main3, cex.main=0.8)
  plot(x4,y4, col="red", xlab="Autosomal heterozygosity",
       ylab="X heterozygosity", main=main4, cex.main=0.8)
  
}

dev.off()

#####################################################################
### exam mis-labeled sex in total
#####################################################################

# annotation file without excluded scans
# scanAnnot2 <- pData(scanAnnot)[-scan.exclude,]
scanAnnot2 <- pData(scanAnnot)

# Assign each sex a color
xcol <- rep(NA, nrow(scanAnnot2))
xcol[scanAnnot2$sex == "M"] <- "blue"
xcol[scanAnnot2$sex == "F"] <- "red"

nx <- sum(snpAnnot$chromosome == 23)
ny <- sum(snpAnnot$chromosome == 25)

#All intensities
x1 <-mean_inten[,"X"] 
y1 <- mean_inten[,"Y"]
main1 <- "Mean X vs \nMean Y Chromosome Intensity"
#Het on X vs X intensity
x2 <- mean_inten[,"X"]
y2 <- scanAnnot2$het.X
main2 <- "Mean X Chromosome Intensity vs
Mean X Chromosome Heterozygosity"
# Het on X vs Y intensity
y3 <- mean_inten[,"Y"]
x3 <- scanAnnot2$het.X
main3 <- "Mean X Chromosome Heterozygosity vs
Mean Y Chromosome Intensity"
# X vs A het
x4 <- scanAnnot2$het.A[scanAnnot2$sex == "F"]
y4 <- scanAnnot2$het.X[scanAnnot2$sex == "F"]
main4 <- "Mean Autosomal Heterozygosity vs
Mean X Chromosome Heterozygosity"
cols <- c("blue","red")
mf <- c("male", "female")
xintenlab <- paste("X intensity (n=", nx, ")", sep="")
yintenlab <- paste("Y intensity (n=", ny, ")", sep="")
pdf("./QC_check/Ch5.1-1 DataCleaning-sex.pdf")
par(mfrow=c(2,2))
plot(x1, y1, xlab=xintenlab, ylab=yintenlab,
     main=main1, col=xcol, cex.main=0.8)
legend("topright",mf,col=cols,pch=c(1,1))
plot(x2, y2, col=xcol, xlab=xintenlab,
     ylab="X heterozygosity", main=main2, cex.main=0.8)
plot(x3, y3, col=xcol, ylab=yintenlab,
     xlab="X heterozygosity", main=main3, cex.main=0.8)
plot(x4,y4, col="red", xlab="Autosomal heterozygosity",
     ylab="X heterozygosity", main=main4, cex.main=0.8)
dev.off()

#####################################################################
### Fine the IDs with mislabeled sex
#####################################################################

ind1 <- which(mean_inten$X>1 & mean_inten$sex=="M")
ind2 <- which(mean_inten$Y>0.5 & mean_inten$sex=="F")

ind <- sort(c(ind1, ind2))

mis_sex <- subset(mean_inten[ind,], select=c("scanName","sex","race","AlzID","STudyID","LabID",
                                               "plate","well", "Case_Control", "Disease_Status",
                                             "APOE","file"))

write.table(mis_sex, "Mislabeled Sex IDs.txt", quote=FALSE, row.names = FALSE, sep="\t")

###

ind3 <- which(mean_inten$het.X<0.05 & mean_inten$sex=="F")
ind4 <- ind3[ind3%nin%ind]

mis_1x <- subset(mean_inten[ind4,], select=c("scanName","sex","race","AlzID","STudyID","LabID",
                                             "plate","well", "Case_Control", "Disease_Status",
                                             "APOE","file"))


# no output here
# write.table(mis_1x, "X singelton female.txt", quote=FALSE, row.names = FALSE, sep="\t")

###################################################
### Chapter 5.2 Relatedness and IBD Estimation  ###
###################################################

library(SNPRelate)

close(genoData)
gdsobj <- snpgdsOpen(geno.file)
ibdobj <- snpgdsIBDKING(gdsobj, sample.id=NULL, snp.id=NULL)

# snpgdsClose(gdsobj)
# names(ibdobj)
# dim(ibdobj$kinship)
# ibdobj$kinship[1:5,1:5]

###################################################
### Select samples with 4 degrees of kindship
###################################################
samp <- pData(scanAnnot)[,c("scanID", "scanName")]
samp <- samp[match(ibdobj$sample.id, samp$scanID),]
names(samp) <- c("scanID", "Individ")
ibd <- snpgdsIBDSelection(ibdobj, kinship.cutoff=1/32)
ibd <- merge(ibd, samp, by.x="ID1", by.y="scanID")
ibd <- merge(ibd, samp, by.x="ID2", by.y="scanID", suffixes=c("1","2"))
ibd$ii <- pasteSorted(ibd$Individ1, ibd$Individ2)

ibd$exp.rel <- "U"
ibd$exp.rel[ibd$Individ1 == ibd$Individ2] <- "Dup"

# table(ibd$exp.rel, useNA="ifany")
# assign observed relationships
ibd$obs.rel <- ibdAssignRelatednessKing(ibd$IBS0, ibd$kinship)
table(ibd$obs.rel, useNA="ifany")

table(ibd$exp.rel, ibd$obs.rel, useNA="ifany")

###################################################
### Plot the estimated kinship  
###################################################
## thresholds for assigning relationships using kinship coefficients
## in table 1 of Manichaikul (2010)
cut.dup <- 1/(2^(3/2))
cut.deg1 <- 1/(2^(5/2))
cut.deg2 <- 1/(2^(7/2))
cut.deg3 <- 1/(2^(9/2))
cols <- c(Dup="magenta", U="black")

# ibd2 <- subset(ibd, ID2%nin%scan.exclude | ID1%nin%scan.exclude)
ibd2 <- ibd

pdf("./QC_check/Ch5.2-1 Kinship Check.pdf")
plot(ibd2$IBS0, ibd2$kinship, col=cols[ibd2$exp.rel],
     xlab="Fraction of IBS=0", ylab="Kinship coefficient")
abline(h=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), lty=2, col="gray")
legend("topright", legend=names(cols), col=cols, pch=1)
dev.off()

write.table(ibd2, "GEM_GWAS_Kinship_Estimate.txt", quote=FALSE, row.names = FALSE, sep="\t")

###################################################
### Chapter 5.3 PCA analysis 
###################################################

filt <- get(data(pcaSnpFilters.hg19))
chrom <- getChromosome(snpAnnot)
pos <- getPosition(snpAnnot)
snpID <- getSnpID(snpAnnot)
snp.filt <- rep(TRUE, length(snpID))
for (f in 1:nrow(filt)) {
  snp.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos
           & pos < filt$end.base[f]] <- FALSE
}
snp.sel <- snpID[snp.filt]
length(snp.sel)
# 1859042

sample.sel <- scanAnnot$scanID[(scanAnnot$scanID%nin%scan.exclude & scanAnnot$duplicated==FALSE)]
length(sample.sel)
# 248

# gdsobj <- snpgdsOpen(geno.file)
snpset <- snpgdsLDpruning(gdsobj, sample.id=sample.sel, snp.id=snp.sel,
                          autosome.only=TRUE, maf=0.05, missing.rate=0.05,
                          method="corr", slide.max.bp=10e6,
                          ld.threshold=sqrt(0.1))

snp.pruned <- unlist(snpset, use.names=FALSE)
length(snp.pruned)

###################################################
### Calculate first 10 PCs  
###################################################

pca <- snpgdsPCA(gdsobj, sample.id=sample.sel, snp.id=snp.pruned)
names(pca)
length(pca$eigenval)
dim(pca$eigenvect)

###################################################
### Calculate the percentage of variance explained
### by first four principal component.
###################################################

pc.frac <- pca$eigenval/sum(pca$eigenval, na.rm=TRUE)
lbls <- paste("EV", 1:4, "\n", format(pc.frac[1:4], digits=2), sep="")
samp <- pData(scanAnnot)[match(pca$sample.id, scanAnnot$scanID),]
cols <- rep(NA, nrow(samp))

samp$re_race <- ifelse(samp$race == "W", "White", 
                       ifelse(samp$race == "B", "Black", 
                              ifelse(samp$race == "A", "Asian", 
                                     ifelse(samp$race %in% c("-11", "-4", "0"), "Unknown", "Other"))))
# Assign colors for each race
cols[samp$re_race == "White"] <- "blue"
cols[samp$re_race == "Black"] <- "red"
cols[samp$re_race == "Other"] <- "black"
cols[samp$re_race == "Asian"] <- "green"
cols[samp$re_race == "Unknown"] <- "grey"

###################################################
### Plot PCA pairs 
###################################################

pdf("./QC_check/Ch5.3-1 PCA Pairs.pdf")
pairs(pca$eigenvect[,1:4], col=cols, labels=lbls,
      main = "White: blue, Black: red, Asian: green, Others: black, Unknown: grey")
dev.off()

pairs(pca$eigenvect[,1:2], col=cols, 
      main = "White: blue, Black: red, Asian: green, Others: black, Unknown: grey")

###################################################
### Plot the PCA Parallel Coorinates
###################################################

par.coord <- pca$eigenvect
rangel <- apply(par.coord, 2, function(x) range(x)[1])
rangeh <- apply(par.coord, 2, function(x) range(x)[2])
std.coord <- par.coord
for (i in 1:14)
  std.coord[,i] <- (par.coord[,i] - rangel[i])/(rangeh[i]-rangel[i])

pdf("./QC_check/Ch5.3-2 PCA Parallel Coorinates.pdf")
plot(c(0,15), c(0,1), type = 'n', axes = FALSE, ylab = "", xlab = "",
     main = "Parallel Coordinates Plot
     White: blue, Black: red, Asian: green, Others: black, Unknown: grey")
for (j in 1:13)
  for (i in sample(1:nrow(std.coord)) )
    lines(c(j,j+1), std.coord[i,c(j,j+1)], col=cols[i], lwd=0.25)
axis(1, at = 1:14, labels = paste("PC",1:14, sep = "."))
dev.off()

###################################################
### Plot the PCA correlations
###################################################

corr <- snpgdsPCACorr(pca, gdsobj, eig.which=1:4)
snpgdsClose(gdsobj)

snp <- snpAnnot[match(corr$snp.id, snpID),]
chrom <- getChromosome(snp, char=TRUE)
pdf("./QC_check/Ch5.3-3 DataCleaning-corr.pdf")
par(mfrow=c(4,1))
for (i in 1:4) {
  snpCorrelationPlot(abs(corr$snpcorr[i,]), chrom,
                     main=paste("Eigenvector",i), ylim=c(0,1))
}
dev.off()

###--------------------------------------------------------------------------###
### (5) Case-Control Confounding -- Chapter 6
###--------------------------------------------------------------------------###

###################################################
### Ch6.2 Missing call rate difference
###################################################

# Collect eigen vectors from the PCS analysis
princomp <- as.data.frame(pca$eigenvect)

# scan Annotation from non-duplicated samples
samples.nodup <- pData(scanAnnot)[sample.sel,]

# combine scanID, case/control, and race information
princomp$scanID <- as.factor(samples.nodup$scanID)
princomp$case.ctrl.status <- as.factor(samples.nodup$Case_Control)
princomp$race <- as.factor(samples.nodup$race)

# Calculate the % variations explained by each PC
pc.percent <- 100 * pca$eigenval[1:32]/sum(pca$eigenval, na.rm = TRUE)
pc.percent

# Make a lable for plot later
lbls <- paste("EV", 1:3, "\n", format(pc.percent[1:3], digits=2), "%", sep="")

# Check numbers of AD and Controls in the Data
table(samples.nodup$Case_Control)

# Assign color codes for AD and Control
cols <- rep(NA, nrow(samples.nodup))
cols[samples.nodup$Case_Control == "Case"] <- "green"
cols[samples.nodup$Case_Control == "Control"] <- "magenta"

# Plot First Three EVs by Case-Control Status
pdf("./QC_check/Ch6.1-1 First Three EVs by Case-Control Status.pdf")
pairs(pca$eigenvect[,1:3], col=cols, labels=lbls,
      main = "First Three EVs by Case-Control Status")
dev.off()

# Plot PC1 vs. Case-control Status
pdf("./QC_check/Ch6.1-2 PC1 vs. Case-control Status.pdf")
boxplot(princomp[, 1] ~ princomp$case.ctrl.status,
        ylab = "PC1", main = "PC1 vs. Case-control Status")
dev.off()

# Plot PC2 vs. Case-control Status
pdf("./QC_check/Ch6.1-3 PC2 vs. Case-control Status.pdf")
boxplot(princomp[, 2] ~ princomp$case.ctrl.status,
        ylab = "PC2", main = "PC2 vs. Case-control Status")
dev.off()

# Plot PC3 vs. Case-control Status
pdf("./QC_check/Ch6.1-4 PC3 vs. Case-control Status.pdf")
boxplot(princomp[, 3] ~ princomp$case.ctrl.status,
        ylab = "PC3", main = "PC3 vs. Case-control Status")
dev.off()

# Anova for first 3 PCs vs Race
# aov.p1 <- aov(princomp[,1] ~ princomp$race *
#                 princomp$case.ctrl.status, princomp)
# summary(aov.p1)
# aov.p2 <- aov(princomp[,2] ~ princomp$race *
#                 princomp$case.ctrl.status, princomp)
# summary(aov.p2)
# aov.p3 <- aov(princomp[,3] ~ princomp$race *
#                 princomp$case.ctrl.status, princomp)
# summary(aov.p3)

###################################################
### Ch6.2 Missing call rate difference
###################################################

# Testing the association between missing call rate vs. status
lm.all <- lm(scanAnnot$missing.e1 ~ scanAnnot$Case_Control)
summary(aov(lm.all))

#                     Df  Sum Sq   Mean Sq F value Pr(>F)
# scanAnnot$AD        1 0.000001 1.031e-06   0.126  0.723
# Residuals         235 0.001924 8.185e-06     
# 11 observations deleted due to missingness

# Plot Missing call rate by status
pdf("./QC_check/Ch6.2-1 Missing call rate by status.pdf")
boxplot(scanAnnot$missing.e1 ~ scanAnnot$Case_Control, 
        ylab="Mean missing call rate", main="Mean missing call rate by case status")
dev.off()

###--------------------------------------------------------------------------###
### (6) Chromosome Anomaly -- Chapter 7
###--------------------------------------------------------------------------###

###################################################
### Re-open the gds geno data
###################################################

gds1 <- GdsGenotypeReader(geno.file)
genoData <- GenotypeData(gds1, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
genoData

###################################################
### Ch7.1 B Allele Frequency filtering
###################################################

# Identify some low quality samples by looking at the standard deviation of BAF

baf.sd <- sdByScanChromWindow(blData, genoData, var="BAlleleFreq")
med.baf.sd <- medianSdOverAutosomes(baf.sd)
low.qual.ids <- med.baf.sd$scanID[med.baf.sd$med.sd > 0.05]

length(low.qual.ids)
# 248

###################################################
### Decide which SNPs to exclude based on genome build
###################################################

chrom <- getChromosome(snpAnnot, char=TRUE)
pos <- getPosition(snpAnnot)
hla.df <- get(data(HLA.hg19))
hla <- chrom == "6" & pos >= hla.df$start.base & pos <= hla.df$end.base
xtr.df <- get(data(pseudoautosomal.hg19))
xtr <- chrom == "X" & pos >= xtr.df["X.XTR", "start.base"] &
  pos <= xtr.df["X.XTR", "end.base"]
centromeres <- get(data(centromeres.hg19))
gap <- rep(FALSE, nrow(snpAnnot))
for (i in 1:nrow(centromeres)) {
  ingap <- chrom == centromeres$chrom[i] & pos > centromeres$left.base[i] &
    pos < centromeres$right.base[i]
  gap <- gap | ingap
}
# ignore includes intensity-only and failed snps
ignore <- snpAnnot$missing.n1 == 1 
snp.exclude <- hla | xtr | gap | ignore
snp.ok <- snpAnnot$snpID[!snp.exclude]

length(snp.ok)
# 1,878,117

###################################################
### circular binary segmentation to find change points in BAF.
###################################################
scan.ids <- scanAnnot$scanID[1:10]
chrom.ids <- 21:23
baf.seg <- anomSegmentBAF(blData, genoData, scan.ids=scan.ids,
                          chrom.ids=chrom.ids, snp.ids=snp.ok, verbose=FALSE)
head(baf.seg)

str(baf.seg)
# 35 obs

###################################################
### Filter segments to detect anomalies, treating the low quality samples differently
###################################################

baf.anom <- anomFilterBAF(blData, genoData, segments=baf.seg,
                          snp.ids=snp.ok, centromere=centromeres, low.qual.ids=low.qual.ids,
                          verbose=FALSE)
names(baf.anom)
baf.filt <- baf.anom$filtered
head(baf.filt)
# NULL

###################################################
### Ch7.2 Loss of Heterozygosity
###################################################

# look for Loss of Heterozygosity (LOH) anomalies by identifying homozygous
# runs with change in LRR. 
# Change points in LRR are found by circular binary segmentation. 
# Known anomalies from the BAF detection are excluded.

loh.anom <- anomDetectLOH(blData, genoData, scan.ids=scan.ids,
                          chrom.ids=chrom.ids, snp.ids=snp.ok, known.anoms=baf.filt,
                          verbose=FALSE)
names(loh.anom)
loh.filt <- loh.anom$filtered
str(loh.filt)
# 21 obs 

###################################################
### Ch7.3 Statistics
###################################################

# Calculate statistics for the anomalous segments found with the BAF and LOH methods.
# create required data frame

# baf.filt$method <- "BAF"
# Since baf.filt has no output

if (!is.null(loh.filt)) {
  #   loh.filt$method <- "LOH"
  cols <- intersect(names(baf.filt), names(loh.filt))
  anoms <- rbind(baf.filt[,cols], loh.filt[,cols])
} else {
  anoms <- baf.filt
}
anoms$anom.id <- 1:nrow(anoms)
anoms$method <- "LOH"

stats <- anomSegStats(blData, genoData, snp.ids=snp.ok, anom=anoms,
                      centromere=centromeres)
# names(stats)

###################################################
### Plot the anomalies with relevant statistics, one anomaly per plot. 
### Each plot has two parts:
###   upper part is a graph of LRR and lower part is a graph of BAF
###################################################

snp.not.ok <- snpAnnot$snpID[snp.exclude]
pdf("./QC_check/Ch7.3-1 Anomalies.pdf")
for(i in 1:nrow(stats)){
  anomStatsPlot(blData, genoData, anom.stats=stats[i,],
                snp.ineligible=snp.not.ok, centromere=centromeres, cex.leg=1)
}
dev.off()

###################################################
### Ch7.4 Identify low quality samples
###################################################

# To identify low quality samples, one measure we use is the 
# standard deviation of BAF and LRR.
# BAF results were found previously, 
# now we find results for LRR. Unlike for BAF, all genotypes are included.
lrr.sd <- sdByScanChromWindow(blData, var="LogRRatio", incl.hom=TRUE)
med.lrr.sd <- medianSdOverAutosomes(lrr.sd)

###################################################
### Find the number of segments found using circular binary segmentation 
### in anomaly detection
###################################################

baf.seg.info <- baf.anom$seg.info
loh.seg.info <- loh.anom$base.info[,c("scanID", "chromosome", "num.segs")]

###################################################
### identify low quality samples separately for BAF and LOH, 
### using different threshold parameters
###################################################
snpAnnot$eligible <- !snp.exclude
baf.low.qual <- anomIdentifyLowQuality(snpAnnot, med.baf.sd, baf.seg.info,
                                       sd.thresh=0.1, sng.seg.thresh=0.0008, auto.seg.thresh=0.0001)
loh.low.qual <- anomIdentifyLowQuality(snpAnnot, med.lrr.sd, loh.seg.info,
                                       sd.thresh=0.25, sng.seg.thresh=0.0048, auto.seg.thresh=0.0006)

length(baf.low.qual)
# 10
length(loh.low.qual)


###################################################
### Ch7.5 Filter anomalies
###################################################
# anomalies to filter
anom.filt <- stats[,c("scanID", "chromosome", "left.base", "right.base")]
# whole.chrom column is required and can be used for sex chromosome
#   anomalies such as XXX
anom.filt$whole.chrom <- FALSE
# select unique subjects
subj <- scanAnnot$scanID[!scanAnnot$duplicated]
subj.filt.file <- "subj_filt.gds"

# close(gds1)
setMissingGenotypes(geno.file, subj.filt.file, anom.filt,
                    file.type="gds", sample.include=subj, verbose=FALSE)

# (gds <- GdsGenotypeReader(subj.filt.file))
# close(gds)

###--------------------------------------------------------------------------###
### (7) SNP Quality Checks -- Chapter 8
###--------------------------------------------------------------------------###

###################################################
### Ch8.1 Duplicate Sample Discordance
###################################################

# Find the bad samples (missing > 0.05) & bad snps (missing = 100%)
scan.excl <- scanAnnot$scanID[scanAnnot$missing.e1 >= 0.05]
length(scan.excl)
# 0
snp.excl <- snpAnnot$snpID[snpAnnot$missing.n1 == 1]
length(snp.excl)
# 9935 

# Reload the genotype data
# gds1 <- GdsGenotypeReader(geno.file)
# genoData <- GenotypeData(gds1, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
# genoData

# Calcualte the discordance 
# do not have duplicate samples, skip this chunk
# dupdisc <- duplicateDiscordance(genoData, subjName.col="scanName",
#                                 scan.exclude=scan.excl, 
#                                 snp.exclude=snp.excl)

# Check the results
# names(dupdisc)
# head(dupdisc$discordance.by.snp)
# length(dupdisc$discordance.by.subject)
# dupdisc$discordance.by.subject[[2]]
# each entry is a 2x2 matrix, but only one value of each
# is important since these are all pairs

# npair <- length(dupdisc$discordance.by.subject)
# disc.subj <- rep(NA, npair)
# subjID <- rep(NA, npair)
# race <- rep(NA, npair)
# for (i in 1:npair) {
#   disc.subj[i] <- dupdisc$discordance.by.subject[[i]][1,2]
#   subjID[i] <- names(dupdisc$discordance.by.subject)[i]
#   race[i] <- scanAnnot$race[scanAnnot$GEMID == subjID[i]][1]
# }
# dat <- data.frame(subjID=subjID, disc=disc.subj, pop=race,
#                   stringsAsFactors=FALSE)

# summary(dat$disc)

# Assign colors for the duplicate samples based on population group.
# dat$col <- NA
# dat$col[dat$pop == "Black"] <- "blue"
# dat$col[dat$pop == "White"] <- "red"
# dat$col[dat$pop == "Asian"] <- "green"
# dat$col[dat$pop == "Others"] <- "purple"

# dat <- dat[order(dat$disc),]
# dat$rank <- 1:npair

# Plot the sample discordance rates color coded by race.

# pdf("./QC_check/Ch8.1-1 Sample discordance rates.pdf")
# plot(dat$disc, dat$rank, col=dat$col, ylab="rank", las=1,
#      xlab="Discordance rate between duplicate samples",
#      main="Duplicate Sample Discordance by Continental Ancestry")
# legend("bottomright", unique(dat$pop), pch=rep(1,2), col=unique(dat$col))
# dev.off()

###################################################
### 
###################################################
# duplicateDiscordanceProbability(npair)

###################################################
### Ch8.2 Mendelian Error Checking 
###
###       Not appliable in our study !!!
###################################################

# men.list <- with(pData(scanAnnot), mendelList(family, subjectID,
#                                               father, mother, sex, scanID))
# res <- mendelListAsDataFrame(men.list)
# head(res)
# dim(res)
# # Only want to use SNPs with missing.n1 < 0.05
# snp.excl <- snpAnnot$snpID[snpAnnot$missing.n1 >= 0.05]
# length(snp.excl)
# mend <- mendelErr(genoData, men.list, snp.exclude=snp.excl)
# names(mend)
# head(mend$trios)
# names(mend$snp)


# ###################################################
# ### code chunk number 94: DataCleaning.Rnw:2102-2105
# ###################################################
# # Calculate the error rate
# err <- mend$snp$error.cnt / mend$snp$check.cnt
# table(err == 0, useNA="ifany")


# ###################################################
# ### code chunk number 95: DataCleaning.Rnw:2108-2110
# ###################################################
# plot(err, rank(err), xlab="Error Rate (fraction)",
#      ylab="rank", main="Mendelian Error Rate per SNP, ranked")


# ###################################################
# ### code chunk number 96: DataCleaning.Rnw:2118-2130
# ###################################################
# fam <- mend$snp$error.cnt
# n <- mend$snp$check.cnt
# summary(fam)
# # SNPs with errors
# length(fam[n > 0 & fam > 0])
# # SNPs for which more than one family has an error
# length(fam[n > 0 & fam > 1])
# # Get the SNPs with valid trios for error detection
# val <- length(fam[n > 0])
# noerr <- length(fam[n > 0 & fam == 0])
# # Divide to get fraction with no errors
# noerr / val


# ###################################################
# ### code chunk number 97: DataCleaning.Rnw:2137-2150
# ###################################################
# snp.sel <- match(names(mend$snp$error.cnt), snpAnnot$snpID)
# snpAnnot$mendel.err.count[snp.sel] <- mend$snp$error.cnt
# snpAnnot$mendel.err.sampsize[snp.sel] <- mend$snp$check.cnt
# allequal(snpAnnot$snpID, sort(snpAnnot$snpID))
# # The high number of NA values is due to the filtering out of SNPs
# #    before the Mendelian error rate calculation
# sum(is.na(snpAnnot$mendel.err.count))
# sum(is.na(snpAnnot$mendel.err.sampsize))
# varMetadata(snpAnnot)["mendel.err.count", "labelDescription"] <-
#   paste("number of Mendelian errors detected in trios averaged over",
#         "multiple combinations of replicate genotyping instances")
# varMetadata(snpAnnot)["mendel.err.sampsize", "labelDescription"] <-
#   "number of opportunities to detect Mendelian error in trios"


# ###################################################
# ### code chunk number 98: DataCleaning.Rnw:2162-2179
# ###################################################
# # Get a vector of SNPs to check
# snp <- pData(snpAnnot)
# snp$err.rate <- snp$mendel.err.count /
#   snp$mendel.err.sampsize
# snp <- snp[order(snp$err.rate, decreasing=TRUE),]
# snp <- snp[1:9,]
# xyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
# xyGDS <- GdsIntensityReader(xyfile)
# xyData <- IntensityData(xyGDS, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
# pdf(file="DataCleaning-mendel.pdf")
# par(mfrow = c(3,3))
# mtxt <- paste("SNP", snp$rsID, "\nMendelian Error Rate",
#               format(snp$err.rate, digits=5))
# genoClusterPlot(xyData, genoData, snpID=snp$snpID, main.txt=mtxt,
#                 cex.main=0.9)
# dev.off()
# close(xyData)


# ###################################################
# ### code chunk number 99: DataCleaning.Rnw:2192-2209
# ###################################################
# # Calculate the fraction of SNPs with an error for each trio
# trios <- mend$trios
# trios$Mend.err <- trios$Men.err.cnt/trios$Men.cnt
# summary(trios$Mend.err)
# # Start by pulling out the vectors needed from `trios'
# tmp <- trios[, c("fam.id", "Mend.err")]; dim(tmp)
# # Change fam.id to match the sample annotation column name
# names(tmp) <- c("family", "Mend.err.rate.fam")
# # Merge the variables into the sample annotation file
# scanAnnot$mend.err.rate.fam <- NA
# for (i in 1:nrow(tmp)) {
#   ind <- which(is.element(scanAnnot$family, tmp$family[i]))
#   scanAnnot$mend.err.rate.fam[ind] <- tmp$Mend.err.rate.fam[i]
# }
# head(scanAnnot$mend.err.rate.fam)
# varMetadata(scanAnnot)["mend.err.rate.fam", "labelDescription"] <-
#   "Mendelian error rate per family"


# ###################################################
# ### code chunk number 100: DataCleaning.Rnw:2217-2228
# ###################################################
# # Get the families that have non-NA values for the family
# #     Mendelian error rate
# fams <- pData(scanAnnot)[!is.na(scanAnnot$mend.err.rate.fam) &
#                            !duplicated(scanAnnot$family), c("family",
#                                                             "mend.err.rate.fam", "race")]
# dim(fams)
# table(fams$race, useNA="ifany")
# # Assign colors for the different ethnicities in these families
# pcol <- rep(NA, nrow(fams))
# pcol[fams$race == "CEU"] <- "blue"
# pcol[fams$race == "YRI"] <- "red"


# ###################################################
# ### code chunk number 101: DataCleaning.Rnw:2231-2236
# ###################################################
# plot(fams$mend.err.rate.fam*100, rank(fams$mend.err.rate.fam),
#      main="Mendelian Error rate per Family, ranked",
#      xlab="Mendelian error rate per family (percent)",
#      ylab="rank", col=pcol)
# legend("bottomright", c("CEU", "YRI"), pch=c(1,1), col=c("blue", "red"))


# ###################################################
# ### code chunk number 102: DataCleaning.Rnw:2254-2269
# ###################################################
# head(pData(scanAnnot)[,c("father", "mother")])
# nonfounders <- scanAnnot$father != 0 &
#   scanAnnot$mother != 0
# table(nonfounders)
# scan.excl <- scanAnnot$scanID[scanAnnot$race != "CEU" |
#                                 nonfounders | scanAnnot$duplicated]
# length(scan.excl)

# chr <- getChromosome(genoData)
# auto <- range(which(chr %in% 1:22))
# X <- range(which(chr == 23))
# hwe <- exactHWE(genoData, scan.exclude=scan.excl, snpStart=auto[1], snpEnd=auto[2])
# hweX <- exactHWE(genoData, scan.exclude=scan.excl, snpStart=X[1], snpEnd=X[2])
# hwe <- rbind(hwe, hweX)
# close(genoData)


# ###################################################
# ### code chunk number 103: DataCleaning.Rnw:2275-2285
# ###################################################
# names(hwe)
# dim(hwe)
# # Check on sample sizes for autosomes and X chromosome
# hwe$N <- hwe$nAA + hwe$nAB + hwe$nBB
# summary(hwe$N[is.element(hwe$chr,1:22)])
# summary(hwe$N[is.element(hwe$chr,23)])
# hwe$pval[1:10]
# sum(is.na(hwe$pval[hwe$chr == 23])) # X
# hwe$MAF[1:10]
# hwe[1:3, c("nAA", "nAB", "nBB")]


# ###################################################
# ### code chunk number 104: DataCleaning.Rnw:2292-2293
# ###################################################
# summary(hwe$f)


# ###################################################
# ### code chunk number 105: DataCleaning.Rnw:2296-2298
# ###################################################
# hist(hwe$f, main="Histogram of the Inbreeding Coefficient
#   For CEU Samples", xlab="Inbreeding Coefficient")


# ###################################################
# ### code chunk number 106: DataCleaning.Rnw:2301-2304
# ###################################################
# # Check the MAF of those SNPs with f=1
# chkf <- hwe[!is.na(hwe$f) & hwe$f==1,]; dim(chkf)
# summary(chkf$MAF)


# ###################################################
# ### code chunk number 107: DataCleaning.Rnw:2311-2322
# ###################################################
# hwe.0 <- hwe[hwe$MAF > 0,]; dim(hwe.0)
# # Only keep the autosomal SNPs for first plot
# pval <- hwe.0$pval[is.element(hwe.0$chr, 1:22)]
# length(pval)
# pval <- pval[!is.na(pval)]
# length(pval)
# # X chromosome SNPs for plot 2
# pval.x <- hwe.0$pval[is.element(hwe.0$chr, 23)]
# length(pval.x)
# pval.x <- pval.x[!is.na(pval.x)]
# length(pval.x)


# ###################################################
# ### code chunk number 108: DataCleaning.Rnw:2324-2331
# ###################################################
# pdf(file = "DataCleaning-hwe.pdf")
# par(mfrow=c(2,2))
# qqPlot(pval=pval, truncate = FALSE, main="Autosomes, all")
# qqPlot(pval=pval, truncate = TRUE, main="Autosomes, truncated")
# qqPlot(pval=pval.x, truncate = FALSE, main="X chromosome, all")
# qqPlot(pval=pval.x, truncate = TRUE, main="X chromosome, truncated")
# dev.off()


# ###################################################
# ### code chunk number 109: DataCleaning.Rnw:2337-2340
# ###################################################
# plot(hwe.0$MAF, -log10(hwe.0$pval),
#      xlab="Minor Allele Frequency", ylab="-log(p-value)",
#      main="Minor Allele Frequency vs\nP-value")# 

###--------------------------------------------------------------------------###
### (8) Preliminary Association Tests -- Chapter 8
###--------------------------------------------------------------------------###

# ###################################################
# ### code chunk number 110: DataCleaning.Rnw:2382-2398
# ###################################################
# genoGDS <- GdsGenotypeReader(subj.filt.file)
# subjAnnot <- scanAnnot[scanAnnot$scanID %in% getScanID(genoGDS),]
# subjAnnot$sex <- as.factor(subjAnnot$sex)
# subjAnnot$EV1 <- pca$eigenvect[match(subjAnnot$scanID, pca$sample.id), 1]
# 
# genoData <- GenotypeData(genoGDS, scanAnnot=subjAnnot)
# chr <- getChromosome(genoData)
# assoc.list <- lapply(unique(chr), function(x) {
#   ## Y chromsome only includes males, cannot have sex as a covariate
#   covar <- ifelse(x == 25, "EV1", c("sex", "EV1"))
#   start <- which(chr == x)[1]
#   assocRegression(genoData, outcome="status", covar=covar, model.type="logistic", 
#                   snpStart=start, snpEnd=start+50)
# })
# assoc <- do.call(rbind, assoc.list)
# close(genoData)
# 
# 
# ###################################################
# ### code chunk number 111: DataCleaning.Rnw:2412-2414
# ###################################################
# qqPlot(pval=assoc$Wald.pval,
#        truncate=TRUE, main="QQ Plot of Wald Test p-values")
# 
# 
# ###################################################
# ### code chunk number 112: DataCleaning.Rnw:2422-2425
# ###################################################
# chrom <- getChromosome(snpAnnot, char=TRUE)
# snp.sel <- getSnpID(snpAnnot) %in% assoc$snpID
# manhattanPlot(assoc$Wald.pval, chromosome=chrom[snp.sel])
# 
# 
# ###################################################
# ### code chunk number 113: DataCleaning.Rnw:2436-2454
# ###################################################
# # Identify SNPs with lowest p-values
# snp <- pData(snpAnnot)[snp.sel, c("snpID", "rsID")]
# snp$pval <- assoc$Wald.pval
# snp <- snp[order(snp$pval),]
# snp <- snp[1:9,]
# xyfile <- system.file("extdata", "illumina_qxy.gds", package="GWASdata")
# xyGDS <- GdsIntensityReader(xyfile)
# xyData <- IntensityData(xyGDS, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
# genofile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
# genoGDS <- GdsGenotypeReader(genofile)
# genoData <- GenotypeData(genoGDS, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
# pdf(file="DataCleaning-cluster.pdf")
# par(mfrow = c(3,3))
# mtxt <- paste("SNP", snp$rsID, "\np =", format(snp$pval, digits=4))
# genoClusterPlot(xyData, genoData, snpID=snp$snpID, main.txt=mtxt)
# dev.off()
# close(xyData)
# close(genoData)
# 
# 
# ###################################################
# ### code chunk number 114: DataCleaning.Rnw:2461-2462
# ###################################################
# unlink(subj.filt.file)


###--------------------------------------------------------------------------###
### (9) Select eligible samples and snps and make the vcf file for imputation
###--------------------------------------------------------------------------###

elig_subj <- scanAnnot$scanID[!scanAnnot$duplicated & !(scanAnnot$scanName %in% c("IGN598", "IGN749", "IGN751")) & scanAnnot$missing.e1 < 0.05 ]
elig_snps <- snpAnnot$snpID[!snp.exclude]

eligible_subj_snps_file <- "eligible_subj_snps.gds"

# showfile.gds(closeall=TRUE)
# Reload the genotype data
# geno.file <- "GWAS5_geno.gds"
# gds1 <- GdsGenotypeReader(geno.file)
# genoData <- GenotypeData(gds1, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
# genoData

gdsSubset(geno.file, eligible_subj_snps_file,
          sample.include=elig_subj, snp.include=elig_snps,
          sub.storage=NULL,
          compress="LZMA_RA",
          block.size=5000,
          verbose=TRUE)

# 
genoGDS <- GdsGenotypeReader(eligible_subj_snps_file)
subjAnnot <- scanAnnot[scanAnnot$scanID %in% getScanID(genoGDS),]
snpsAnnot <- snpAnnot[snpAnnot$snpID %in% getSnpID(genoGDS),]

genoData <- GenotypeData(genoGDS, snpAnnot=snpsAnnot, scanAnnot=subjAnnot)
genoData

# update ref and alt allele information from reference file
ref.allele=snpsAnnot@data$a1
alt.allele=snpsAnnot@data$a2

vcfWrite(genoData, vcf.file="GWAS_5_raw.vcf", sample.col="scanName",
         id.col="snpName", qual.col=NULL, filter.cols=NULL,
         info.cols=NULL, scan.exclude=NULL, snp.exclude=NULL,
         scan.order=NULL, ref.allele=snpsAnnot@data$a1, block.size=1000, verbose=TRUE)

save.image("./QC_check/GWAS_5_QC.RData")
# load("./QC_check/GWAS_5_QC.RData")

# d1 <- merge(ibd[,c("Individ1","ii","exp.rel","obs.rel")], 
#             scanAnnot2[c("scanName","sex", "race", "Case_Control", "Disease_Status", "APOE", "plate","well")], 
#             by.x = "Individ1", by.y = "scanName")

# colnames(d1) <- c("Individ1", "ii", "exp.rel", "obs.rel", paste0(colnames(d1)[5:11], "_Ind1"))

# d2 <- merge(ibd[,c("Individ2","ii","exp.rel","obs.rel")], 
#             scanAnnot2[c("scanName","sex", "race", "Case_Control", "Disease_Status", "APOE", "plate","well")], 
#             by.x = "Individ2", by.y = "scanName")

# colnames(d2) <- c("Individ2", "ii", "exp.rel", "obs.rel", paste0(colnames(d2)[5:11], "_Ind2"))

# d <- merge(d1, d2, by = intersect(colnames(d1), colnames(d2)))
# View(d)
# write.table(d, "d.csv", 
#             col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")
