
###--------------------------------------------------------------------------###
### environment settings
###--------------------------------------------------------------------------###

rm(list=ls())
options(stringsAsFactors = FALSE)

# set the working directory to ***your own*** directory 
setwd("/zfs1/kfan/Ruyu/harmonization_Sep5/GEM/")

# laod the packages that we are going to use later
# This is the main package we used in GWAS QC
library(GWASTools)
# This package can load Excel (.xlsx) files into R 
library(gdata)

# a small function to select which "x" is not in "y"
# which will be using lot of times later
"%nin%" <- function(x,y){
  !x%in%y
}

###--------------------------------------------------------------------------###
### The following R script is creating a scan annotation DF required by 
###   createDataFile funciton in the GWASTools package
### The following comments are copied from the createDataFile funciton 
###   in the GWASTools namual 
###       scan.annotation
###			      Scan annotation data.frame with columns "scanID" (unique id of 
###           genotyping instance), "scanName", (sample name inside the raw 
###           data file) and "file" (corresponding raw data file name).
###
### We might need the sex for quality check
###   we need a phenotype file for other information
###--------------------------------------------------------------------------###

# read-in the scanName
SN <- read.table("GEM_scanName.txt", header=FALSE)

# create a text string with all the individual file names
SF <- paste("GEM_Sample_",SN$V1,".txt", sep="")

# make a main dataframe to merge other information later
d <- data.frame(scanName=SN$V1, file=SF)

###--------------------------------------------------------------------------###
### Other clinical information such as age, race, and gender
###   as well as the technical information like the batch and the plate
###--------------------------------------------------------------------------###

#############################
### technical informaiton ###
#############################
# create an empty list to store all 30 sheets
plate <- vector("list",30)

# load the data from excel sheets and save it in the same list object
# and also select the information that we need  
for(i in 1:30){
  plate[[i]] <- as.data.frame(readxl::read_excel("../GEM_Reference_Data/Kamboh GEM-GWAS 2019-01-07.xlsx", sheet=i)[,1:3])
}

# collapse the list in to a single data frame & rename the column
tech <- do.call("rbind", plate)
names(tech) <- c("plate", "well_position", "scanName")

# remove the last one row
# tail(tech)
tech <- tech[1:2832,]

# Remove dulpicate ID (The two IDs accenditally duplicated, but not marked -R in ID)
# Sample_ID	SentrixPosition_A	Sample_Plate	Sample_Well
# 3200965	R02C01	25	B09
# 3200965	R01C01	30	A01
# 3206165	R08C01	24	H01
# 3206165	R08C01	30	H04
# By checking the scanName file, the earlier duplicate genotyping was removed 
# therefore, but duplicated samples are genpotyped on plate 30 !!!

ind1 <- which(tech$scanName==3200965 & tech$plate==25)
ind2 <- which(tech$scanName==3206165 & tech$plate==24)

tech <- tech[-c(ind1,ind2),]

# str(tech)

########################################################
### merge the clinical and the technical information ###
########################################################

sub_info <- merge(d, tech, all.x=TRUE, by="scanName")

# Transfer some Cbal ID back to GEM ID
sub_info$scanName[sub_info$scanName=="1193"] <- "6208185"
sub_info$scanName[sub_info$scanName=="1254"] <- "6208959"
sub_info$scanName[sub_info$scanName=="2297"] <- "5206154"
sub_info$scanName[sub_info$scanName=="2532"] <- "4210301"
sub_info$scanName[sub_info$scanName=="600"] <- "4201159"

# Create a corrected GEM ID 
sub_info$plate <- as.integer(sub_info$plate)
sub_info$corrected_scanName <- sub_info$scanName

# str(sub_info)

tech[which(tech$scanName%nin%d$scanName),]
# plate well_position scanName
# 2736    29           H06  6209831
# This sample is failed to genotyping

###--------------------------------------------------------------------------###
### Swap phenotype data between plate3 and plate 4
###--------------------------------------------------------------------------###

# create a mapping dataframe
swapd <- subset(sub_info, plate%in%3:4)

# odered rows by plate and well position
swapd <- swapd[order(swapd$plate, swapd$well_position),]

# swap the plate number
swapd$old_plate <- swapd$plate
swapd$plate <- ifelse(swapd$plate==3,4,3)
# table(swapd$plate, swapd$old_plate)

# swap the scan name according the well position

oldplate3 <- subset(swapd, old_plate==3)
oldplate4 <- subset(swapd, old_plate==4)

for(i in 1:nrow(oldplate3)){
  
  ind <- which(oldplate4$well_position==oldplate3$well_position[i])
  if(length(ind)==1){
    oldplate3$corrected_scanName[i] <- oldplate4$scanName[ind]
  }
}

for(i in 1:nrow(oldplate4)){
  
  ind <- which(oldplate3$well_position==oldplate4$well_position[i])
  if(length(ind)==1){
    oldplate4$corrected_scanName[i] <- oldplate3$scanName[ind]
  }
}

swapd <- rbind(oldplate3, oldplate4)
save(swapd, file="Plate_3_4_Mapping_Data.RData")

new_plate34 <- subset(swapd, select=c("scanName", "file", "plate","well_position","corrected_scanName"))
d2 <- subset(sub_info, plate%nin%3:4)

# str(new_plate34)
# str(d2)
# str(swapd)

d3 <- rbind(d2, new_plate34)
# str(d3)
# test <- d3[d3$plate%in%3:4,]
# test

# Remove Raj's Samples
# Regular samples
raj <- as.data.frame(readxl::read_excel("../GEM_Reference_Data/Raj Sample List.xlsx", sheet=1, skip=1))

# There is one repeat samples 
# 11430-R
# However, this sample doesn't have regular ID in the GWAS
# Therefore, it's not a repeat!!!
# Remove all Raj's Samples
d4 <- d3[d3$corrected_scanName%nin%c(raj$ID,"11430-R"),]

# str(d4)

############################
### clinical information ###
############################
# load the RData from the Dr. Beth Snitz's files
load("../GEM_Reference_Data/Basic_Clinical_Info_for_GEM.RData")

# rename the ID column for merging
names(GEM_CLIN)[1] <- "corrected_scanName"

# merge the clinical data and the experiment data
dd <- merge(d4, GEM_CLIN, by="corrected_scanName", all.x=TRUE)

# get the replicated IDs (-R)
ind <- grep("-R", dd$corrected_scanName)

# overwrite clinical informations if available
for(i in ind){
  rep_id <- dd$corrected_scanName[i]
  nrid <- sub("-R","", rep_id)
  nrclin <- dd[dd$corrected_scanName==nrid,]
  if(nrow(nrclin)==1){
    dd[i,c("sex","race","age","AD","dementia","edu",
           "center","cdr","weight","height","bmi")] <- nrclin[,c("sex","race","age","AD","dementia","edu",
                                                                 "center","cdr","weight","height","bmi")]
  }
}

# Double check
# rep_ids <- grep("-R", dd$corrected_scanName, value=TRUE)
# ind2 <- which(dd$corrected_scanName%in%gsub("-R","",rep_ids))
# ind3 <- sort(c(ind,ind2))
# dd[ind3,]
# looks good ! 

###--------------------------------------------------------------------------###
### Create a Scan Annotation Data Frame
###--------------------------------------------------------------------------###

# Re-order the rows
dd2 <- dd[order(dd$plate, dd$well_position),]

# Since we remove some samples, I would like to re-assign scanID in to continuous integer 
# because it also representing as row number 
dd2$scanID <- 1:nrow(dd2)

# Save the GEM ID and remove the -R 
dd2$GEMID <- as.integer(gsub("-R", "", dd2$corrected_scanName))

# convert GEM ID back to Cbal ID 
#   The GWASTools will match the sample ID in the raw text files
#   each sample file is named as the ID show in the GWAS data
dd2$scanName[dd2$scanName=="6208959"] <- "1254"
dd2$scanName[dd2$scanName=="5206154"] <- "2297"
dd2$scanName[dd2$scanName=="6208185"] <- "1193"
dd2$scanName[dd2$scanName=="4210301"] <- "2532"
dd2$scanName[dd2$scanName=="4201159"] <- "600"

# Re-order the columns
dd2 <- subset(dd2, select=c("scanID","scanName","file","sex","race","age","AD","dementia","plate","well_position",
                            "edu","center","cdr","weight","height","bmi","GEMID","corrected_scanName"))

write.table(as.data.frame(dd2$scanName), "GEM_scanName_update.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# un-named the rows
row.names(dd2) <- NULL

# Create a ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(dd2)

# names of columns
# varLabels(scanAnnot)
# data
# head(pData(scanAnnot))

# Add metadata to describe the columns
meta <- varMetadata(scanAnnot)
meta[c("scanID","scanName","file","sex","race","age","AD","dementia","plate","well_position",
       "edu","center","cdr","weight","height","bmi","GEMID","corrected_scanName"), "labelDescription"] <-
  c("unique ID for scans",
    "subject identifier",
    "raw data file",
    "Sex",
    "Race",
    "Age at the enrollment",
    "AD case or Control",
    "AD Dementia or Non-AD Dementia or Control",
    "Plate number",
    "Well Position on the plate",
    "Education in years",
    "Clinical center",
    "Baseline CDR value",
    "Weight at the first visit",
    "Height at the first visit",
    "BMI at the fist visit",
    "GEMID for subject (may have multiple for replicates)",
    "Actual scan name (may different from scanName)")
varMetadata(scanAnnot) <- meta

# str(scanAnnot)

###--------------------------------------------------------------------------###
### Thr following R script is creating a SNP annotation DF required by 
###   createDataFile funciton in the GWASTools package
### The following comments are copied from the createDataFile funciton 
###   in the GWASTools namual 
###       snp.annotation 
###           Snp annotation dataframe with columns "snpID", "chromosome", "position" and
###           "snpName". snpID should be a unique integer vector, sorted with respect to
###           chromosome and position. snpName should match the snp identifiers inside the
###           raw genoypic data files If file.type="gds", optional columns "alleleA", and
###           "alleleB" will be written if present.
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
### Important Note: 
###     By checking the dbSNP, the position in followed by build 37
###--------------------------------------------------------------------------###

# The reference file provided by Illumina 
# https://support.illumina.com/downloads/infinium-multi-ethnic-global-8-v1-support-files.html
# Infinium Multi-Ethnic Global-8 v1.0 Locus Report 
# multi-ethnic-global-8-d1-locus-report.zip
ref <- read.csv("../GEM_Reference_Data/Multi-EthnicGlobal_D1_LocusReport.txt", sep="\t")

# select the columns we need, discard others
d1 <- subset(ref,select=c("Index","Name", "Chr", "Position"))

# rename the columns 
names(d1) <- c("Index", "snpName", "chromosome", "position")

# make the chromosome index as an integer: (following statement is from "DataCleaning.pdf") 
#   chromosome, an integer mapping for each chromosome, with values 1-27, mapped in order
#   from 1-22, 23=X, 24=XY (the pseudoautosomal region), 25=Y, 26=M (the mitochondrial probes), 
#   and 27=U (probes with unknown positions).                                                                      caution. See the manual pages for more details.)

d1$chromosome[d1$chromosome=="X"] <- 23
d1$chromosome[d1$chromosome=="Y"] <- 25
d1$chromosome[d1$chromosome=="MT"] <- 26
d1$chromosome[d1$chromosome=="0"] <- 27

# If the gene annotation file ("Multi-EthnicGlobal_D1.annotated.txt") is going to be used, 
# then use 23 for XY instead of 24
# d1$chromosome[d1$chromosome=="XY"] <- 23
d1$chromosome[d1$chromosome=="XY"] <- 24

# make the chromosome as a interger
d1$chromosome <- as.integer(d1$chromosome)

###--------------------------------------------------------------------------###
### We are not applying this chunck anymore because it cause some problems!!!
###   (1) It can't deal with flipping allele issues
###   (2) It can only recognize A,T,C,G, and assign all the InDels into missing
###  
### We also need the A/B allele for genotype calling 
### The following comments are copied from the createDataFile funciton 
###   in the GWASTools namual 
###     If the genotype calls are nucleotides (A,C,G,T), the columns "alleleA" and 
###     "alleleB" in snp.annotation are used to map to AB format
###--------------------------------------------------------------------------###

# The reference file provided by Illumina 
# https://support.illumina.com/downloads/infinium-multi-ethnic-global-8-v1-support-files.html
# Infinium Multi-Ethnic Global-8 v1.0 Gene Annotation File (Build 37)
# multi-ethnic-global-8-d1-annotated.zip
# ref_AB <- read.csv("Multi-EthnicGlobal_D1.annotated.txt", sep="\t")

# create two new columns with names alleleA & alleleB 
#-------------------------------------------------#
# checking the format for "Alleles"  
# table(sapply(ref_AB$Alleles,nchar))
# all cell has exactly 5 characters
#-------------------------------------------------#
# ref_AB$alleleA <- sapply(ref_AB$Alleles, substr, start=2, stop=2)
# ref_AB$alleleB <- sapply(ref_AB$Alleles, substr, start=4, stop=4)

# select the columns we need, discard others

# d2 <- subset(ref_AB, select=c("Name","Chr","MapInfo","alleleA", "alleleB"))

# rename the columns

# names(d2) <- c("snpName","chromosome","position","alleleA", "alleleB")

# make the chromosome index as an integer

# d2$chromosome[d2$chromosome=="X"] <- 23
# d2$chromosome[d2$chromosome=="Y"] <- 25
# d2$chromosome[d2$chromosome=="M"] <- 26
# d2$chromosome[d2$chromosome=="0"] <- 27

# d2$chromosome <- as.integer(d2$chromosome)

###--------------------------------------------------------------------------###
### Merge two datasets to create snp annotation file  
###
### merge by snpName, chromosome, and the position to make sure 
### everything is aligned properly
###
### can't not use chromosome because it's XY, MT in one file, but 
### it's X, Y, M in the other. Since we are using Build 37 in both files,
### then we should be able to match variants by position
###--------------------------------------------------------------------------###

# d3 <- merge(d1, d2, by=c("snpName", "chromosome", "position"), all=TRUE)

# nrow(d3)==nrow(d1)
# TRUE

###--------------------------------------------------------------------------###
### The rs number is also avaliable on the Illumina reference 
### https://support.illumina.com/downloads/infinium-multi-ethnic-global-8-v1-support-files.html
### Infinium Multi-Ethnic Global-8 v1.0 Loci Name to rsID Conversion File (Build 37)
###--------------------------------------------------------------------------###

ref_rs <- read.csv("../GEM_Reference_Data/Multi-EthnicGlobal_D1_b150_rsids.txt", sep="\t")

# rename the column for data merging
names(ref_rs) <- c("snpName", "rsID")

# Merge the data via snpName
dd <- merge(d1, ref_rs, by="snpName", all = TRUE)

# nrow(dd)==nrow(d1)
# TRUE
###--------------------------------------------------------------------------###
### load locusReport and geneAnnotation file
###
### Note: file obtained from https://support.illumina.com/downloads/infinium-multi-ethnic-global-8-v1-support-files.html
###   1) locusReport
###   2) geneAnnotation
###--------------------------------------------------------------------------###

chip_dir <- "/zfs1/kfan/Ruyu/harmonization_Sep5/GEM_Reference_Data"

# GEM chip info: 
#   1) Multi-EthnicGlobal_D1_LocusReport.txt
#   2) Multi-EthnicGlobal_D1.annotated.txt

# read in locus report file
locusReport <- read.table(paste0(chip_dir, "/Multi-EthnicGlobal_D1_LocusReport.txt"), fill = TRUE, header = TRUE)
# read in gene annotation file
geneAnnot <- read.table(paste0(chip_dir, "/Multi-EthnicGlobal_D1.annotated.txt"), fill = TRUE, header = TRUE)

# split Alleles into a1 and a2
# note: checked with dbSNP, the [a1/a2] format in the gene annotation file does not correponding to fixed [alt/ref]
alleles <- separate(geneAnnot, Alleles, into = c("a1", "a2"), sep = "/")

# If you have brackets around the values like [A/G], you might need to remove them
geneAnnot$a1 <- gsub("\\[|\\]", "", alleles$a1)
geneAnnot$a2 <- gsub("\\[|\\]", "", alleles$a2)

# length(unique(locusReport$Name))
# 1748250
# length(intersect(locusReport$Name, geneAnnot$Name))
# 1748250
locusAnnot <- merge(locusReport[,c("Name", "Chr", "Position")], 
                    geneAnnot[,c("Name", "a1", "a2")], by.x = "Name")
# length(intersect(locusAnnot$Name, dd$snpName))
# 1748250

dd <- merge(dd, locusAnnot[,c("Name", "a1", "a2")], by.x = "snpName", by.y = "Name", all.x = TRUE)

# sort the SNPs by the chromosome and position
# "snpID should be a unique integer vector, sorted with respect to chromosome and position"
ind <- order(dd$chromosome, dd$position)
d <- dd[ind,]

# The create a snpID and it has to be in ordered 
d$snpID <- 1:nrow(d)

# Reorder the columns
d <- subset(d, select=c("snpID", "Index", "snpName", "chromosome", "position", "rsID", "a1", "a2"))

# Create a SnpAnnotationDataFrame
snpAnnot <- SnpAnnotationDataFrame(d)

# names of columns
# varLabels(snpAnnot)
# data
# head(pData(snpAnnot))

# Add metadata to describe the columns
meta <- varMetadata(snpAnnot)
meta[c("snpID", "Index", "snpName", "chromosome", "position", "rsID", "a1", "a2"),
     "labelDescription"] <- c("unique integer ID for SNPs",
                              "BeadSet SNP Index from Illumina",
                              "BeadSet SNP ID from Illumina",
                              paste("integer code for chromosome: 1:22=autosomes,",
                                    "23=X, 24=pseudoautosomal, 25=Y, 26=Mitochondrial, 27=Unknown"),
                              "base pair position on chromosome (build 37)",
                              "RS identifier",
                              "a1 allele from reference", 
                              "a2 allele from reference")
varMetadata(snpAnnot) <- meta

# str(snpAnnot)

###--------------------------------------------------------------------------###
### remove all other dataframes other than annotation files 
###   and then save the annotation data as an R image 
###--------------------------------------------------------------------------###

rm(list=setdiff(ls(),c("snpAnnot","scanAnnot")))
# ls()

save.image("GEM_GWAS_Annotation.RData")

