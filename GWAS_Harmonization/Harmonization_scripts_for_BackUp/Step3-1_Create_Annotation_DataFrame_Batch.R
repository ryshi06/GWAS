
###--------------------------------------------------------------------------###
### environment settings
###--------------------------------------------------------------------------###

rm(list=ls())
options(stringsAsFactors = FALSE)

# set the working directory to ***your own*** directory 
setwd("/zfs1/User/Batch/")

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

# The overlapped reference file processed before
ref <- read.table("Batch_Ref_Overlap.txt", sep="\t", header = TRUE)

# make the chromosome as a interger
ref$chromosome <- as.integer(ref$chromosome)

###--------------------------------------------------------------------------###
### load locusReport and geneAnnotation file
###
### Note: file obtained from https://support.illumina.com/downloads/infinium-global-diversity-array-support-files.html
###   1) locusReport
###   2) geneAnnotation
###--------------------------------------------------------------------------###

chip_dir <- "/zfs1/User/GWAS_Chip_Info"

#############################################
### reference file obtained from Illumina ###
#############################################

# Batch chip info: 
#   1) infinium-global-diversity-array-8-v1-0_D1_LocusReport.txt
#   2) infinium-global-diversity-array-8-v1-0_D1.hg19.annotated.txt

# read in locus report file
locusReport <- read.table(paste0(chip_dir, "/infinium-global-diversity-array-8-v1-0_D1_LocusReport.txt"), fill = TRUE, header = TRUE)
# read in gene annotation file
geneAnnot <- read.table(paste0(chip_dir, "/infinium-global-diversity-array-8-v1-0_D1.hg19.annotated.txt"), fill = TRUE, header = TRUE)

# split Alleles into a1 and a2
# note: checked with dbSNP, the [a1/a2] format in the gene annotation file does not correponding to fixed [alt/ref]
alleles <- separate(geneAnnot, Alleles, into = c("a1", "a2"), sep = "/")

# If you have brackets around the values like [A/G], you might need to remove them
geneAnnot$a1 <- gsub("\\[|\\]", "", alleles$a1)
geneAnnot$a2 <- gsub("\\[|\\]", "", alleles$a2)

# length(unique(locusReport$Name))
# length(intersect(locusReport$Name, geneAnnot$Name))

locusAnnot <- merge(locusReport[,c("Name", "Chr", "Position")], 
                    geneAnnot[,c("Name", "a1", "a2")], by.x = "Name")
# length(intersect(locusAnnot$Name, ref$snpName))

ref <- merge(ref, locusAnnot[,c("Name", "a1", "a2")], by.x = "snpName", by.y = "Name", all.x = TRUE)

d <- ref[order(ref$chromosome, ref$position), ]
# add a snpID column as required by GWASTools
d$snpID <- 1:nrow(d)

dd <- subset(d, select=c("snpID", "Index", "snpName", "chromosome", "position", "a1", "a2"))

###############################################################
### reference file obtained from https://www.strand.org.uk/ ###
###############################################################

# Batch chip info: 
#   1) HumanOmni2-5-8-v1-2-A-b37.strand
#   2) HumanOmni2-5-8-v1-2-A-b37.strand.RefAlt

# read in strand file
# strand <- read.table(paste0(chip_dir, "/HumanOmni2-5-8-v1-2-A-b37-strand/HumanOmni2-5-8-v1-2-A-b37.strand"))
# colnames(strand) <- c("SNP", "Chr", "Pos", "percentMatch", "Strand", "alleleTOP")
# strand$allele1 <- substr(strand$alleleTOP, 1, 1)
# strand$allele2 <- substr(strand$alleleTOP, 2, 2)

# strand$a1 <- ifelse(strand$Strand == "-" & strand$allele1 == "A", "T", 
#                           ifelse(strand$Strand == "-" & strand$allele1 == "T", "A", 
#                                  ifelse(strand$Strand == "-" & strand$allele1 == "C", "G", 
#                                         ifelse(strand$Strand == "-" & strand$allele1 == "G", "C", strand$allele1))))

# strand$a2 <- ifelse(strand$Strand == "-" & strand$allele2 == "A", "T", 
#                           ifelse(strand$Strand == "-" & strand$allele2 == "T", "A", 
#                                  ifelse(strand$Strand == "-" & strand$allele2 == "C", "G", 
#                                         ifelse(strand$Strand == "-" & strand$allele2 == "G", "C", strand$allele2))))

# read in Ref/Alt file
# ref <- read.table(paste0(chip_dir, "/HumanOmni2-5-8-v1-2-A-b37-strand/HumanOmni2-5-8-v1-2-A-b37.strand.RefAlt"))
# colnames(ref) <- c("SNP", "Ref")

# manually assign Alt based on Ref 
# snpRef <- merge(strand[,c("SNP", "Chr", "Pos", "a1", "a2")], ref, by = 'SNP')
# snpRef$Alt <- ifelse(snpRef$Ref == snpRef$a1, snpRef$a2, snpRef$a1)

# length(intersect(snpRef$SNP, genoData@snpAnnot@data$snpName))

# ref <- merge(ref, snpRef[,c("SNP", "Ref", "Alt")], by.x = "snpName", by.y = "SNP", all.x = TRUE)

# d <- ref[order(ref$chromosome, ref$position), ]

# add a snpID column as required by GWASTools
# d$snpID <- 1:nrow(d)

dd <- subset(d, select=c("snpID", "Index", "snpName", "chromosome", "position", "Ref", "Alt"))

# create a SnpAnnotationDataFrame
snpAnnot <- SnpAnnotationDataFrame(dd)

# names of columns
# varLabels(snpAnnot)
# data
# head(pData(snpAnnot))

# Add metadata to describe the columns
meta <- varMetadata(snpAnnot)
meta[c("snpID", "Index", "snpName", "chromosome", "position", "a1", "a2"),
     "labelDescription"] <- c("unique integer ID for SNPs",
                              "BeadSet SNP Index from Illumina",
                              "BeadSet SNP ID from Illumina",
                              paste("integer code for chromosome: 1:22=autosomes,",
                                    "23=X, 24=pseudoautosomal, 25=Y, 26=Mitochondrial, 27=Unknown"),
                              "base pair position on chromosome (build 37)",
                              "a1 allele from reference", 
                              "a2 allele from reference")
varMetadata(snpAnnot) <- meta 

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
### replace the actual Batch name in the following command line
SN <- read.csv("Batch_scanName.txt", header=FALSE)

# create a text string with all the individual file names
SF <- paste("Batch_Sample_",SN$V1,".txt", sep="")

# make a main dataframe to merge other information later
d <- data.frame(scanName=SN$V1, file=SF)

# read-in phenotype information
### replace the actual phenotype file in the following command line
clinical_info <- read.table("/zfs1/User/pheno.csv",
                            sep = ",", header = TRUE)
Batch_info <- clinical_info[, c("AlzID", "STudyID", "LabID", "Sex", "Race", "Case_Control", "Disease_Status", "APOE", "Batch")]

# head(Batch_info)
# d1 <- d[tolower(d$scanName) %in% tolower(Batch_info$STudyID),]
# d2 <- d[tolower(d$scanName) %in% tolower(Batch_info$LabID),]
# setdiff(tolower(d$scanName), tolower(Batch_info$LabID))

Batch_info$LABID <- toupper(Batch_info$LabID)
dd <- merge(d, Batch_info, by.x = "scanName", by.y = "LABID", all.x = TRUE)

d <- unique(dd)
# length(unique(d$scanName))

# add the scanID column as required by GWASTools
d$scanID <- 1:nrow(d)
dd <- subset(d, select=c("scanID", "scanName", "file", "Sex", "Race", "AlzID", "STudyID", "LabID", "Case_Control", "Disease_Status", "APOE"))
colnames(dd) <- c("scanID", "scanName", "file", "sex", "race", "AlzID", "STudyID", "LabID", "Case_Control", "Disease_Status", "APOE")

# table(dd$sex)
# have NA in sex, change to "F" first, can be fixed in QC
dd$sex <- ifelse(dd$sex %in% c("F", "M"), dd$sex, "F")

# read-in plate information
plate_info <- read.table("/zfs1/User/harmonization_Sep5/GWAS_Reference_Data/2023-02-06_DNA_Manifest-Kamboh-GWAS5.csv", 
                          sep = ",", header = TRUE, skip = 5)
colnames(plate_info)
Batch_plate <- plate_info[, c("Institute.Plate.Label", "Well", "Institute.Sample.Label")]
colnames(Batch_plate) <- c("plate", "well", "scanName")

# length(intersect(Batch_plate$Institute.Sample.Label, dd$scanName))
d <- merge(dd, Batch_plate, by = "scanName", all.x=TRUE)
# colnames(d)

# create a ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(d)

# names of columns
# varLabels(scanAnnot)
# data
# head(pData(scanAnnot))

# Add metadata to describe the columns
meta <- varMetadata(scanAnnot)
meta[c("scanID","scanName","file","sex","race",
       "AlzID", "STudyID", "LabID", "Case_Control", 
       "Disease_Status", "APOE", "plate", "well"), "labelDescription"] <-
  c("unique ID for scans",
    "subject identifier",
    "raw data file",
    "Sex",
    "Race",
    "Alz identifier",
    "Study identifier",
    "Lab identifier",
    "Case or control status",
    "Alzheimer status",
    "APOE genotype",
    "Institute plate info",
    "Institute well info")
varMetadata(scanAnnot) <- meta

# str(scanAnnot)

###--------------------------------------------------------------------------###
### remove all other dataframes other than annotation files 
###   and then save the annotation data as an R image 
###--------------------------------------------------------------------------###

rm(list=setdiff(ls(),c("snpAnnot","scanAnnot")))
# ls()

save.image("Batch_Annotation.RData")
