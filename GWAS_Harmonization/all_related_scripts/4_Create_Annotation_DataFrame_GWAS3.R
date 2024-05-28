
###--------------------------------------------------------------------------###
### environment settings
###--------------------------------------------------------------------------###

rm(list=ls())
options(stringsAsFactors = FALSE)

# set the working directory to ***your own*** directory 
setwd("/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS3/")

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

# The reference file provided by Illumina 
ref <- read.table("GWAS_3_Ref_Overlap.txt", sep="\t", header = TRUE)

# make the chromosome as a interger
ref$chromosome <- as.integer(ref$chromosome)

###--------------------------------------------------------------------------###
### load strand and ref file
###
### Note: file obtained from https://strand.org.uk/
###   1) .strand: This contains all the variants where the match to the relevant genomic sequence >90%. The strand file contains six columns;
###               SNP ID	Chromosome	Position	%match to genome	Strand	TOP Alleles
###   2) .RefAlt: Files listing the Reference and Alternate allele mappings for the Illumina arrays (contains ref allele only) 
###--------------------------------------------------------------------------###

chip_dir <- "/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS_Chip_Info"

# gwas3 chip info: 
#   1) HumanOmni2-5-8-v1-2-A-b37.strand
#   2) HumanOmni2-5-8-v1-2-A-b37.strand.RefAlt

# read in strand file
gwas3_strand <- read.table(paste0(chip_dir, "/HumanOmni2-5-8-v1-2-A-b37-strand/HumanOmni2-5-8-v1-2-A-b37.strand"))
colnames(gwas3_strand) <- c("SNP", "Chr", "Pos", "percentMatch", "Strand", "alleleTOP")
gwas3_strand$allele1 <- substr(gwas3_strand$alleleTOP, 1, 1)
gwas3_strand$allele2 <- substr(gwas3_strand$alleleTOP, 2, 2)

gwas3_strand$a1 <- ifelse(gwas3_strand$Strand == "-" & gwas3_strand$allele1 == "A", "T", 
                          ifelse(gwas3_strand$Strand == "-" & gwas3_strand$allele1 == "T", "A", 
                                 ifelse(gwas3_strand$Strand == "-" & gwas3_strand$allele1 == "C", "G", 
                                        ifelse(gwas3_strand$Strand == "-" & gwas3_strand$allele1 == "G", "C", gwas3_strand$allele1))))

gwas3_strand$a2 <- ifelse(gwas3_strand$Strand == "-" & gwas3_strand$allele2 == "A", "T", 
                          ifelse(gwas3_strand$Strand == "-" & gwas3_strand$allele2 == "T", "A", 
                                 ifelse(gwas3_strand$Strand == "-" & gwas3_strand$allele2 == "C", "G", 
                                        ifelse(gwas3_strand$Strand == "-" & gwas3_strand$allele2 == "G", "C", gwas3_strand$allele2))))

# read in Ref/Alt file
gwas3_ref <- read.table(paste0(chip_dir, "/HumanOmni2-5-8-v1-2-A-b37-strand/HumanOmni2-5-8-v1-2-A-b37.strand.RefAlt"))
colnames(gwas3_ref) <- c("SNP", "Ref")

# manually assign Alt based on Ref 
gwas3_snpRef <- merge(gwas3_strand[,c("SNP", "Chr", "Pos", "a1", "a2")], gwas3_ref, by = 'SNP')
gwas3_snpRef$Alt <- ifelse(gwas3_snpRef$Ref == gwas3_snpRef$a1, gwas3_snpRef$a2, gwas3_snpRef$a1)

# length(intersect(gwas3_snpRef$SNP, genoData@snpAnnot@data$snpName))

ref <- merge(ref, gwas3_snpRef[,c("SNP", "Ref", "Alt")], by.x = "snpName", by.y = "SNP", all.x = TRUE)

d <- ref[order(ref$chromosome, ref$position), ]
# add a snpID column as required by GWASTools
d$snpID <- 1:nrow(d)

dd <- subset(d, select=c("snpID", "Index", "snpName", "chromosome", "position", "Ref", "Alt"))

# create a SnpAnnotationDataFrame
snpAnnot <- SnpAnnotationDataFrame(dd)

# names of columns
# varLabels(snpAnnot)
# data
# head(pData(snpAnnot))

# Add metadata to describe the columns
meta <- varMetadata(snpAnnot)
meta[c("snpID", "Index", "snpName", "chromosome", "position", "Ref", "Alt"),
     "labelDescription"] <- c("unique integer ID for SNPs",
                              "BeadSet SNP Index from Illumina",
                              "BeadSet SNP ID from Illumina",
                              paste("integer code for chromosome: 1:22=autosomes,",
                                    "23=X, 24=pseudoautosomal, 25=Y, 26=Mitochondrial, 27=Unknown"),
                              "base pair position on chromosome (build 37)", 
                              "reference allele based on strand Ref/Alt file", 
                              "alternative allele based on strand Ref/Alt file")
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
SN <- read.csv("GWAS_3_scanName.txt", header=FALSE)

# create a text string with all the individual file names
SF <- paste("GWAS_3_Sample_",SN$V1,".txt", sep="")

# make a main dataframe to merge other information later
d <- data.frame(scanName=SN$V1, file=SF)

# read in sex and race information
clinical_info <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS_Reference_Data/GWAS_ALL_2023_07_24.csv",
                            sep = ",", header = TRUE)
GWAS3_info <- clinical_info[clinical_info$GWAS_3 == "Y", c("AlzID", "STudyID", "LabID", "Sex", "Race", "Case_Control", "Disease_Status", "APOE", "GWAS_3")]

# head(GWAS3_info)
# d1 <- d[tolower(d$scanName) %in% tolower(GWAS3_info$STudyID),]
# d2 <- d[tolower(d$scanName) %in% tolower(GWAS3_info$LabID),]

GWAS3_info$labid <- tolower(GWAS3_info$LabID)
dd <- merge(d, GWAS3_info, by.x = "scanName", by.y = "labid", all.x = TRUE)

d <- unique(dd)
# length(unique(d$scanName))

# add the scanID column as required by GWASTools
d$scanID <- 1:nrow(d)
dd <- subset(d, select=c("scanID", "scanName", "file", "Sex", "Race", "AlzID", "STudyID", "LabID", "Case_Control", "Disease_Status", "APOE"))
colnames(dd) <- c("scanID", "scanName", "file", "sex", "race", "AlzID", "STudyID", "LabID", "Case_Control", "Disease_Status", "APOE")

# table(dd$sex)
# have "-12" in sex, change to "F" first, can be fixed in QC
dd$sex <- ifelse(dd$sex %in% c("F", "M"), dd$sex, "F")
# create a ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(dd)

# names of columns
# varLabels(scanAnnot)
# data
# head(pData(scanAnnot))

# Add metadata to describe the columns
meta <- varMetadata(scanAnnot)
meta[c("scanID","scanName","file","sex","race",
       "AlzID", "STudyID", "LabID", "Case_Control", 
       "Disease_Status", "APOE"), "labelDescription"] <-
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
    "APOE genotype")
varMetadata(scanAnnot) <- meta

# str(scanAnnot)

###--------------------------------------------------------------------------###
### remove all other dataframes other than annotation files 
###   and then save the annotation data as an R image 
###--------------------------------------------------------------------------###

rm(list=setdiff(ls(),c("snpAnnot","scanAnnot")))
# ls()

save.image("GWAS_3_Annotation.RData")
