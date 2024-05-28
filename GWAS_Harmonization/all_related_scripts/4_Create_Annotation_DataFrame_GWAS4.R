
###--------------------------------------------------------------------------###
### environment settings
###--------------------------------------------------------------------------###

rm(list=ls())
options(stringsAsFactors = FALSE)

# set the working directory to ***your own*** directory 
setwd("/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS4/")

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
SN <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS4/GWAS_4_scanName.txt", sep = "\t", 
                             header = FALSE)

# create a text string with all the individual file names
SF <- paste("GWAS_4_Sample_",SN$V1,".txt", sep="")

# make a main dataframe to merge other information later
d <- data.frame(scanName=SN$V1, file=SF)
# length(unique(d$scanName))

###--------------------------------------------------------------------------###
### Other clinical information such as age, race, and gender
###   as well as the technical information like the batch and the plate
###--------------------------------------------------------------------------###

#############################
### technical informaiton ###
#############################
# not available

#############################
### clinical informaiton ###
#############################
# GWAS_Ref <- read.table("../GWAS_Reference_Data/GWAS_ALL_2023_07_24.csv", sep = ",", header = TRUE)
GWAS_4_Ref <- read.csv("GWAS_4_Clinical_Info_updated_2023_12_13.csv", header = TRUE)
colnames(GWAS_4_Ref) <- c("AlzID","scanName", "Study","Sex","Race","Disease_Status","Case_Control_Baseline","Case_Control_Current" )
# length(unique(GWAS_4_Ref$scanName))

# unmatch <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS4/GWAS_4_unmatch.txt", header = FALSE)
# GWAS_4_unmatch <- GWAS_4_Ref[GWAS_4_Ref$scanName %in% unmatch$V1, ]
# GWAS_4_unmatch0$scanName <- paste0(GWAS_4_unmatch0$scanName, " ", GWAS_4_unmatch0$V2)

# GWAS_4_Ref <- rbind(GWAS_4_Ref[!(GWAS_4_Ref$scanName %in% unmatch$V1),], GWAS_4_unmatch0[,-ncol(GWAS_4_unmatch0)])
# length(unique(GWAS_4_Ref$scanName))

GWAS_4_Manifest <- read.csv("GWAS4_Plate_Well_Manifest_from_coreLab_11172023.csv", header = TRUE)
# colnames(GWAS_4_Manifest)
table(GWAS_4_Manifest$updated_Plate)

# check sample match
d$scanName_Dup <- sapply(strsplit(d$scanName, "_"), function(x) x[1])

d_overLap <- d[d$scanName %in% GWAS_4_Manifest$updated_ID,]
# d_setdiff <- d[!(d$scanName %in% GWAS_4_Manifest$updated_ID),]

# length(intersect(GWAS_4_Manifest$updated_ID,d$scanName))

###--------------------------------------------------------------------------###
### Create ScanAnnotation object
###--------------------------------------------------------------------------###
# check ID match

dd <- merge(d, GWAS_4_Ref, by = "scanName")
# duplicated_rows <- dd[duplicated(dd$scanName) | duplicated(dd$scanName, fromLast = TRUE), ]
# dd$scanID <- 1:nrow(dd)

ddd <- merge(dd, GWAS_4_Manifest[,c("updated_ID", "updated_Plate", "Sample_Well")], 
             by.x = "scanName", by.y = "updated_ID", all.x = TRUE)

# duplicated_rows <- ddd[duplicated(ddd$scanName) | duplicated(ddd$scanName, fromLast = TRUE), ]

###-----------------------------------------------------------------------------
# after pre-processing and integrating plate/well information, extract sample IDs by plate and

d <- subset(ddd, select=c("scanName", "file", "Sex", "Race", "Case_Control_Baseline", "Case_Control_Current", "Disease_Status",
                                   "updated_Plate", "Sample_Well"))
colnames(d) <- c("scanName", "file", "sex", "race", "Case_Control_Baseline", "Case_Control_Current", "Disease_Status",
                 "plate", "well_position")


# create scanAnnot
# head(d)
d$sex <- ifelse(d$sex %in% c("F", "M"), d$sex, "F")
d$scanID <- 1:nrow(d)
dd <- subset(d, select=c("scanID", "scanName", "file", "sex", "race", "Case_Control_Baseline", "Case_Control_Current", "Disease_Status",
                         "plate", "well_position"))


# table(d$plate)
scanAnnot <- ScanAnnotationDataFrame(d)

meta <- varMetadata(scanAnnot)
meta[c("scanID","scanName","file","sex","race",
       "Case_Control_Baseline", "Case_Control_Current", "Disease_Status",
       "plate", "well_position"), "labelDescription"] <-
  c("unique ID for scans","subject identifier","raw data file","Sex","Race",
    "Case or control status at baseline","Case or control status at current","Disease status",
    "Plate information", "Well Position")
varMetadata(scanAnnot) <- meta

# create scanAnnot by plate
plate_list <- split(d, d$plate)
names(plate_list) <- paste0("plate", names(plate_list))
# list2env(plate_list, envir = .GlobalEnv)

scanAnnot_list <- list()
for (plate in 1:length(plate_list)){
  
  d_plate <- as.data.frame(plate_list[[plate]])
  d_plate$scanID <- 1:nrow(d_plate)
  
  dd <- subset(d_plate, select=c("scanID", "scanName", "file", "sex", "race", "Case_Control_Baseline", "Case_Control_Current", "Disease_Status",
                                     "plate", "well_position"))
  # have NA in sex, change to "F" first, can be fixed in QC
  # otherwise, error: Error in validObject(.Object) : 
  # invalid class “ScanAnnotationDataFrame” object: sex should have values M/F
  dd$sex <- ifelse(dd$sex %in% c("F", "M"), dd$sex, "F")
  scanAnnot_list[[plate]] <- ScanAnnotationDataFrame(dd)
 
  # add metadata to describe the columns
  meta <- varMetadata(scanAnnot_list[[plate]])
  meta[c("scanID","scanName","file","sex","race",
         "Case_Control_Baseline", "Case_Control_Current", "Disease_Status",
         "plate", "well_position"), "labelDescription"] <-
    c("unique ID for scans","subject identifier","raw data file","Sex","Race",
      "Case or control status at baseline","Case or control status at current","Disease status",
      "Plate information", "Well Position")
  varMetadata(scanAnnot_list[[plate]]) <- meta
}
names(scanAnnot_list) <- names(plate_list)

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
ref <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS_Chip_Reference/GWAS_4_5_GDA-8-v1-0_D1/infinium-global-diversity-array-8-v1-0_D1_LocusReport.txt", sep="\t", header = TRUE)

# select the columns we need, discard others
d1 <- subset(ref,select=c("Index","Name", "Chr", "Position"))

# rename the columns 
names(d1) <- c("Index", "snpName", "chromosome", "position")

# make the chromosome index as an integer: (following statement is from "DataCleaning.pdf") 
#   chromosome, an integer mapping for each chromosome, with values 1-27, mapped in order
#   from 1-22, 23=X, 24=XY (the pseudoautosomal region), 25=Y, 26=M (the mitochondrial probes), 
#   and 27=U (probes with unknown positions).                                                                      caution. See the manual pages for more details.)

d1$chromosome[d1$chromosome=="X"] <- 23
d1$chromosome[d1$chromosome=="XY"] <- 24
d1$chromosome[d1$chromosome=="Y"] <- 25
d1$chromosome[d1$chromosome=="MT"] <- 26
d1$chromosome[d1$chromosome=="0"] <- 27

d1$chromosome <- as.integer(d1$chromosome)
d <- d1[order(d1$chromosome, d1$position), ]

###--------------------------------------------------------------------------###
### load locusReport and geneAnnotation file
###
### Note: file obtained from https://support.illumina.com/downloads/infinium-global-diversity-array-support-files.html
###   1) locusReport
###   2) geneAnnotation
###--------------------------------------------------------------------------###

chip_dir <- "/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS_Chip_Reference"

# GWAS4 chip info: 
#   1) infinium-global-diversity-array-8-v1-0_D1_LocusReport.txt
#   2) infinium-global-diversity-array-8-v1-0_D1.hg19.annotated.txt

# read in locus report file
locusReport <- read.table(paste0(chip_dir, "/GWAS_4_5_GDA-8-v1-0_D1/infinium-global-diversity-array-8-v1-0_D1_LocusReport.txt"), fill = TRUE, header = TRUE)
# read in gene annotation file
geneAnnot <- read.table(paste0(chip_dir, "/GWAS_4_5_GDA-8-v1-0_D1/infinium-global-diversity-array-8-v1-0_D1.hg19.annotated.txt"), fill = TRUE, header = TRUE)

# split Alleles into a1 and a2
# note: checked with dbSNP, the [a1/a2] format in the gene annotation file does not correponding to fixed [alt/ref]
alleles <- tidyr::separate(geneAnnot, Alleles, into = c("a1", "a2"), sep = "/")

# If you have brackets around the values like [A/G], you might need to remove them
geneAnnot$a1 <- gsub("\\[|\\]", "", alleles$a1)
geneAnnot$a2 <- gsub("\\[|\\]", "", alleles$a2)

# length(unique(locusReport$Name))
# 1904599
# length(intersect(locusReport$Name, geneAnnot$Name))
# 1904599
locusAnnot <- merge(locusReport[,c("Name", "Chr", "Position")], 
                    geneAnnot[,c("Name", "a1", "a2")], by.x = "Name")
# length(intersect(locusAnnot$Name, d$snpName))
# 1904599

ref <- merge(d, locusAnnot[,c("Name", "a1", "a2")], by.x = "snpName", by.y = "Name", all.x = TRUE)

d <- ref[order(ref$chromosome, ref$position), ]
# add a snpID column as required by GWASTools
d$snpID <- 1:nrow(d)

dd <- subset(d, select=c("snpID", "Index", "snpName", "chromosome", "position", "a1", "a2"))

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
### remove all other dataframes other than annotation files 
###   and then save the annotation data as an R image 
###--------------------------------------------------------------------------###

rm(list=setdiff(ls(),c("snpAnnot","scanAnnot")))
save.image("GWAS_4_Annotation.RData")

# rm(list=setdiff(ls(),c("snpAnnot","scanAnnot_list")))
# save.image("GWAS_4_Annotation_byPlate.RData")
