
###--------------------------------------------------------------------------###
### environment settings
###--------------------------------------------------------------------------###

rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/zfs1/kfan/Ruyu/harmonization_Sep5/")

# BiocManager::install("GWASTools", force=TRUE)
library(GWASTools)
# BiocManager::install("gdsfmt", force=TRUE)
library(gdsfmt)

###--------------------------------------------------------------------------###
### This script is testing the createDataFile function in the GWASTools package
###   with smaller dataset 
### The goal is making sure everything is running perfectly, and the
###   data is stored in the proper way  
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
### (1) Annotation files
###--------------------------------------------------------------------------###

# Loading Annotation files 
load("./GEM/GEM_GWAS_Annotation.RData")

# retrieve the annotation files and save as dataframes
snp_annotation <- getAnnotation(snpAnnot)
scan_annotation <- getAnnotation(scanAnnot)

###--------------------------------------------------------------------------###
### (2) Loading Genotypes 
###--------------------------------------------------------------------------###

# assign the path to the data
path_to_data <- "./GEM/Raw_Data/"

# create a gds file to store genotypes in the same path
geno.file <- "./GEM/geno.gds"

# create an R image for the diagnostic file
diag.geno.file <- "diag.geno.RData"

# Assign the column number to the corresponding feilds and name
col.nums <- as.integer(c(1,2,10,11))
names(col.nums) <- c("snp", "sample", "a1", "a2")

# Reading genotype data and create a gds file name geno.gds

diag.geno <- createDataFile(path=path_to_data, geno.file, file.type="gds",
                            variables="genotype",
                            snp.annotation=snp_annotation,
                            scan.annotation=scan_annotation, 
                            sep.type="\t",
                            skip.num=10, col.total=11, col.nums=col.nums,
                            scan.name.in.file=1, 
                            diagnostics.filename=diag.geno.file,
                            verbose=TRUE)
###--------------------------------------------------------------------------###
### (3) checkGenotypeFile function  
###--------------------------------------------------------------------------###

# create an R image for checking genotype data
# check.geno.file <- "check.geno.RData"

# Select a subset of the gds file to check the genotype
# check.geno <- checkGenotypeFile(path = path_to_data,
#                                 filename = geno.file,
#                                 file.type = "gds",
#                                 snp.annotation = snp.annotation,
#                                 scan.annotation = scan.annotation,
#                                 sep.type = "\t",
#                                 skip.num = 10,
#                                 col.total = 11,
#                                 col.nums = col.nums,
#                                 scan.name.in.file = 1,
#                                 check.scan.index = 1:2,
#                                 n.scans.loaded = 100,
#                                 diagnostics.filename = check.geno.file,
#                                 verbose = FALSE)

###--------------------------------------------------------------------------###
### (4) Loading Intensity Variables 
###--------------------------------------------------------------------------###

# create a gds file to store x, y, and quality scores in the same path
qxy.file <- "./GEM/qxy.gds"

# create an R image for the diagnostic file
diag.qxy.file <- "diag.qxy.RData"

# Assign the column number to the corresponding feilds and name
col.nums <- as.integer(c(1,2,3,4,9))
names(col.nums) <- c("snp", "sample", "X", "Y", "quality")

# Reading x, y, and quality scores data and create a gds file name qxy.gds
# DON'T RUN 
# It will take 10+ hours to run 
diag.qxy <- createDataFile(path = path_to_data,
                           filename = qxy.file,
                           file.type = "gds",
                           variables = c("X", "Y", "quality"),
                           snp.annotation = snp_annotation,
                           scan.annotation = scan_annotation,
                           sep.type = "\t",
                           skip.num = 10,
                           col.total = 11,
                           col.nums = col.nums,
                           scan.name.in.file = 1,
                           diagnostics.filename = diag.qxy.file,
                           verbose = FALSE)

###--------------------------------------------------------------------------###
### (5) checkIntensityFile function  
###--------------------------------------------------------------------------###

# create an R image for checking quality data
# check.qxy.file <- "check.qxy.RData"
# check.qxy <- checkIntensityFile(path = path_to_data,
#                                 filename = qxy.file,
#                                 file.type = "gds",
#                                 snp.annotation = snp.annotation,
#                                 scan.annotation = scan.annotation,
#                                 sep.type = "\t",
#                                 skip.num = 10,
#                                 col.total = 11,
#                                 col.nums = col.nums,
#                                 scan.name.in.file = 1,
#                                 check.scan.index = 1:2,
#                                 n.scans.loaded = 2,
#                                 diagnostics.filename = check.qxy.file,
#                                 verbose = FALSE)

###--------------------------------------------------------------------------###
### (6) Loading BAF &n LRR Variables
###--------------------------------------------------------------------------###

# create a gds file to store BAF and LRR in the same path
bl.file <- "./GEM/bl.gds"

# create an R image for the diagnostic file
diag.bl.file <- "diag.bl.RData"

# Assign the column number to the corresponding feilds and name
col.nums <- as.integer(c(1,2,5,6))
names(col.nums) <- c("snp", "sample", "BAlleleFreq", "LogRRatio")

# Reading BAF and LRR and create a gds file name bl.gds
# DON'T RUN 
# It will take 10+ hours to run 
diag.bl <- createDataFile(path = path_to_data,
                          filename = bl.file,
                          file.type = "gds",
                          variables = c("BAlleleFreq","LogRRatio"),
                          snp.annotation = snp_annotation,
                          scan.annotation = scan_annotation,
                          sep.type = "\t",
                          skip.num = 10,
                          col.total = 11,
                          col.nums = col.nums,
                          scan.name.in.file = 1,
                          diagnostics.filename = diag.bl.file,
                          verbose = FALSE)

###--------------------------------------------------------------------------###
### (7) cleanup data files
###--------------------------------------------------------------------------###

### DON'T RUN 
### It will delete the files that cost days to read
### file.remove(geno.file, qxy.file, bl.file)
### file.remove(diag.geno.file, diag.qxy.file, diag.bl.file)
### file.remove(check.geno.file, check.qxy.file)

###--------------------------------------------------------------------------###
### (8) combining data files 
###--------------------------------------------------------------------------###

# showfile.gds(closeall=TRUE)
# get the genotype from the gds file created previously
# gds1 <- GdsGenotypeReader(geno.file)
# gds1 <- openfn.gds(geno.file, allow.error = TRUE)
# annotated with scan and SNP annotation
# genoData <- GenotypeData(gds1, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
# genoData
# close(genoData)

# get the quality score, X and X intensity from the gds file created previously
# gds2 <- GdsIntensityReader(qxy.file)
# annotated with scan and SNP annotation
# qxyData <- IntensityData(gds2, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
# qxyData
# close(qxyData)

# # get the BAF and LRR from the gds file created previously
# gds3 <- GdsIntensityReader(bl.file)
# annotated with scan and SNP annotation
# blData <- IntensityData(gds3, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
# blData
# close(blData)

