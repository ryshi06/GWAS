#!/bin/bash

###--------------------------------------------------------------------------###
###
### This is a Linux script to load Illumina FinalReport.txt in to GWASTools 
###		It has several steps
###		(1) Chop one huge file into 2829 smaller text files, 
###			while each file has the genotype for only one sample.
###			Also, create AlleleA & AlleleB columns in the smaller raw data
###			for GWASTools to load currect genotype.
### 	(2) prepare the Snp and Scan Annotation dataframe which are necessary for 
###			GWASTools to work properly.
###		(3) Load the genotype, x, y, BAF, LRR, quality scores in the GWASTools, 
###			and save all the information as an R image
###
###																Frank Fan
###																6/5/2019
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
### (1) setup environment
###--------------------------------------------------------------------------###

# Assuming you are on the Geneva server
#  then move to your own directory in the GEM_GWAS project folder

cd /home/Data/GEM_GWAS/Frank

# Then you have to copy all the scripts in the /Public/Frank_Scripts into your own folder

cp ../Public/Frank_Scripts/* ./

###--------------------------------------------------------------------------###
### (2) Prepare the raw data
###--------------------------------------------------------------------------###

# Run a Linux script to prepare the raw data into proper format for later use.
# This script will take a long time to run, but it only need to be run once. 
# Therefore, this execution line is commented out (#). 
# See details comments inside the script.

./2_Raw_Data_Prep.sh

###--------------------------------------------------------------------------###
### (3) Prepare the phenotype data
###--------------------------------------------------------------------------###

# Run an R script to prepare necessary phenotype data for the QC process
# The source excel files are obtained from Dr. Snitz
# This R script gether necessary felds like sex and case/control for QC 
#	as well as some basic phenotype data for association analysis such as age
# The output is saved as an R object and will be used in creating scan annotaiton file

R CMD BATCH --vanilla 3_Phenotype_Data_from_Dr_Snitz_Arrangement.R


###--------------------------------------------------------------------------###
### (4) Running an R script to generate snp and scan Annotation dataframe
###--------------------------------------------------------------------------###

# Run an R script to create both Annotation files 
# The output is saved in the same directory named "GEM_GWAS_Annotation.RData"
# It has one SnpAnnotationDataFrame and one ScanAnnotationDataFrame GWASTools objects

R CMD BATCH --vanilla 4_Create_Annotation_DataFrame.R

###--------------------------------------------------------------------------###
### (5) Apply createDataFiles function 
###--------------------------------------------------------------------------###

# Run an R script to load the geneotype data in GWASTools

R CMD BATCH --vanilla 5_CreateDataFile.R

###--------------------------------------------------------------------------###
### (6) Run the actual QC 
###--------------------------------------------------------------------------###

# Run an R script to do QC with GWASTools

R CMD BATCH --vanilla 6_GEM_GWAS_QC.R


###--------------------------------------------------------------------------###
### (7) QC result check 
###--------------------------------------------------------------------------###

# Run an R script to check kinship 

R CMD BATCH --vanilla 7_GEM_GWAS_Kinship_Check.R





