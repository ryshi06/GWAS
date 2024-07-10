#!/bin/bash

###--------------------------------------------------------------------------###
###
### This is a README file for Step2: Prepare genotype and phenotype data for within-batch QC
###
###	Two R scripts were used for this step:
### (1) Prepare genotype file
### (2) Prepare phenotype file
###																 Ruyu Shi
###															   05/26/2024
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
### (1) Prepare genotype file
###--------------------------------------------------------------------------###

### The reference file of each batch is obtained from the Illumina Webpage according to the genotype chip array version
### Available online resources:
###	https://support.illumina.com/array/downloads.html
###	https://www.strand.org.uk/
### 	Each webpage has detailed instruction and explanation

### Run a bash script to get overlapped SNPs if the reference file is not a perfect match with the genotype data
###	Detailed explanation is available in the R script

Rscript Step2-1_GenotypeRef_Prep.R

###--------------------------------------------------------------------------###
### (2) Prepare phenotype file
###--------------------------------------------------------------------------###

### Required columns including: sample_name (scanName), sample_id, Sex, Race for QC purpose
### 	as well as some basic phenotype data for association analysis such as age
### The requested phenotype file can be loaded directly into R in the next step
