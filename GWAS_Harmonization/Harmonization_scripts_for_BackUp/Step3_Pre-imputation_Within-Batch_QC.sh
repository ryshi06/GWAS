#!/bin/bash

###--------------------------------------------------------------------------###
###
### This is a README file for Step3: Perform within-batch QC using R GWASTools and HRC check tool
###
###	Three R scripts were used for this step:
### (1) Generate snp and scan Annotation dataframe
### (2) Apply CreateDataFile function 
### (3) Perform within-batch QC
### 
### One extra step (4) to format VCF files before uploading to Michigan Imputation Server
###
###																 Ruyu Shi
###															   05/26/2024
###--------------------------------------------------------------------------###

module load gcc/6.3.0 vcftools/0.1.16 plink/1.90b6.7 bcftools/1.15.1

###--------------------------------------------------------------------------###
### (1) Generate snp and scan Annotation dataframe
###--------------------------------------------------------------------------###

### Run an R script to create both Annotation files 
### The output is saved in the same directory named "cohort_Annotation.RData"
### It has one SnpAnnotationDataFrame and one ScanAnnotationDataFrame GWASTools objects

# modify the R script as needed before running 

Rscript Step3-1_Create_Annotation_DataFrame_Batch.R

###--------------------------------------------------------------------------###
### (2) Apply CreateDataFile function 
###--------------------------------------------------------------------------###

# modify the R script as needed before running 

Rscript Step3-2_CreateDataFile_Batch.R

###--------------------------------------------------------------------------###
### (3) Perform within-batch QC
###--------------------------------------------------------------------------###

# modify the R script as needed before running 

# note: certain functions in the GWASTools required in this step may be updated by the tool developer and generates no output
# 		after checking the raw script posted on Github by the developer, made some chanages accordingly and need to run the modified functions

Rscript supp_anomDetectLOH_modified.R
Rscript supp_anomIdentifyLowQuality_modified.R
Rscript supp_vcfWrite_modified.R

# run the QC

Rscript Step3-3_Within-Batch_QC.R

###--------------------------------------------------------------------------###
### (4) Pre-imputation Check
###--------------------------------------------------------------------------###

### download the HRC panel pre-imputation check tool 

# V4.3 version (the most up-to-date: chromosome X included)

wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

gzip -d HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

### convert the .vcf file to PLINK binary file sets and generate a .frq (frequency file) as required by the check tool

vcftools --vcf batch_raw.vcf --out batch_raw --plink
plink --bfile batch_raw --freq --make-bed --out batch_raw

### run the pre-imputation check tool

perl HRC-1000G-check-bim-v4.3.0/HRC-1000G-check-bim.pl -b batch_raw.bim -f batch_raw.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

# the above command will output the Run-plink.sh automatically
sh Run-plink.sh

### generate the sorted vcf.gz files by chromosome as required by Michigan Imputation Server

mkdir -p ./impu_prep

for chrom in $(seq 1 23); do

    input_file="batch_raw-updated-chr${chrom}.vcf"
    output_file="${input_file}.gz"

    # sort and compress VCF
    bcftools sort "${input_file}" -Oz -o ./impu_prep/"${output_file}"
    # create index file
    bcftools index -t ./impu_prep/"${output_file}"

    echo "Sorting or compressing for ${input_file} finished."

done









