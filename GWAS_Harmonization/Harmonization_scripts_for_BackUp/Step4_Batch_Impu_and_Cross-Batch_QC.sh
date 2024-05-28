#!/bin/bash

###--------------------------------------------------------------------------###
###
### This is a README file for Step4: Upload files to Michigan Imputation Server and process imputation output
###
### This file contains the main parameters used to do the Within-Batch imputation
### 
### Several steps were mentioned in this step:
### (1) Submit to the Michigan Impuation Server
### (2) Download the imputation outputs from the Michigan Impuation Server
### (3) Concatenate each batch imputation output
### (4) Cross-batch sample filtering by exact same scanName
### (5) Merge vcf.gz of all available batches and convert to PLINK binary file sets
###
###																 Ruyu Shi
###															   05/26/2024
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
### (1) Submit to the Michigan Impuation Server
###--------------------------------------------------------------------------###

# Michigan Imputation Server: https://imputationserver.sph.umich.edu/index.html#!
# 
# 	Register to use the server
# 	Michigan Imputation Server imputation uses Minimac4 
#
# 	Once data is ready (variants data from the pre-imputation step), under the ‘Run’ tab, select "Genotype Imputation":
#
# 	Enter Job name (optional)
#   Reference panel: HRC r1.1 2016 (GRCh37/hg19)
# 	Input Files (VCF)
# 	Array Build:  GRCh37/hg19
# 	rsq Filter: off
# 	Phasing: Eagle v2.4 (phased output) (if chrX is included)
# 	Population: Others/Mixed
# 	Mode: Quality Control & Imputation
# 	AES 256 encryption checked
# 	Generate Meta-imputation file checked
# 	Select - I will not attempt to re-identify or contact research participants.
# 	Select - I will report any inadvertent data release, security breach or other data management incident of which I become aware.
# 	Click “Submit Job”
# 
# 	Wait for files to upload correctly (without reported error)
# 		This step may take a while if the input files are large
# 	Errors during imputation will be reported onscreen and sent to your login email
# 	Troubleshooting is available through the ‘help’ page (https://imputationserver.readthedocs.io/en/latest/prepare-your-data/), check ‘Data Preparation > Quality Control for HRC, 1000G and CAAPA imputation’ 
# 	When imputation is completed, a notification e-mail will be sent with password.

#   Important Note: One account can only upload 3 jobs each time, please submitting jobs accordingly

###--------------------------------------------------------------------------###
### (2) Download the imputation outputs from the Michigan Impuation Server
###--------------------------------------------------------------------------###

module load p7zip/16.02

mkdir -p ./impu_output/Batch
cd ./impu_output/Batch

### The Michigan Imputation Server provides three sets of outputs as listed below and corresponding command to download
### 	the password of the imputation results will be sent to the email used to login

# QC statistics

# Imputation Results: the password will be mentioned in the email, repace the 'password' with the actual password in the command below

for file in chr_*.zip; do     
    7z e -p'password' "$file"; 
done

rm chr_*.zip

# Logs

###--------------------------------------------------------------------------###
### (3) Concatenate each batch imputation output
###--------------------------------------------------------------------------###

module load gcc/6.3.0 vcftools/0.1.16 plink/1.90b6.7 bcftools/1.15.1 
module load gcc/8.2.0 r/4.0.0 
module load snpeff/4.3t

### set a batch variable for further use

batches=("GWAS1" "GWAS2" "GWAS3" "GWAS4" "GWAS5" "GEM")

### process each batch using a loop 

for batch in "${batches[@]}"; do
    
    echo "Processing Batch ${batch}..."

    cd ./impu_output/${batch}/
    mkdir -p ../after_imputation

    for i in {1..22} X; do  
    
    echo "Processing Chromosome ${i}..."

    vcftools --gzvcf ./chr${i}.dose.vcf.gz --recode --recode-INFO-all --out ../after_imputation/${batch}_Chr${i}_All_SNPs

    bgzip ../after_imputation_02_01_2024/${batch}_Chr${i}_All_SNPs.recode.vcf  
    tabix -p vcf ../after_imputation_02_01_2024/${batch}_Chr${i}_All_SNPs.recode.vcf.gz
    
    done

    echo "Processing Batch ${batch}..."
    
    ### concatenate all the chromosomes together

    cd ./${batch}/after_imputation
    
    echo '${batch}_Chr1_All_SNPs.recode.vcf.gz
    ${batch}_Chr2_All_SNPs.recode.vcf.gz
    ${batch}_Chr3_All_SNPs.recode.vcf.gz
    ${batch}_Chr4_All_SNPs.recode.vcf.gz
    ${batch}_Chr5_All_SNPs.recode.vcf.gz
    ${batch}_Chr6_All_SNPs.recode.vcf.gz
    ${batch}_Chr7_All_SNPs.recode.vcf.gz
    ${batch}_Chr8_All_SNPs.recode.vcf.gz
    ${batch}_Chr9_All_SNPs.recode.vcf.gz
    ${batch}_Chr10_All_SNPs.recode.vcf.gz
    ${batch}_Chr11_All_SNPs.recode.vcf.gz
    ${batch}_Chr12_All_SNPs.recode.vcf.gz
    ${batch}_Chr13_All_SNPs.recode.vcf.gz
    ${batch}_Chr14_All_SNPs.recode.vcf.gz
    ${batch}_Chr15_All_SNPs.recode.vcf.gz
    ${batch}_Chr16_All_SNPs.recode.vcf.gz
    ${batch}_Chr17_All_SNPs.recode.vcf.gz
    ${batch}_Chr18_All_SNPs.recode.vcf.gz
    ${batch}_Chr19_All_SNPs.recode.vcf.gz
    ${batch}_Chr20_All_SNPs.recode.vcf.gz
    ${batch}_Chr21_All_SNPs.recode.vcf.gz
    ${batch}_Chr22_All_SNPs.recode.vcf.gz
    ${batch}_ChrX_All_SNPs.recode.vcf.gz' > ${batch}_allChr_VCFs.txt

    # concatenate 1-22,X chromosomes of each Batch
    bcftools concat -f ${batch}_allChr_VCFs.txt -Oz -o ${batch}_allChr.vcf.gz
    tabix -p vcf ${batch}_allChr.vcf.gz
    
 done

###--------------------------------------------------------------------------###
### (4) Cross-batch sample filtering 
###--------------------------------------------------------------------------###

### the exact same sample name will lead to error when converting VCF to PLINK in the further steps
### if duplicated scanName is found, remove the older one 

### generate list of samples with sex mismatch information

Rscript Step4-1_Cross-batch_Check.R

###--------------------------------------------------------------------------###
### (5) Merge vcf.gz of all available batches and convert to PLINK binary file sets
###--------------------------------------------------------------------------###

mkdir -p ./pseudoGWAS_input
cd ./pseudoGWAS_input

### set a batch variable for further use

batches=("GWAS1" "GWAS2" "GWAS3" "GWAS4" "GWAS5" "GEM")

### merge each batch together

for batch in "${batch[@]}"; do

    bcftools view -S ../after_imputation/${batch}_plink_IDs_noDup.txt ../after_imputation/${batch}_allChr.vcf.gz > ./${batch}_allChr_noDup.vcf.gz

    # get the scanName of each batch vcf.gz
    bcftools query -l ./${batch}_allChr_scanName.txt
done

echo 'GEM_allChr_noDup.vcf.gz
    GWAS1_allChr_noDup.vcf.gz
    GWAS2_allChr_noDup.vcf.gz
    GWAS3_allChr_noDup.vcf.gz
    GWAS4_allChr_noDup.vcf.gz
    GWAS5_allChr_noDup.vcf.gz' > All_Cohorts_allChr_VCFs.txt

bcftools merge --force-samples --file-list All_Cohorts_allChr_VCFs.txt -Oz -o All_Cohorts_allChr.vcf.gz
tabix -p vcf All_Cohorts_allChr.vcf.gz

bcftools query -l All_Cohorts_allChr.vcf.gz | wc -l 

echo "Merging done for all batches!"

### convert vcf to PLINK binary file sets

plink --vcf All_Cohorts_allChr.vcf.gz --double-id --make-bed --out All_Cohorts_allChr

### remove samples with sex mismatch

plink --bfile All_Cohorts_allChr --remove All_Cohorts_ScanName_Sex_Mismatch_plink_ex.txt --make-bed --out All_Cohorts_allChr_exSamples





