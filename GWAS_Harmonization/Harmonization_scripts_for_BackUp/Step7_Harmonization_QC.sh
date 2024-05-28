#!/bin/bash

###--------------------------------------------------------------------------###
###
### This is a README file for Step7: Final QC for Harmonized File
###
###	Several steps were mentioned in this step:
### (1) Update IDs, Sex, and Case/Control status
### (2) Obtain information on individuals with missing SNP data
### (3) SNP level and individual level QC
### (4) update APOE genotype with Lab APOE
### (5) sample check on race/apoe discrepancies
### 
###																 Ruyu Shi
###															   05/26/2024
###--------------------------------------------------------------------------###

mkdir -p QC_Check

###-------------------------------------------------------------------------------------------------###
### (1) Update IDs, Sex, and Case/Control status
###-------------------------------------------------------------------------------------------------###

Rscript Step7-1_ID_Sex_Race_Update.R

### update IDs

# note: when converting files between VCF and PLINK, the underscores are causing problems, need to update IDs to the original scanNames

##plink --bfile AllCohorts_AllChr_Annot --update-ids All_Cohorts_2nd_Imputation_ID_Update_forPlink.txt --make-bed --out AllCohorts_AllChr_IDcorrected_Annot

### update Sex in the .fam file

##plink --bfile AllCohorts_AllChr_IDcorrected_Annot --update-sex All_Cohorts_2nd_Imputation_Sex_Update_forPlink.txt --make-bed --out AllCohorts_AllChr_IDcorrected_SexAdded_Annot

### update Case/Control status

Rscript Step7-2_Case-Control_Update.R

###-------------------------------------------------------------------------------------------------###
### (2) Obtain information on individuals with missing SNP data
###-------------------------------------------------------------------------------------------------###

##plink --bfile AllCohorts_AllChr_IDcorrected_SexAdded_Annot --autosome --missing --allow-extra-chr --allow-no-sex --out ./QC_Check/AllCohorts_AllChr_Annot_Missing

###-------------------------------------------------------------------------------------------------###
### (3) SNP level and individual level QC
###-------------------------------------------------------------------------------------------------###

### SNP level QC
##plink --bfile AllCohorts_AllChr_IDcorrected_SexAdded_Annot --hwe 1e-6 include-nonctrl --geno 0.05 --allow-extra-chr --allow-no-sex --make-bed --out ./QC_Check/AllCohorts_AllChr_Annot_QC1
##plink --bfile ./QC_Check/AllCohorts_AllChr_Annot_QC1 --test-missing --out ./QC_Check/AllCohorts_AllChr_Case-Control

### Individual level QC
##plink --bfile ./QC_Check/AllCohorts_AllChr_Annot_QC1 --mind 0.05 --allow-extra-chr --allow-no-sex --make-bed --out ./QC_Check/AllCohorts_AllChr_Annot_QC2

### IBD Check

##plink --bfile ./QC_Check/AllCohorts_AllChr_Annot_QC2 --autosome --hwe 0.05 --indep-pairphase 50 5 0.2 --mind 0.05 --geno 0.05 --maf 0.15 --allow-extra-chr --biallelic-only --out ./QC_Check/STEP1_IBD
##plink --bfile ./QC_Check/AllCohorts_AllChr_Annot_QC2 --autosome --allow-extra-chr --extract ./QC_Check/STEP1_IBD.prune.in --make-bed --out ./QC_Check/STEP2_IBD
##plink --bfile ./QC_Check/STEP2_IBD --genome --autosome --min 0.40 --out ./QC_Check/AllCohorts_AllChr_IBD

# IBD removal 

Rscript Step7-3_IBD_Removal.R

plink --bfile ./QC_Check/AllCohorts_AllChr_Annot_QC2 --remove ./QC_Check/IBD_Remove_list.txt --make-bed --out ./QC_Check/AllCohorts_AllChr_Annot_QC3

# IBD Check again

##plink --bfile ./QC_Check/AllCohorts_AllChr_Annot_QC3 --autosome --hwe 0.05 --indep-pairphase 50 5 0.2 --mind 0.05 --geno 0.05 --maf 0.15 --allow-extra-chr --biallelic-only --out ./QC_Check/STEP1_IBD
##plink --bfile ./QC_Check/AllCohorts_AllChr_Annot_QC3 --autosome --allow-extra-chr --extract ./QC_Check/STEP1_IBD.prune.in --make-bed --out ./QC_Check/STEP2_IBD
##plink --bfile ./QC_Check/STEP2_IBD --genome --autosome --min 0.40 --out ./QC_Check/AllCohorts_AllChr_IBD2

### PCA

plink --bfile ./QC_Check/AllCohorts_AllChr_Annot_QC3 --autosome --hwe 1E-05 --indep-pairphase 20000 200 0.5 --mind 0.05 --geno 0.05 --maf 0.05 --allow-extra-chr --biallelic-only --out ./QC_Check/STEP1_PCA
plink --bfile ./QC_Check/AllCohorts_AllChr_Annot_QC3 --autosome --allow-extra-chr --biallelic-only --extract ./QC_Check/STEP1_PCA.prune.in --make-bed --out ./QC_Check/STEP2_PCA
plink --bfile ./QC_Check/STEP2_PCA --autosome --cluster --genome gz --pca header --neighbour 1 10 --out ./QC_Check/AllCohorts_AllChr_PCA

rm ./QC_Check/STEP*

Rscript Step7-4_PCA_Visualization.R

###-------------------------------------------------------------------------------------------------###
### (4) update APOE genotype with Lab APOE
###-------------------------------------------------------------------------------------------------###

### update SNP IDs using the format of chr:pos:ref:alt[b37] to get rid of the duplicated SNPs error (different alleles on the same rsID)

# rename SNPs 

module load plink/2.3_alpha
plink2 --bfile AllCohorts_AllChr_Annot_QC3 --set-all-var-ids @:#:\$r:\$a[b37] --make-bed --out AllCohorts_AllChr_Annot_QC4

module load plink/1.90b6.7
plink --bfile AllCohorts_AllChr_Annot_QC4 --write-snplist --out AllCohorts_AllChr_Annot_QC4_all_snps

### get APOE E2 and E4 allele genotype from the imputed file

# no APOE E4, found 19:45426792:G:A[b37] (rs141622900) in high LD (R2 > 0.6)

plink --bfile AllCohorts_AllChr_Annot_QC4 --snps 19:45411941:T:C[b37],19:45426792:G:A[b37] --recodeA --out AllCohorts_AllChr_Annot_APOE_SNPs
plink --bfile AllCohorts_AllChr_Annot_QC4 --recode --out AllCohorts_AllChr_Annot_QC4

### use PLINK --merge command to update APOE genotype

Rscript Step7-5_APOE_Update.R
# --merge-mode 5 to overwrite

plink --bfile AllCohorts_AllChr_Annot_QC4 --merge AllCohorts_AllChr_Annot_APOE_SNPs.ped AllCohorts_AllChr_Annot_APOE_SNPs.map --merge-mode 5 --merge-equal-pos --make-bed --out AllCohorts_AllChr_Annot_QCed

# double-check the update works
plink --bfile AllCohorts_AllChr_Annot_QCed --snps 19:45411941:T:C[b37],19:45426792:G:A[b37] --recodeA --out AllCohorts_AllChr_Annot_APOE_SNPs_updated



