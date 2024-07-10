#!/bin/bash

###--------------------------------------------------------------------------###
###
### This is a README file for Step5: Perform pseudoGWAS in Controls from all batches
###
###	Several steps were mentioned in this step:
### (1) Prepare for running pseudoGWAS in Controls
### (2) Run PCA in Controls only and update the 
### (3) Run pseudoGWAS in Controls
### (4) Remove SNPs with p-value less than 1E-05
###																 Ruyu Shi
###															   05/26/2024
###--------------------------------------------------------------------------###

module load gcc/6.3.0 vcftools/0.1.16 plink/1.90b6.7 bcftools/1.15.1 
module load gcc/8.2.0 r/4.0.0 
module load snpeff/4.3t

###--------------------------------------------------------------------------###
### (1) Prepare the list of Cases and Controls
###--------------------------------------------------------------------------###

# this script is used to specify cases and controls of all available samples
# also, binary pseudo-case/control phenotypes is created

Rscript Step5-1_Cross-Batch_Case_Control_Prep.R

###--------------------------------------------------------------------------###
### (2) Run PCA in Controls only and update the covariate file
###--------------------------------------------------------------------------###

### subset to Controls only 
plink --bfile All_Cohorts_allChr_exSamples --keep AllCohorts_Controls_ID_forPlink.txt --make-bed --out All_Cohorts_allChr_exSamples_Controls

### calculate PCA in Controls

plink --bfile All_Cohorts_allChr_exSamples_Controls --hwe 1E-05 --indep-pairphase 20000 200 0.5 --mind 0.05 --geno 0.05 --maf 0.01 --allow-extra-chr --biallelic-only --out STEP1_PCA

plink --bfile All_Cohorts_allChr_exSamples_Controls --allow-extra-chr --biallelic-only --extract STEP1_PCA.prune.in --make-bed --out STEP2_PCA

plink --bfile STEP2_PCA --cluster --genome gz --pca header --neighbour 1 10 --out All_Cohorts_allChr_exSamples_Controls_PCA

rm STEP*

### update the covariate file using the same R script as above

Rscript Step5-1_Cross-Batch_Case_Control_Prep.R

###--------------------------------------------------------------------------###
### (3) Run pseudoGWAS in Controls
###--------------------------------------------------------------------------###

mkdir -p ../pseudoGWAS_output

### create a variable to hold all pseudoGWAS phenotypes

phenotype=("GWAS1_2" "GWAS1_3" "GWAS1_4" "GWAS1_5" "GWAS1_GEM" "GWAS2_3" "GWAS2_4" "GWAS2_5" "GWAS2_GEM" "GWAS3_5" "GWAS3_4" "GWAS3_GEM" "GWAS4_5" "GWAS4_GEM" "GWAS5_GEM")

### run pseudoGWAS (sex, and PC1-PC8 as covariates) and extract SNPs with p-value less than 1e-05

for pheno in "${phenotype[@]}"; do

	echo "Running pseudoGWAS for between ${pheno}..."
	
	plink --bfile All_Cohorts_allChr_Controls \
	--maf 0.05 \
	--pheno All_Cohorts_Control_Phenotypes_for_PLINK.txt \
	--pheno-name ${pheno} \
	--allow-no-sex \
	--covar All_Cohorts_Control_Covariates_for_PLINK.txt \
	--covar-name sex_num_recode-PC8 \
	--logistic hide-covar \
	--out ../pseudoGWAS_output/${pheno}_pseudoGWAS

	# extract variants with p-value less than 1E-05
	awk 'BEGIN{OFS="\t"} NR==1 || ($9<1E-05) {print $0, FILENAME}' ../pseudoGWAS_output/${pheno}_pseudoGWAS.assoc.logistic > ../pseudoGWAS_output/${pheno}_pseudoGWAS_assoc_logistic_pval_threshold_1E05.txt
	# output the number of SNPs with p-value less than 1E-05
	wc -l ../pseudoGWAS_output/${pheno}_pseudoGWAS_assoc_logistic_pval_threshold_1E05.txt

done

### concatenate all SNPs and save unique SNPs into a new file for further exclusion
cd ../pseudoGWAS_output

cat $(find . -type f -name "*_pseudoGWAS_assoc_logistic_pval_threshold_1E05.txt") > pseudoGWAS_pval_threshold_1E05_SNPs.txt

awk '!seen[$2]++ {print $2}' pseudoGWAS_pval_threshold_1E05_SNPs.txt > unique_pseudoGWAS_pval_threshold_1E05_SNPs.txt
wc -l unique_pseudoGWAS_pval_threshold_1E05_SNPs.txt 

###--------------------------------------------------------------------------###
### (4) Remove SNPs with p-value less than 1E-05
###--------------------------------------------------------------------------###

### remove the header line of the SNP file
tail -n +2 ./unique_pseudoGWAS_pval_threshold_1E05_SNPs.txt > ./unique_pseudoGWAS_pval_threshold_1E05_SNPs_noheader.txt 

### remove SNPs

plink --bfile ../pseudoGWAS_input/All_Cohorts_allChr_exSamples --exclude ./unique_pseudoGWAS_pval_threshold_1E05_SNPs_noheader.txt --make-bed --out ./All_Cohorts_allChr_cleaned






