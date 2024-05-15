#!/bin/bash

module load plink/1.90b6.7

# change the directory 

cd /ix/kfan/Ruyu/multi-trait_GWAS

### ------------------------------------------------------------------------ ###
### Set 1: Abeta group (N = 199)
### ------------------------------------------------------------------------ ###

### process genotype file 

# subset based on available phenotype
plink --bfile ./Data/Genotype/ABCDS_Cases_NonHis_WH --keep ./Data/Phenotype/ABCDS_ab_FID_IID.txt --allow-no-sex --biallelic-only --make-bed --out ./Data/Genotype/ABCDS_multitrait_Set1_AB

# calculate IBD of the original 
plink --bfile ./Data/Genotype/ABCDS_Cases_NonHis_WH --hwe 0.05 --indep-pairphase 50 5 0.2 --mind 0.05 --geno 0.05 --maf 0.15 --allow-extra-chr --biallelic-only --out ./Data/Genotype/STEP1_IBD
plink --bfile ./Data/Genotype/ABCDS_Cases_NonHis_WH --allow-extra-chr --extract ./Data/Genotype/STEP1_IBD.prune.in --make-bed --out ./Data/Genotype/STEP2_IBD
plink --bfile ./Data/Genotype/STEP2_IBD --genome --min 0.40 --out ./Data/Genotype/ABCDS_Cases_NonHis_WH_IBD

# Remove intermidiate files
rm ./Data/Genotype/STEP2_IBD.* ./Data/Genotype/STEP1_IBD.*

# convert BIMBAM
plink --bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB --snps-only --recode bimbam --out ./Data/Genotype/ABCDS_multitrait_Set1_AB

# LD pruning 
plink --bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB --indep-pairwise 20000 200 0.5 --out ./Data/Genotype/ABCDS_multitrait_Set1_AB_20000_200_05
wc -l ./Data/Genotype/ABCDS_multitrait_Set1_AB_20000_200_05.prune.in

plink --bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB --extract ./Data/Genotype/ABCDS_multitrait_Set1_AB_20000_200_05.prune.in --make-bed --out ./Data/Genotype/ABCDS_multitrait_Set1_AB_20000_200_05

# calculate PCA using plink
plink --bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated --hwe 1E-05 --indep-pairphase 20000 200 0.5 --mind 0.05 --geno 0.05 --maf 0.05 --allow-extra-chr --biallelic-only --out ./Data/Genotype/STEP1_PCA
plink --bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated --allow-extra-chr --extract ./Data/Genotype/STEP1_PCA.prune.in --make-bed --out ./Data/Genotype/STEP2_PCA
plink --bfile ./Data/Genotype/STEP2_PCA --cluster --genome gz --pca 199 header --neighbour 1 20 --out ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated_PCA

# sort the eigenval 
sort -g ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated_PCA.eigenval -o ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated_PCA_eigenval.txt
# remove the ID columns and header row of the eigenvec file
tail -n +2 ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated_PCA.eigenvec | cut -d ' ' -f 3- > ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated_PCA_eigenvec.txt

# after checking IBD, remove samples
plink --bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated --remove ./Data/Phenotype/set1AB_remove-list.txt --make-bed --out ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated_removed

# manually assign 1s to the phenotype in the fam file
awk '{$6 = ($6 == -9) ? 1 : $6}1' ./Data/Genotype/ABCDS_multitrait_Set1_AB.fam > temp && mv temp ./Data/Genotype/ABCDS_multitrait_Set1_AB.fam
awk '{$6 = ($6 == -9) ? 1 : $6}1' ./Data/Genotype/ABCDS_multitrait_Set1_AB_20000_200_05.fam > temp && mv temp ./Data/Genotype/ABCDS_multitrait_Set1_AB_20000_200_05.fam

### Estimate Relatedness Matrix from Genotypes

/ix/kfan/My_Bioinf_Apps/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated -gk 2 -n 1 -o ABCDS_multitrait_Set1_AB_updated
/ix/kfan/My_Bioinf_Apps/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated -k ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated.sXX.txt -eigen -o ABCDS_multitrait_Set1_AB_updated

/ix/kfan/My_Bioinf_Apps/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_20000_200_05 -gk 2 -n -o ABCDS_multitrait_Set1_AB_20000_200_05

/ix/kfan/My_Bioinf_Apps/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated_removed -gk 2 -n -o ABCDS_multitrait_Set1_AB_updated_removed
/ix/kfan/My_Bioinf_Apps/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated_removed -k ./output/ABCDS_multitrait_Set1_AB_updated_removed.sXX.txt -eigen -o ABCDS_multitrait_Set1_AB_updated_removed

mv ./output/* ./Data/Genotype/

### Predict Phenotype Values
# no need for set 1

### Run Analysis 

# remove header of the phenotype file
tail -n +2 ./Data/Phenotype/ABCDS_ab_pheno_residuals.txt > ./Data/Phenotype/ABCDS_ab_pheno_residuals_no_header.txt

# run analysis 
/ix/kfan/My_Bioinf_Apps/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB -maf 0.02 -k ./Data/Genotype/ABCDS_multitrait_Set1_AB.sXX.txt -lmm 4 -p ./Data/Phenotype/ABCDS_ab_pheno_residuals_no_header.txt -n 3 4 5 6 -o ABCDS_multitrait_Set1_AB_BFILE_output.txt
/ix/kfan/My_Bioinf_Apps/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated -maf 0.02 -k ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated.sXX.txt -lmm 4 -n 1 2 3 4 -o ABCDS_multitrait_Set1_AB_BFILE_sXX

/ix/kfan/My_Bioinf_Apps/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated -maf 0.02 -d ./output/ABCDS_multitrait_Set1_AB_updated.eigenD.txt -u ./output/ABCDS_multitrait_Set1_AB_updated.eigenU.txt -lmm 4 -n 1 2 3 4 -o ABCDS_multitrait_Set1_AB_BFILE_moduleLoadGEMMA
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated -maf 0.02 -d ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated_PCA_eigenval.txt -u ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated_PCA_eigenvec.txt -lmm 4 -n 1 2 3 4 -o ABCDS_multitrait_Set1_AB_BFILE_PCA_output.txt

/ix/kfan/My_Bioinf_Apps/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated -maf 0.02 -k ./Data/Genotype/ABCDS_multitrait_Set1_AB_updated.nXX.txt -lmm 4 -n 1 2 3 4 -o ABCDS_multitrait_Set1_AB_BFILE_output.txt

/ix/kfan/My_Bioinf_Apps/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_20000_200_05 -maf 0.02 -k ./Data/Genotype/ABCDS_multitrait_Set1_AB_20000_200_05.sXX.txt -lmm 4 -p ./Data/Phenotype/ABCDS_ab_pheno_residuals_no_header.txt -n 3 4 5 6 -o ABCDS_multitrait_Set1_AB_20000_200_05_output.txt

### ------------------------------------------------------------------------ ###
### Set 3: Tau group (N = 320)
### ------------------------------------------------------------------------ ###

### subset 50 individuals testing for sample size
mkdir -p ./Data/Genotype/sampleSize_test
mkdir -p ./Data/Phenotype/sampleSize_test

head -n 50 ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_prdt.fam | cut -f 1,2 > ./Data/Genotype/sampleSize_test/ABCDS_multitrait_Set3_Tau_residuals_prdt_50_Samples.txt
plink --bfile ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_prdt --keep ./Data/Genotype/sampleSize_test/ABCDS_multitrait_Set3_Tau_residuals_prdt_50_Samples.txt --make-bed --out ./Data/Genotype/sampleSize_test/ABCDS_multitrait_Set3_Tau_residuals_prdt_50_Samples

head -n 50 ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_prdt.fam > ./Data/Genotype/sampleSize_test/ABCDS_multitrait_Set3_Tau_residuals_prdt_50_Samples.fam

### calculate the kinship matrix of 50 samples
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/sampleSize_test/ABCDS_multitrait_Set3_Tau_residuals_prdt_50_Samples -gk 2 -n 1 -o ABCDS_multitrait_Set3_Tau_50_Samples
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/sampleSize_test/ABCDS_multitrait_Set3_Tau_residuals_prdt_50_Samples -k ./output/ABCDS_multitrait_Set3_Tau_50_Samples.sXX.txt -eigen -o ABCDS_multitrait_Set3_Tau_50_Samples

### run with predicted residuals on 50 samples
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/sampleSize_test/ABCDS_multitrait_Set3_Tau_residuals_prdt_50_Samples -maf 0.01 -d ./output/ABCDS_multitrait_Set3_Tau_50_Samples.eigenD.txt -u ./output/ABCDS_multitrait_Set3_Tau_50_Samples.eigenU.txt -lmm 4 -n 1 2 3 4 -o ABCDS_multitrait_Set3_Tau_residuals_prdt_50_Samples

### try phenotype imputation with AutoComplete
module load gcc/8.2.0 python/anaconda3.10-2022.10

wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py

git clone https://github.com/sriramlab/AutoComplete
cd AutoComplete
pip install -r requirements.txt

python ./tools/AutoComplete/fit.py ./Data/Phenotype/ABCDS_multitrait_Set3_Tau_residuals.csv --id_name IID --batch_size 320 --epochs 50 --lr 0.1 --seed 215 --device cpu:0 --quality --save_imputed --output ./Data/Phenotype/sampleSize_test/ABCDS_multitrait_Set3_Tau_residuals_impu.csv

## fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set3_Tau.fam")
## table(fam$V5)
## head(fam)
## impu <- fread("./Data/Phenotype/sampleSize_test/ABCDS_multitrait_Set3_Tau_residuals_impu_seed215.csv")
## impu.fam <- merge(fam, impu, by.x = "V2", by.y = "IID")
## write.table(impu.fam[, c("FID", "V2", "V3.x", "V4.x", "V5", "V1.y", "V2.y", "V3.y", "V4.y")], "./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_impu.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

cp ./Data/Genotype/ABCDS_multitrait_Set3_Tau.bim ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_impu.bim
cp ./Data/Genotype/ABCDS_multitrait_Set3_Tau.bed ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_impu.bed

### Run Analysis 
### run with imputed residuals
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_impu -maf 0.02 -d ./output/ABCDS_multitrait_Set3_Tau.eigenD.txt -u ./output/ABCDS_multitrait_Set3_Tau.eigenU.txt -lmm 4 -n 1 2 3 4 -o ABCDS_multitrait_Set3_Tau_residuals_impu

### 


### ------------------------------------------------------------------------ ###
### Set 5: Tau (plus PET) N = 326
### ------------------------------------------------------------------------ ###

awk -F',' 'BEGIN {OFS = ","} NF > 2 && ($3 != "" || $4 != "" || $5 != "" || $6 != "" || $7 != "") {print}' ./Data/Phenotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals.csv > ./Data/Phenotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals_cleaned.csv

python ./tools/AutoComplete/fit.py ./Data/Phenotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals.csv --id_name IID --batch_size 50 --epochs 50 --lr 0.1 --seed 215 --device cpu:0 --quality --save_imputed --output ./Data/Phenotype/sampleSize_test/ABCDS_multitrait_Set5_Tau_plus_PET_residuals_impu.csv

###
### ------------------------------------------------------------------------ ###
### Set 6: Tau group (N = 71)
### ------------------------------------------------------------------------ ###

### process genotype file 

# subset based on available phenotype
plink --bfile ./Data/Genotype/ABCDS_Cases_NonHis_WH --keep ./Data/Phenotype/ABCDS_multitrait_set6_Tau_plus_PET_FID_IID.txt --allow-no-sex --biallelic-only --freq --make-bed --out ./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET

# change the phenotype column in fam to 1 for GEMMA to calculate kinship matrix 
awk '{$6 = ($6 == -9) ? 1 : $6}1' ./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET.fam > temp && mv temp ./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET.fam

# calculate relatedness matrix from genotype file 

/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET -gk 2 -n 1 -o ABCDS_multitrait_set6_Tau_plus_PET
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET -k ./output/ABCDS_multitrait_set6_Tau_plus_PET.sXX.txt -eigen -o ABCDS_multitrait_set6_Tau_plus_PET

# change the 0 value to 0.00001 manually in ./output/ABCDS_multitrait_set6_Tau_plus_PET.eigenD.txt
awk 'NR==1 && $1==0 { $1=0.0000001 } { print }' ./output/ABCDS_multitrait_set6_Tau_plus_PET.eigenD.txt > temp && mv temp ./output/ABCDS_multitrait_set6_Tau_plus_PET.eigenD.txt

### process .fam phenotype in R 

## fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set3_Tau.fam")
## table(fam$V5)
## head(fam)
## residuals <- fread("./Data/Phenotype/ABCDS_multitrait_Set3_Tau_residuals.txt")
## residuals.fam <- merge(fam, residuals, by = c("V1","V2"))
## write.table(residuals.fam[,which(colnames(residuals.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

## phenos <- fread("./Data/Phenotype/ABCDS_multitrait_Set3_Tau_phenos.txt")
## phenos.fam <- merge(fam, phenos, by = c("V1","V2"))
## write.table(phenos.fam[,which(colnames(phenos.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set3_Tau_phenos.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

cp ./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET.bim ./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET_residuals.bim
cp ./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET.bed ./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET_residuals.bed

#############################################################################################
### example files 
#############################################################################################

cd /ihome/kfan/rus39/gemma-example-
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -g mouse_hs1940.geno.txt -p mouse_hs1940.pheno.txt -n 1 6 -k ./output/mouse_hs1940.cXX.txt -predict -o mouse_hs1940_prdt

plink --bfile mouse_hs1940 --recode vcf --out mouse_hs1940
tar -zxvf qctool.tgz
chmod +x ./qctool/qctool

../qctool -g mouse_hs1940.vcf -ofiletype bimbam_dosage -og mouse_hs1940_qctool.txt

/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -g mouse_hs1940_qctool.txt -p mouse_hs1940.pheno.txt -n 1 6 -k ./output/mouse_hs1940.cXX.txt -predict -o mouse_hs1940_prdt

/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile mouse_hs1940 -p mouse_hs1940.pheno.txt -n 1 6 -k ./output/mouse_hs1940.cXX.txt -predict -o mouse_hs1940_predict

