#!/bin/bash

module load plink/1.90b6.7

# change the directory 

cd /ix/kfan/Ruyu/multi-trait_GWAS

# calculate IBD of the original 
plink --bfile ./Data/Genotype/ABCDS_Cases_NonHis_WH --hwe 0.05 --indep-pairphase 50 5 0.2 --mind 0.05 --geno 0.05 --maf 0.15 --allow-extra-chr --biallelic-only --out ./Data/Genotype/STEP1_IBD
plink --bfile ./Data/Genotype/ABCDS_Cases_NonHis_WH --allow-extra-chr --extract ./Data/Genotype/STEP1_IBD.prune.in --make-bed --out ./Data/Genotype/STEP2_IBD
plink --bfile ./Data/Genotype/STEP2_IBD --genome --min 0.40 --out ./Data/Genotype/ABCDS_Cases_NonHis_WH_IBD
# Remove intermidiate files
rm ./Data/Genotype/STEP2_IBD.* ./Data/Genotype/STEP1_IBD.*

### use an older Version of GEMMA: GEMMA Ver.0.94
# downloaded from https://xiangzhou.github.io/software/gemma-0.94.tar.gz
# tar -xvzf gemma-0.94.tar.gz

### download qctools
# downloaded from https://www.chg.ox.ac.uk/~gav/qctool_v2/
# tar -zxvf qctool.tgz
# chmod +x ./qctool/qctool

### ------------------------------------------------------------------------ ###
### Set 1: Abeta group (complete) (N = 199)
### ------------------------------------------------------------------------ ###

### process genotype file 

# subset based on available phenotype
plink --bfile ./Data/Genotype/ABCDS_Cases_NonHis_WH --keep ./Data/Phenotype/ABCDS_multitrait_Set1_AB_FID_IID.txt --allow-no-sex --biallelic-only --make-bed --out ./Data/Genotype/ABCDS_multitrait_Set1_AB

# change the phenotype column in fam to 1 for GEMMA to calculate kinship matrix 
awk '{$6 = ($6 == -9) ? 1 : $6}1' ./Data/Genotype/ABCDS_multitrait_Set1_AB.fam > temp && mv temp ./Data/Genotype/ABCDS_multitrait_Set1_AB.fam

# calculate relatedness matrix from genotype file 

/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB -gk 2 -n 1 -o ABCDS_multitrait_Set1_AB
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB -k ./output/ABCDS_multitrait_Set1_AB.sXX.txt -eigen -o ABCDS_multitrait_Set1_AB

# change the 0 value to 0.00001 manually in ./output/ABCDS_multitrait_Set1_AB.eigenD.txt
awk 'NR==1 && $1==0 { $1=0.0000001 } { print }' ./output/ABCDS_multitrait_Set1_AB.eigenD.txt > temp && mv temp ./output/ABCDS_multitrait_Set1_AB.eigenD.txt

### process .fam phenotype in R 

## fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set1_AB.fam")
## table(fam$V5)
## head(fam)
## residuals <- fread("./Data/Phenotype/ABCDS_multitrait_Set1_AB_residuals.txt")
## residuals.fam <- merge(fam, residuals, by = c("V1","V2"))
## write.table(residuals.fam[,which(colnames(residuals.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set1_AB_residuals.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

## phenos <- fread("./Data/Phenotype/ABCDS_multitrait_Set1_AB_phenos.txt")
## phenos.fam <- merge(fam, phenos, by = c("V1","V2"))
## write.table(phenos.fam[,which(colnames(phenos.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set1_AB_phenos.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

cp ./Data/Genotype/ABCDS_multitrait_Set1_AB.bim ./Data/Genotype/ABCDS_multitrait_Set1_AB_residuals.bim
cp ./Data/Genotype/ABCDS_multitrait_Set1_AB.bed ./Data/Genotype/ABCDS_multitrait_Set1_AB_residuals.bed

cp ./Data/Genotype/ABCDS_multitrait_Set1_AB.bim ./Data/Genotype/ABCDS_multitrait_Set1_AB_phenos.bim
cp ./Data/Genotype/ABCDS_multitrait_Set1_AB.bed ./Data/Genotype/ABCDS_multitrait_Set1_AB_phenos.bed

### generate BIMBAM format files 
# convert plink to VCF
plink --bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_residuals --recode vcf --out ./Data/Genotype/ABCDS_multitrait_Set1_AB_residuals
# convert VCF to BIMBAM
./qctool -g ./Data/Genotype/ABCDS_multitrait_Set1_AB_residuals.vcf -ofiletype bimbam_dosage -og ./Data/Genotype/ABCDS_multitrait_Set1_AB_residuals.txt

### Run Analysis 
### run with residuals
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_residuals -maf 0.02 -d ./output/ABCDS_multitrait_Set1_AB.eigenD.txt -u ./output/ABCDS_multitrait_Set1_AB.eigenU.txt -lmm 4 -n 1 2 3 4 -o ABCDS_multitrait_Set1_AB_residuals
### run with raw phenotypes
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set1_AB_phenos -maf 0.02 -d ./output/ABCDS_multitrait_Set1_AB.eigenD.txt -u ./output/ABCDS_multitrait_Set1_AB.eigenU.txt -n 1 2 3 4 -c ./Data/Phenotype/ABCDS_multitrait_Set1_AB_covars.txt -lmm 4 -o ABCDS_multitrait_Set1_AB_phenos

### Subset SNPs with p-value less than 5e-08
awk -F'\t' '($22 < 5e-08 || $23 < 5e-08 || $24 < 5e-08)' ./output/ABCDS_multitrait_Set1_AB_BFILE_output.txt.assoc.txt > ./output/ABCDS_multitrait_Set1_AB_BFILE_output.txt.assoc.sig.txt

### ------------------------------------------------------------------------ ###
### Set 2: Tau group (complete) (N = 76)
### ------------------------------------------------------------------------ ###

### process genotype file 

# subset based on available phenotype
plink --bfile ./Data/Genotype/ABCDS_Cases_NonHis_WH --keep ./Data/Phenotype/ABCDS_multitrait_Set2_Tau_FID_IID.txt --allow-no-sex --biallelic-only --make-bed --out ./Data/Genotype/ABCDS_multitrait_Set2_Tau

# change the phenotype column in fam to 1 for GEMMA to calculate kinship matrix 
awk '{$6 = ($6 == -9) ? 1 : $6}1' ./Data/Genotype/ABCDS_multitrait_Set2_Tau.fam > temp && mv temp ./Data/Genotype/ABCDS_multitrait_Set2_Tau.fam

# calculate relatedness matrix from genotype file 

/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set2_Tau -gk 2 -n 1 -o ABCDS_multitrait_Set2_Tau
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set2_Tau -k ./output/ABCDS_multitrait_Set2_Tau.sXX.txt -eigen -o ABCDS_multitrait_Set2_Tau

# change the 0 value to 0.00001 manually in ./output/ABCDS_multitrait_Set2_Tau.eigenD.txt
awk 'NR==1 && $1==0 { $1=0.0000001 } { print }' ./output/ABCDS_multitrait_Set2_Tau.eigenD.txt > temp && mv temp ./output/ABCDS_multitrait_Set2_Tau.eigenD.txt

### process .fam phenotype in R 

## fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set2_Tau.fam")
## table(fam$V5)
## head(fam)
## residuals <- fread("./Data/Phenotype/ABCDS_multitrait_Set2_Tau_residuals.txt")
## residuals.fam <- merge(fam, residuals, by = c("V1","V2"))
## write.table(residuals.fam[,which(colnames(residuals.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set2_Tau_residuals.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

## phenos <- fread("./Data/Phenotype/ABCDS_multitrait_Set2_Tau_phenos.txt")
## phenos.fam <- merge(fam, phenos, by = c("V1","V2"))
## write.table(phenos.fam[,which(colnames(phenos.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set2_Tau_phenos.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

cp ./Data/Genotype/ABCDS_multitrait_Set2_Tau.bim ./Data/Genotype/ABCDS_multitrait_Set2_Tau_residuals.bim
cp ./Data/Genotype/ABCDS_multitrait_Set2_Tau.bed ./Data/Genotype/ABCDS_multitrait_Set2_Tau_residuals.bed

cp ./Data/Genotype/ABCDS_multitrait_Set2_Tau.bim ./Data/Genotype/ABCDS_multitrait_Set2_Tau_phenos.bim
cp ./Data/Genotype/ABCDS_multitrait_Set2_Tau.bed ./Data/Genotype/ABCDS_multitrait_Set2_Tau_phenos.bed

### generate BIMBAM format files 
# convert plink to VCF
plink --bfile ./Data/Genotype/ABCDS_multitrait_Set2_Tau_residuals --recode vcf --out ./Data/Genotype/ABCDS_multitrait_Set2_Tau_residuals
# convert VCF to BIMBAM
./qctool -g ./Data/Genotype/ABCDS_multitrait_Set2_Tau_residuals.vcf -ofiletype bimbam_dosage -og ./Data/Genotype/ABCDS_multitrait_Set2_Tau_residuals.txt

### Run Analysis 
### run with residuals
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set2_Tau_residuals -maf 0.02 -d ./output/ABCDS_multitrait_Set2_Tau.eigenD.txt -u ./output/ABCDS_multitrait_Set2_Tau.eigenU.txt -lmm 4 -n 1 2 3 4 -o ABCDS_multitrait_Set2_Tau_residuals
### run with raw phenotypes
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set2_Tau_phenos -maf 0.02 -d ./output/ABCDS_multitrait_Set2_Tau.eigenD.txt -u ./output/ABCDS_multitrait_Set2_Tau.eigenU.txt -n 1 2 3 4 -c ./Data/Phenotype/ABCDS_multitrait_Set2_Tau_covars.txt -lmm 4 -o ABCDS_multitrait_Set2_Tau_phenos

### ------------------------------------------------------------------------ ###
### Set 3: Tau group (N = 320)
### ------------------------------------------------------------------------ ###

### process genotype file 

# subset based on available phenotype
plink --bfile ./Data/Genotype/ABCDS_Cases_NonHis_WH --keep ./Data/Phenotype/ABCDS_multitrait_Set3_Tau_FID_IID.txt --allow-no-sex --biallelic-only --freq --make-bed --out ./Data/Genotype/ABCDS_multitrait_Set3_Tau

# change the phenotype column in fam to 1 for GEMMA to calculate kinship matrix 
awk '{$6 = ($6 == -9) ? 1 : $6}1' ./Data/Genotype/ABCDS_multitrait_Set3_Tau.fam > temp && mv temp ./Data/Genotype/ABCDS_multitrait_Set3_Tau.fam

# calculate relatedness matrix from genotype file 

/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set3_Tau -gk 2 -n 1 -o ABCDS_multitrait_Set3_Tau
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set3_Tau -k ./output/ABCDS_multitrait_Set3_Tau.sXX.txt -eigen -o ABCDS_multitrait_Set3_Tau

# change the 0 value to 0.00001 manually in ./output/ABCDS_multitrait_Set3_Tau.eigenD.txt
awk 'NR==1 && $1==0 { $1=0.0000001 } { print }' ./output/ABCDS_multitrait_Set3_Tau.eigenD.txt > temp && mv temp ./output/ABCDS_multitrait_Set3_Tau.eigenD.txt

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

cp ./Data/Genotype/ABCDS_multitrait_Set3_Tau.bim ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals.bim
cp ./Data/Genotype/ABCDS_multitrait_Set3_Tau.bed ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals.bed

cp ./Data/Genotype/ABCDS_multitrait_Set3_Tau.bim ./Data/Genotype/ABCDS_multitrait_Set3_Tau_phenos.bim
cp ./Data/Genotype/ABCDS_multitrait_Set3_Tau.bed ./Data/Genotype/ABCDS_multitrait_Set3_Tau_phenos.bed

### generate BIMBAM format files 
# convert plink to VCF
plink --bfile ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals --recode vcf --out ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals
# convert VCF to BIMBAM
./qctool -g ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals.vcf -ofiletype bimbam_dosage -og ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals.txt

### predict phenotype value 
# ./gemma -bfile [prefix] -k [filename] -predict -n [num1] [num2] [num3] -o [prefix]
cut -d $'\t' -f 3- ./Data/Phenotype/ABCDS_multitrait_Set3_Tau_residuals.txt > ./Data/Phenotype/ABCDS_multitrait_Set3_Tau_residuals_noID.txt
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -g ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals.txt -p ./Data/Phenotype/ABCDS_multitrait_Set3_Tau_residuals_noID.txt -n 1 2 3 4 -k ./output/ABCDS_multitrait_Set3_Tau.sXX.txt -predict -o ABCDS_multitrait_Set3_Tau_residuals_prdt

### update .fam phenotype in R again with predicted residuals

### update .fam again with predicted residuals 
## fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set3_Tau.fam")
## table(fam$V5)
## head(fam)
## residuals <- fread("./output/ABCDS_multitrait_Set3_Tau_residuals_prdt.prdt.txt")
## colnames(residuals) <- paste0("R",1:5)
## residuals.fam <- cbind(fam, residuals)
## colnames(residuals.fam)
## write.table(residuals.fam[,which(!(colnames(residuals.fam) %in% c("V6", "R5")))], "./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_prdt.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

cp ./Data/Genotype/ABCDS_multitrait_Set3_Tau.bim ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_prdt.bim
cp ./Data/Genotype/ABCDS_multitrait_Set3_Tau.bed ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_prdt.bed

### Run Analysis 
### run with predicted residuals
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_prdt -maf 0.01 -d ./output/ABCDS_multitrait_Set3_Tau.eigenD.txt -u ./output/ABCDS_multitrait_Set3_Tau.eigenU.txt -lmm 2 -n 1 2 3 4 -o ABCDS_multitrait_Set3_Tau_residuals_prdt

### ------------------------------------------------------------------------ ###
### Set 4: Abeta (no PET) N = 275
### ------------------------------------------------------------------------ ###

### process genotype file 

# subset based on available phenotype
plink --bfile ./Data/Genotype/ABCDS_Cases_NonHis_WH --keep ./Data/Phenotype/ABCDS_multitrait_Set4_AB_noPET_FID_IID.txt --allow-no-sex --biallelic-only --freq --make-bed --out ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET

# change the phenotype column in fam to 1 for GEMMA to calculate kinship matrix 
awk '{$6 = ($6 == -9) ? 1 : $6}1' ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET.fam > temp && mv temp ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET.fam

# calculate relatedness matrix from genotype file 

/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET -gk 2 -n 1 -o ABCDS_multitrait_Set4_AB_noPET
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET -k ./output/ABCDS_multitrait_Set4_AB_noPET.sXX.txt -eigen -o ABCDS_multitrait_Set4_AB_noPET

# change the 0 value to 0.00001 manually in ./output/ABCDS_multitrait_Set4_AB_noPET.eigenD.txt
awk 'NR==1 && $1==0 { $1=0.0000001 } { print }' ./output/ABCDS_multitrait_Set4_AB_noPET.eigenD.txt > temp && mv temp ./output/ABCDS_multitrait_Set4_AB_noPET.eigenD.txt

### process .fam phenotype in R 

## fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET.fam")
## table(fam$V5)
## head(fam)
## residuals <- fread("./Data/Phenotype/ABCDS_multitrait_Set4_AB_noPET_residuals.txt")
## residuals.fam <- merge(fam, residuals, by = c("V1","V2"))
## write.table(residuals.fam[,which(colnames(residuals.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_residuals.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

## phenos <- fread("./Data/Phenotype/ABCDS_multitrait_Set4_AB_noPET_phenos.txt")
## phenos.fam <- merge(fam, phenos, by = c("V1","V2"))
## write.table(phenos.fam[,which(colnames(phenos.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_phenos.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

cp ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET.bim ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_residuals.bim
cp ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET.bed ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_residuals.bed

cp ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET.bim ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_phenos.bim
cp ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET.bed ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_phenos.bed

### generate BIMBAM format files 
# convert plink to VCF
plink --bfile ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_residuals --recode vcf --out ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_residuals
# convert VCF to BIMBAM
./qctool -g ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_residuals.vcf -ofiletype bimbam_dosage -og ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_residuals.txt

### Run Analysis 
### run with residuals
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_residuals -maf 0.02 -d ./output/ABCDS_multitrait_Set4_AB_noPET.eigenD.txt -u ./output/ABCDS_multitrait_Set4_AB_noPET.eigenU.txt -lmm 4 -n 1 2 3 -o ABCDS_multitrait_Set4_AB_noPET_residuals

### ------------------------------------------------------------------------ ###
### Set 5: Tau (plus PET) N = 326
### ------------------------------------------------------------------------ ###

### process genotype file 

# subset based on available phenotype
plink --bfile ./Data/Genotype/ABCDS_Cases_NonHis_WH --keep ./Data/Phenotype/ABCDS_multitrait_Set5_Tau_plus_PET_FID_IID.txt --allow-no-sex --biallelic-only --freq --make-bed --out ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET

# change the phenotype column in fam to 1 for GEMMA to calculate kinship matrix 
awk '{$6 = ($6 == -9) ? 1 : $6}1' ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.fam > temp && mv temp ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.fam

# calculate relatedness matrix from genotype file 

/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET -gk 2 -n 1 -o ABCDS_multitrait_Set5_Tau_plus_PET
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET -k ./output/ABCDS_multitrait_Set5_Tau_plus_PET.sXX.txt -eigen -o ABCDS_multitrait_Set5_Tau_plus_PET

# change the 0 value to 0.00001 manually in ./output/ABCDS_multitrait_Set5_Tau_plus_PET.eigenD.txt
awk 'NR==1 && $1==0 { $1=0.0000001 } { print }' ./output/ABCDS_multitrait_Set5_Tau_plus_PET.eigenD.txt > temp && mv temp ./output/ABCDS_multitrait_Set5_Tau_plus_PET.eigenD.txt

### process .fam phenotype in R 

## fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.fam")
## table(fam$V5)
## head(fam)
## residuals <- fread("./Data/Phenotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals.txt")
## residuals.fam <- merge(fam, residuals, by = c("V1","V2"))
## write.table(residuals.fam[,which(colnames(residuals.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

## phenos <- fread("./Data/Phenotype/ABCDS_multitrait_Set5_Tau_plus_PET_phenos.txt")
## phenos.fam <- merge(fam, phenos, by = c("V1","V2"))
## write.table(phenos.fam[,which(colnames(phenos.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_phenos.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

cp ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.bim ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals.bim
cp ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.bed ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals.bed

cp ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.bim ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_phenos.bim
cp ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.bed ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_phenos.bed

### generate BIMBAM format files 
# convert plink to VCF
plink --bfile ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals --recode vcf --out ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals
# convert VCF to BIMBAM
./qctool -g ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals.vcf -ofiletype bimbam_dosage -og ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals.txt

### predict phenotype value 
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -g ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals.txt -p ./Data/Phenotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals_noID.txt -n 1 2 3 4 5 -k ./output/ABCDS_multitrait_Set5_Tau_plus_PET.sXX.txt -predict -o ABCDS_multitrait_Set5_Tau_plus_PET_residuals_prdt

### update .fam again with predicted residuals 
## fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.fam")
## table(fam$V5)
## head(fam)
## residuals <- fread("./output/ABCDS_multitrait_Set5_Tau_plus_PET_residuals_prdt.prdt.txt")
## colnames(residuals) <- paste0("R",1:6)
## residuals.fam <- cbind(fam, residuals)
## colnames(residuals.fam)
## write.table(residuals.fam[,which(!(colnames(residuals.fam) %in% c("V6", "R6")))], "./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals_prdt.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

cp ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.bim ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals_prdt.bim
cp ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.bed ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals_prdt.bed

### Run Analysis 
### run with predicted residuals
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals_prdt -maf 0.02 -d ./output/ABCDS_multitrait_Set5_Tau_plus_PET.eigenD.txt -u ./output/ABCDS_multitrait_Set5_Tau_plus_PET.eigenU.txt -lmm 2 -n 1 2 3 4 5 -o ABCDS_multitrait_Set5_Tau_plus_PET_residuals_prdt

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


