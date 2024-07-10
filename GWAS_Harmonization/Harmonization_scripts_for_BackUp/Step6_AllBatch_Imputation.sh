#!/bin/bash

###--------------------------------------------------------------------------###
###
### This is a README file for Step6: 2nd Imputation for all batches
###
###	Several steps were mentioned in this step:
### (1) Perform LD pruning on autosomes
### (2) Pre-imputation check
### (3) Submit to the Michigan Impuation Server
### (4) Download the imputation outputs from the Michigan Impuation Server and annotate each chromosome
### (5) Concatenate each annotated chromosome into one file
### 
###																 Ruyu Shi
###															   05/26/2024
###--------------------------------------------------------------------------###

module load gcc/6.3.0 vcftools/0.1.16 plink/1.90b6.7 bcftools/1.15.1 
module load gcc/8.2.0 r/4.0.0 
module load snpeff/4.3t

cd pseudoGWAS_output 

###--------------------------------------------------------------------------###
### (1) Perform LD pruning on autosomes
###--------------------------------------------------------------------------###

### Michigan Imputation Server cannot handle huge size files, consider perform LD pruning to keep reasonable number of SNPs across the chromosome

plink --bfile All_Cohorts_allChr_cleaned --chr 1-22 --hwe 1E-05 --indep-pairphase 20000 200 0.9 --mind 0.05 --geno 0.05 --maf 0.01 --allow-extra-chr --biallelic-only --out All_Cohorts_allChr_cleaned_09

plink --bfile ./All_Cohorts_allChr_cleaned --allow-extra-chr --biallelic-only --extract ./All_Cohorts_allChr_cleaned_09.prune.in --freq --make-bed --out ./All_Cohorts_allChr_cleaned_09

###--------------------------------------------------------------------------###
### (2) Pre-imputation check
###--------------------------------------------------------------------------###

mkdir -p ../Michigan_Imputation_Server_LD09/
mkdir -p ../Michigan_Imputation_Server_ChrX_noPruning/

### extract chrX 
plink --bfile ./All_Cohorts_allChr_cleaned --chr 23 --allow-extra-chr --biallelic-only --freq --make-bed --out ./All_Cohorts_allChr_cleaned_chr23 

### perform pre-imputation check

### on autosomes

perl HRC-1000G-check-bim-v4.3.0/HRC-1000G-check-bim.pl -b ./All_Cohorts_allChr_cleaned_09.bim -f ./All_Cohorts_allChr_cleaned_09.frq -r /ix/kfan/Ruyu/harmonization/tools/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h 
sh Run-plink.sh

mv All_Cohorts_allChr_cleaned_09-updated-chr* ../Michigan_Imputation_Server_LD09/
cd ../Michigan_Imputation_Server_LD09

# generate sorted vcf.gz files
for chrom in $(seq 1 22); do

    input_file="All_Cohorts_allChr_cleaned_09-updated-chr${chrom}.vcf"
    output_file="${input_file}.gz"

    # sort and compress VCF
    bcftools sort "${input_file}" -Oz -o ../impu_prep/"${output_file}"
    # create index file
    bcftools index -t ../impu_prep/"${output_file}"

    echo "Sorting or compressing for ${input_file} finished."

done

### on chr23

perl HRC-1000G-check-bim-v4.3.0/HRC-1000G-check-bim.pl -b ./All_Cohorts_allChr_cleaned_chr23.bim -f ./All_Cohorts_allChr_cleaned_chr23.frq -r /ix/kfan/Ruyu/harmonization/tools/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h 
sh Run-plink.sh

mv All_Cohorts_allChr_cleaned_09-updated-chr23* ../Michigan_Imputation_Server_ChrX_noPruning/

input_file="All_Cohorts_allChr_cleaned_chr23-updated-chr23.vcf"
output_file="${input_file}.gz"

# sort and compress VCF
bcftools sort "${input_file}" -Oz -o ../impu_prep/"${output_file}"
# create index file
bcftools index -t ../impu_prep/"${output_file}"

echo "Sorting or compressing for ${input_file} finished."

###--------------------------------------------------------------------------###
### (3) Submit to the Michigan Impuation Server
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
# 	rsq Filter: 0.3
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
### (4) Download the imputation outputs from the Michigan Impuation Server
###--------------------------------------------------------------------------###

module load p7zip/16.02

### note: a total of 9 jobs were submitted to the Michigan Imputation Server 

### chr1

mkdir -p Chr1

cd ./Chr1

# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3512539/bc488b43b39bcd9444e2f319a9b89af32c0d4210a6f17aa4166615a6148cca7d | bash
7z e -p'8|K3ouDSogCsXk' chr_1.zip
# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3512540/e7b4304a5ff2e481b7f2760da6c300eb6cd729e6ed37971687966cbd25d91c8b | bash

bcftools annotate -x ID -Oz -o AllCohorts_Chr1_noID.vcf.gz ./chr1.dose.vcf.gz
tabix -p vcf AllCohorts_Chr1_noID.vcf.gz

# annotate rs number 
java -jar /ix/kfan/My_Bioinf_Apps/snpEff/SnpSift.jar annotate -id /ix/kfan/Ruyu/harmonization/tools/byChr/All_20180423.vcf.gz AllCohorts_Chr1_noID.vcf.gz > AllCohorts_Chr1_annotated.vcf
bgzip AllCohorts_Chr1_annotated.vcf

tabix -p vcf AllCohorts_Chr1_annotated.vcf.gz

### chr2

mkdir -p Chr2 

cd ./Chr2
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3512561/7457d70e93cacc3172373725453d27082d46a32af60fbcc6b4f6cee284393d46 | bash
7z e -p'wF7ACutb9DyB|X' chr_2.zip
# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3512562/f35d4561505fe80b505ae2a6d4df45da24566c5b407242e3b3827e86809da6ea | bash

bcftools annotate -x ID -Oz -o AllCohorts_Chr2_noID.vcf.gz ./chr2.dose.vcf.gz
tabix -p vcf AllCohorts_Chr2_noID.vcf.gz

# annotate rs number 
java -jar /ix/kfan/My_Bioinf_Apps/snpEff/SnpSift.jar annotate -id /ix/kfan/Ruyu/harmonization/tools/byChr/All_20180423.vcf.gz AllCohorts_Chr2_noID.vcf.gz > AllCohorts_Chr2_annotated.vcf
bgzip AllCohorts_Chr2_annotated.vcf

### chr3

mkdir -p Chr3

cd ./Chr3
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3512583/08f5374d4b4431b64e2d32de5f5f29fcc358466e4f5660475c13ede4310f5d46 | bash
7z e -p'QWf8D6Vke<pjCz' chr_3.zip
# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3512584/01256ced669be332bf1b43843dc79937142c4846f15316a2deac097df4005abc | bash

bcftools annotate -x ID -Oz -o AllCohorts_Chr3_noID.vcf.gz ./chr3.dose.vcf.gz
tabix -p vcf AllCohorts_Chr3_noID.vcf.gz

# annotate rs number 
java -jar /ix/kfan/My_Bioinf_Apps/snpEff/SnpSift.jar annotate -id /ix/kfan/Ruyu/harmonization/tools/byChr/All_20180423.vcf.gz AllCohorts_Chr3_noID.vcf.gz > AllCohorts_Chr3_annotated.vcf
bgzip AllCohorts_Chr3_annotated.vcf

### chr4
mkdir -p Chr4

cd ./Chr4
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3512605/597c8d18ec71e24140d2dc19719f9d226c9f7fdf3be17f4842eb40b30625c6ca | bash
7z e -p'mvZOGiP9ynMq9Q' chr_4.zip

curl -sL https://imputationserver.sph.umich.edu/get/3512606/965fbf085660dbb5d5c49fd7ebfa669b3cc548166823d2e2014963348c9acb20 | bash

bcftools annotate -x ID -Oz -o AllCohorts_Chr4_noID.vcf.gz ./chr4.dose.vcf.gz
tabix -p vcf AllCohorts_Chr4_noID.vcf.gz

# annotate rs number 
java -jar /ix/kfan/My_Bioinf_Apps/snpEff/SnpSift.jar annotate -id /ix/kfan/Ruyu/harmonization/tools/byChr/All_20180423.vcf.gz AllCohorts_Chr4_noID.vcf.gz > AllCohorts_Chr4_annotated.vcf
bgzip AllCohorts_Chr4_annotated.vcf

tabix -p vcf AllCohorts_Chr4_annotated.vcf.gz

### chr5
mkdir -p Chr5

cd ./Chr5
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3512627/48a7c180fe3ca78c4418a8cabbc39326379bd58a477bc4522caabbf7093faadf | bash
7z e -p'LGOjQoFJz48qiD' chr_5.zip

curl -sL https://imputationserver.sph.umich.edu/get/3512628/78c473098dc4999a5932dccb1f281cd79c7eaf8f86277bdea200095014319485 | bash

bcftools annotate -x ID -Oz -o AllCohorts_Chr15_noID.vcf.gz ./chr15.dose.vcf.gz
tabix -p vcf AllCohorts_Chr15_noID.vcf.gz

# annotate rs number 
java -jar /ix/kfan/My_Bioinf_Apps/snpEff/SnpSift.jar annotate -id /ix/kfan/Ruyu/harmonization/tools/byChr/All_20180423.vcf.gz AllCohorts_Chr15_noID.vcf.gz > AllCohorts_Chr15_annotated.vcf
bgzip AllCohorts_Chr15_annotated.vcf

tabix -p vcf AllCohorts_Chr15_annotated.vcf.gz

### chr6-9

mkdir -p Chr6-9

cd ./Chr6-9
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3512649/59b2dbb9a2d02d1c063aff62aa9b8caeedff9cb62e3d0684c71c82a3f41bd1d5 | bash

for file in chr_*.zip; do     
    7z e -p'Y9hNx0WajRlMmV' "$file"; 
done

# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3512650/85acd912d41bae3e7f2772f5bafd9758bed9ef1a7bd426e4b4268ada141ede1e | bash

for i in {6..9} 
do

	bcftools annotate -x ID -Oz -o AllCohorts_Chr${i}_noID.vcf.gz ./chr${i}.dose.vcf.gz
	tabix -p vcf AllCohorts_Chr${i}_noID.vcf.gz

	# annotate rs number 
	java -jar /ix/kfan/My_Bioinf_Apps/snpEff/SnpSift.jar annotate -id /ix/kfan/Ruyu/harmonization/tools/byChr/All_20180423.vcf.gz AllCohorts_Chr${i}_noID.vcf.gz > AllCohorts_Chr${i}_annotated.vcf
	bgzip AllCohorts_Chr${i}_annotated.vcf

done

### chr10-15

mkdir -p Chr10-15

cd ./Chr10-15
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3512671/de81ebb9ea276860f006535861d4458fd41677edc63bfce269dc9f3a9f737608 | bash
for file in chr_*.zip; do     
    7z e -p'T0kCAsUAifQx7G' "$file"; 
done
# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3512672/b7f7585dcc4e4bfeb2b58c52d4af0f5a093288005fadbc21fd3b205639d0bba8 | bash

for i in {10..15} 
do

	bcftools annotate -x ID -Oz -o AllCohorts_Chr${i}_noID.vcf.gz ./chr${i}.dose.vcf.gz
	tabix -p vcf AllCohorts_Chr${i}_noID.vcf.gz

	# annotate rs number 
	java -jar /ix/kfan/My_Bioinf_Apps/snpEff/SnpSift.jar annotate -id /ix/kfan/Ruyu/harmonization/tools/byChr/All_20180423.vcf.gz AllCohorts_Chr${i}_noID.vcf.gz > AllCohorts_Chr${i}_annotated.vcf
	bgzip AllCohorts_Chr${i}_annotated.vcf
	tabix -p vcf AllCohorts_Chr${i}_annotated.vcf.gz

done

### chr16-22

mkdir -p Chr16-22

cd ./Chr16-22
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3512693/5aee127634640c29b0b4f518043e4664fd191f9d7918ade12f39b417c0b3a425 | bash

for file in chr_*.zip; do     
    7z e -p'lRpRDtDWf8n2z9' "$file"; 
done

# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3512694/3eed3aabe5eec831673d643e4d2156b2fc9468182c5107fc2ed77990fcbffabc | bash

for i in {16..22} 
do

	bcftools annotate -x ID -Oz -o AllCohorts_Chr${i}_noID.vcf.gz ./chr${i}.dose.vcf.gz
	tabix -p vcf AllCohorts_Chr${i}_noID.vcf.gz

	# annotate rs number 
	java -jar /ix/kfan/My_Bioinf_Apps/snpEff/SnpSift.jar annotate -id /ix/kfan/Ruyu/harmonization/tools/byChr/All_20180423.vcf.gz AllCohorts_Chr${i}_noID.vcf.gz > AllCohorts_Chr${i}_annotated.vcf
	bgzip AllCohorts_Chr${i}_annotated.vcf
	tabix -p vcf AllCohorts_Chr${i}_annotated.vcf.gz

done

### chrX

cd /ix/kfan/Ruyu/harmonization/All_Cohorts/Michigan_Imputation_Server_ChrX_noPruning/
# QC Statistics
curl -sL https://imputationserver.sph.umich.edu/get/3511150/d39de8ccd4934306b537dd9144277916297bd03b0e2c267aaf2c4323d58b979e | bash

# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3511152/06fc65031891e099861e1c65095fd41fb56fd679b0c4b3d4f397346981849187 | bash
7z e -p'OFLQg54OkSt:jj' chr_X.zip

# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3511153/0d098b714d15a007a455153286e1c10262454cfd58c2a3558163c1bf05b78e93 | bash

bcftools annotate -x ID -Oz -o AllCohorts_ChrX_noID.vcf.gz ./chrX.dose.vcf.gz
tabix -p vcf AllCohorts_ChrX_noID.vcf.gz

# annotate rs number 
java -jar /ix/kfan/My_Bioinf_Apps/snpEff/SnpSift.jar annotate -id /ix/kfan/Ruyu/harmonization/tools/byChr/All_20180423.vcf.gz AllCohorts_ChrX_noID.vcf.gz > AllCohorts_ChrX_annotated.vcf
bgzip AllCohorts_ChrX_annotated.vcf

###--------------------------------------------------------------------------###
### (5) Concatenate each annotated chromosome into one file
###--------------------------------------------------------------------------###

echo './Michigan_Imputation_Server_LD09/Chr1/AllCohorts_Chr1_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr2/AllCohorts_Chr2_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr3/AllCohorts_Chr3_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr4/AllCohorts_Chr4_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr5/AllCohorts_Chr5_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr6-9/AllCohorts_Chr6_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr6-9/AllCohorts_Chr7_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr6-9/AllCohorts_Chr8_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr6-9/AllCohorts_Chr9_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr10-15/AllCohorts_Chr10_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr10-15/AllCohorts_Chr11_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr10-15/AllCohorts_Chr12_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr10-15/AllCohorts_Chr13_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr10-15/AllCohorts_Chr14_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr10-15/AllCohorts_Chr15_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr16-22/AllCohorts_Chr16_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr16-22/AllCohorts_Chr17_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr16-22/AllCohorts_Chr18_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr16-22/AllCohorts_Chr19_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr16-22/AllCohorts_Chr20_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr16-22/AllCohorts_Chr21_annotated.vcf.gz
./Michigan_Imputation_Server_LD09/Chr16-22/AllCohorts_Chr22_annotated.vcf.gz
./Michigan_Imputation_Server_ChrX_noPruning/AllCohorts_ChrX_noID.vcf.gz' > AllCohorts_AllChr_annotated_VCFs.txt

bcftools concat -f AllCohorts_AllChr_annotated_VCFs.txt -Oz -o AllCohorts_AllChr_annotated.vcf.gz
tabix -p vcf AllCohorts_AllChr_annotated.vcf.gz

# convert the vcf files to plink file
plink --vcf AllCohorts_AllChr_annotated.vcf.gz --make-bed --allow-extra-chr --biallelic-only --set-missing-var-ids @:#[b37]\$1,\$2 --out AllCohorts_AllChr_Annot




