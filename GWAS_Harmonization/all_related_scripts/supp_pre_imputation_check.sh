#!/bin/bash                                     

################################################################
### this file is used to format VCF files for Michigan Imputation Server
### 	HRC Pre-imputation Checks
###
### Note: 
### 	from https://imputationserver.readthedocs.io/en/latest/prepare-your-data/
################################################################

module load gcc/6.3.0 vcftools/0.1.16 plink/1.90b6.7 bcftools/1.15.1

################################
### Download tools and sites ###
################################

cd /zfs1/kfan/Ruyu/harmonization_Sep5/tools

wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

gzip -d HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

###############################################################
### Convert vcf or ped/map to bed & Create a frequency file ###
###############################################################

### GWAS1 
cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS1

vcftools --vcf GWAS1_plink_processed.vcf --out GWAS_1_raw --plink
plink --file GWAS_1_raw --freq --geno 0.1 --mind 0.1 --make-bed --out GWAS_1_raw_01

### GWAS3
# cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS3
# sed -n '161p;162q' GWAS_3_raw.vcf | cut -f1-8

vcftools --vcf GWAS_3_raw.vcf --out GWAS_3_raw --plink
plink --file GWAS_3_raw --freq --make-bed --out GWAS_3_raw

### GWAS4
# cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS4

### GWAS5
# cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS5

###########################################
### Execute pre-imputation check script ###
###########################################

# V4.3 version (the most up-to-date: chromosome X included)
perl ../tools/HRC-1000G-check-bim-v4.3.0/HRC-1000G-check-bim.pl -b GWAS_1_raw_01.bim -f GWAS_1_raw_01.frq -r ../tools/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

sh Run-plink.sh

##################################
### Create sorted vcf.gz files ###
##################################

# bcftools sort GWAS1_plink_processed.vcf -Oz -o GWAS1_plink_processed.vcf.gz
# bcftools index -t GWAS1_plink_processed.vcf.gz

mkdir -p ./impu_prep_11262023
for chrom in $(seq 1 23); do

    input_file="GWAS_1_raw_01-updated-chr${chrom}.vcf"
    output_file="${input_file}.gz"

    # sort and compress VCF
    bcftools sort "${input_file}" -Oz -o ./impu_prep_11262023/"${output_file}"
    # create index file
    bcftools index -t ./impu_prep_11262023/"${output_file}"

    echo "Sorting or compressing for ${input_file} finished."

done

###################################
### Download imputation results ###
###################################

module load p7zip/16.02

###----- GWAS 1
cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS1
mkdir -p ./impu_output_11272023
cd ./impu_output_11272023
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3343331/b90316f8aae9669fc4ccc8617c09c3e154665d5b8fb5d7a6db8d1a8eb1763dc1 | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3343333/b6c6da9bcfd429490a5ffac707d33da30eaa86cfbb35d4e617b6cc8a6672ccb3 | bash

for file in chr_*.zip; do     
    7z e -p9uF2HXvGAhvGms "$file"; 
done

# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3343334/35bf86f25657a2d98a6a335c06ab69184d65aa39a35dc7115bbe6a6086f9a029 | bash

###----- GWAS 1 with chrX
cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS1
mkdir -p ./impu_output_11272023_withX
cd ./impu_output_11272023_withX
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3344918/5af6ef4d79544362a2093fe07663ba64c9ee1717edacbbfb34a45e61f539549e | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3344920/6a74b6773b004b914f0da83a16d46ce042edfb2cedfc699099545f3b11fbcda8 | bash

for file in chr_*.zip; do     
    7z e -pmMZ{5YMO6Umeag "$file"; 
done

rm chr_*.zip

# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3344921/7aed9822c4600660d6386b4479fb5409d44b26ec93ae2013f7474672cb08f680 | bash

###----- GWAS 2
cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS2
mkdir -p ./impu_output_11272023
cd ./impu_output_11272023
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3342588/331e60176e9346e06b8d1a24f31c0c37b3e9882d1cf58b85f8f3d94bc9925265 | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3342590/7db7e7d618325eb8ff9831a3c614ed16a1bed606d2bcfaa185e7e623009e86af | bash

for file in chr_*.zip; do     
    7z e -pJQaFddn2Vc4Xv2 "$file"; 
done

# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3342591/f8664b830a7cd98314eb84e3fbccfb512c0680fcd44903136554bde634b28562 | bash

###----- GWAS 2 with chrX
cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS2
mkdir -p ./impu_output_11272023_withX
cd ./impu_output_11272023_withX
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3344852/ad77efd64b7645d9ebcce224948b054474ea90f5306ccb009053954d63c0fa8a | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3344854/7628121a55952b76935ee51b75879300e8d46266d27cc93ef631a6cc8e45d4fb | bash

for file in chr_*.zip; do     
    7z e -pGqUPF8Bvt4WHnl "$file"; 
done

rm chr_*.zip

# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3344855/4ec6ca87817d41a060ff494a29cd55c844e6dd5cca22fa9e68731d925437bfa6 | bash

# updated chrX from GWAS 2
cd /ix/kfan/Ruyu/harmonization/GWAS2
mkdir -p chrX

cd ./chrX
# QC statistics 
curl -sL https://imputationserver.sph.umich.edu/get/3473087/815350e5aa3082daa5d01194c20b309abc4e75d921925a88eade9cb853ac2ef0 | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3473089/e466bcb8b1233de69eb8f8287a17737c67de499c38a40ce65d3195671f4e4fe1 | bash

for file in chr_*.zip; do     
    7z e -p'Tcz4cBdCMP2nU(' "$file"; 
done
# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3473090/4f92925a0f062cc0e29df8a6ac1c29e8e7ac0fa4aa6b7421532dea0651254719 | bash

###----- GWAS 3
cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS3
mkdir -p ./impu_output_11272023
cd ./impu_output_11272023
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3342610/524601951ade3d5ea206a351d7e60934ae4342611586d32ad20de311d1e9f41a | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3342612/85c2196f09ae74dbeb21644530acf50c43e90b548d6696ae020ad1c14544dcfa | bash

for file in chr_*.zip; do     
    7z e -pgfUWhAz8M1CyzQ "$file"; 
done

# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3342613/9b878e6f5a7d1fd76de90aa91be9d345bdcd63ed28fe48df30e92ff96557e7c8 | bash

###----- GWAS 3 with chrX
cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS3
mkdir -p ./impu_output_11272023_withX
cd ./impu_output_11272023_withX
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3344874/5ec606866b6f65e2c8d5c3e5b97bedd846be0e7e9c3b6ae999ad04f2509dacd8 | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3344876/84cc98f75e315e145dee66e64e6132bb1b17a95a00d8521245b6ebe51b822fb7 | bash

for file in chr_*.zip; do     
    7z e -pkoINRozgRyM08Z "$file"; 
done

rm chr_*.zip

# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3344877/0c6543b55525c376ae7ea069ef487ffcdad49ac245b28c076deb4831f5fa334e | bash

###----- GWAS 4
cd /ix/kfan/Ruyu/harmonization/GWAS4
mkdir -p ./impu_output_withX

cd ./impu_output_withX
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3472450/c7c8f3c6e49bc08e235a8f88864e4a44355d037f1a7284ca6f77e75b0e9a69d8 | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3472452/97a34f80df1734095d49447bdfebfac34fbc71c18063ca84d0d98e5857e1b9f6 | bash

for file in chr_*.zip; do     
    7z e -pjRYxg9SSA7mw "$file"; 
    rm "$file";
done

###----- GWAS 5
cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS5
mkdir -p ./impu_output_11272023
cd ./impu_output_11272023
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3344786/420dc17370cfb6124a2a9a773810f6800caccf1320a1c15541b59640d61209c7 | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3344788/c5981c61874dcfbc296e1cdb6606edc6876259492753f5a45193f3c38e207d27 | bash

for file in chr_*.zip; do     
    7z e -pb6WMcIFu8OH0mp "$file"; 
done

rm chr_*.zip

# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3344789/57681edfe8b45abfcd16e23958f494562b09d87f6098ff02c433f2c9baeb39c1 | bash

###----- GWAS 5 with chrX
cd /zfs1/kfan/Ruyu/harmonization_Sep5/GWAS5
mkdir -p ./impu_output_11272023_withX
cd ./impu_output_11272023_withX
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3344830/011147b1cae693659a37b412896924db2359923a8d2d7af8322cf216d3bcc566 | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3344832/d148d8d65a844fbc8323b210445bfec55b1e503d72afeefce4a08b7b9474d39f | bash

for file in chr_*.zip; do     
    7z e -pJe6q_1kGcSsXmP "$file"; 
    rm "$file";
done

# rm chr_*.zip
 
# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3344833/83857f117db9c0f8dc634cd0d8fbb02d103a4abb629c51c872f74b7406ce7107 | bash

###----- GEM
cd /zfs1/kfan/Ruyu/harmonization_Sep5/GEM
mkdir -p ./impu_output_11272023
cd ./impu_output_11272023
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3344830/011147b1cae693659a37b412896924db2359923a8d2d7af8322cf216d3bcc566 | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3344832/d148d8d65a844fbc8323b210445bfec55b1e503d72afeefce4a08b7b9474d39f | bash

for file in chr_*.zip; do     
    7z e -pJe6q_1kGcSsXmP "$file"; 
done

rm chr_*.zip

# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3344833/83857f117db9c0f8dc634cd0d8fbb02d103a4abb629c51c872f74b7406ce7107 | bash

###----- GEM with chrX
cd /zfs1/kfan/Ruyu/harmonization_Sep5/GEM
mkdir -p ./impu_output_11272023_withX
cd ./impu_output_11272023_withX
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3345897/5e3185728dfae283e4805ac158746fa2aeea6505f0f19641dbf39e28df390c70 | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3345899/4c13dbf9703ede4dde73bb836e1d6083dfac3644f8fd195ed86fb9e64c023bcc | bash

for file in chr_*.zip; do     
    7z e -p'SVA>4qeHt9zIsm' "$file"; 
    rm "$file";
done

# rm chr_*.zip
 
# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3345900/f986ad2e02430eb13b8ce380179b1b0f3b66c375df7944637784bf02df2004b2 | bash

###----- SLE with chrX
cd /ix/kfan/Ruyu/harmonization/SLE
mkdir -p ./impu_output_withX
cd ./impu_output_withX
# QC statistics
curl -sL https://imputationserver.sph.umich.edu/get/3476400/90e38baaf5dd030bd20d1cf2c3a3dfe950cf09a1bd9af075e1f3e52a7c91d244 | bash
# Imputation Results
curl -sL https://imputationserver.sph.umich.edu/get/3476402/2d3f4a75203f0e80fdbb443e7adb27223d2fa79f2d29368be342bf8de5895148 | bash

for file in chr_*.zip; do     
    7z e -p'5g7eCQjTgz<5BX' "$file"; 
    rm "$file";
done

# rm chr_*.zip
 
# Logs
curl -sL https://imputationserver.sph.umich.edu/get/3476403/1fcf635febd5dab34e8532c21762458e5b3f39623ee96cc15f2faaa7daa566dc | bash


