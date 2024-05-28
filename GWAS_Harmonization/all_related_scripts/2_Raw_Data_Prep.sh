#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J GEM_sample_prep
#SBATCH --mem-per-cpu=64GB                                      # Memory per node
#SBATCH --cpus-per-task=1                                       # number of node
#SBATCH --mail-type=ALL               
#SBATCH --mail-user=kaf166@pitt.edu

###--------------------------------------------------------------------------###
### This script is Transforming the raw data into a proper format 
###		for R/GWASTools to read each person's genotype accuractly 
###
### IMPORTANT!!
###		Please DO NOT run this script 
### 	This script should be only run one time, and never run again
###
###	Note:
###		In this script, all the comments are hashng with 3 #
###		and the commend will have only one #
###
###														Frank Fan
###														6/19/2019
###--------------------------------------------------------------------------###

cd /zfs1/kfan/GEM

###--------------------------------------------------------------------------###
### Transform the raw text file from Windows format to Linux format
###--------------------------------------------------------------------------###

### remove ^M in the Windows text format 
# tr -d '\r' <  GEM_GWAS_May2019_FinalReport.txt > tmp
# mv tmp GEM_GWAS_Linux_Jun2019_FinalReport.txt

### create a file with the list of sample IDs (scan names)
tail -n +11 GEM_GWAS_Linux_Jun2019_FinalReport.txt | awk '{print $2}' | uniq > GEM_ScanName.txt

### Save the file header into a separate file
head -n10 GEM_GWAS_Linux_Jun2019_FinalReport.txt > GEM_header.txt

###--------------------------------------------------------------------------###
### Breakdown huge file into several smaller sub-files
###--------------------------------------------------------------------------###

### Breakdown 2829 samples into 28 files with 100 smples each and the 29th file with 29 samples
### without header will be easier 

### First 28 files; 100 sample and 1748250*100 rows for each file 
for ind in {1..28}; do
	# Calculate the first line of the each file 
	((start=11+1748250*\($ind-1\)*100))
	# Subset the one huge into 28 smaller one
	tail -n +${start} GEM_GWAS_Linux_Jun2019_FinalReport.txt | head -n 174825000 > GEM_SubFile_${ind}.txt
done;


### The one last file, only has 29 samples
((last=11+1748250*2800))
tail -n +${last} GEM_GWAS_Linux_Jun2019_FinalReport.txt > GEM_SubFile_29.txt

###--------------------------------------------------------------------------###
### Generate genotype data for each individual 
###--------------------------------------------------------------------------###

### A nested for loop for the first 28 sub files & 2800 samples
### 

for file in {1..28}; do
	echo $file
	for rows in {1..100}; do
		
		((sample=\($file-1\)*100+$rows))
		scan=`sed -n ${sample}p GEM_ScanName.txt`
		
		((start=1+1748250*\($rows-1\)))
		tail -n +${start} GEM_SubFile_${file}.txt | head -n 1748250 > temp.txt
		cat GEM_header.txt temp.txt > tmp		
		awk 'BEGIN{OFS="\t"}{
			if(NR<10)
				print $0;
			else if(NR==10)
				print "snp", "sample", "X", "Y", "BAlleleFreq", "LogRRatio", "rawX", "rawY", "quality", "AlleleA", "AlleleB";
			else if ($3=="-" || $4=="-")
				print $1, $2, $5, $6, $7, $8, $9, $10, $11, "-", "-";
			else if($7>=0 && $7<=0.4 && $3==$4)
				print $1, $2, $5, $6, $7, $8, $9, $10, $11, "A", "A";
			else if($7>0.3 && $7<=0.7 && $3!=$4)
				print $1, $2, $5, $6, $7, $8, $9, $10, $11, "A", "B";
			else if($7>0.6 && $7<=1 && $3==$4)
				print $1, $2, $5, $6, $7, $8, $9, $10, $11, "B", "B";
			else if($7=="NaN")
				print $1, $2, $5, $6, $7, $8, $9, $10, $11, "-", "-";
			else
				print $1, $2, $5, $6, $7, $8, $9, $10, $11, "-", "-";
			}' tmp > GEM_Sample_${scan}.txt
			rm tmp temp.txt
	done;
done;


### The last 29 samples in the 29th subfile

for rows in {1..29}; do
	echo $rows
	((sample=2800+$rows))
	scan=`sed -n ${sample}p GEM_ScanName.txt`
	echo $scan
	((start=1+1748250*\($rows-1\)))
	tail -n +${start} GEM_SubFile_29.txt | head -n 1748250 > temp.txt
	cat GEM_header.txt temp.txt > tmp
	awk 'BEGIN{OFS="\t"}{
		if(NR<10)
			print $0;
		else if(NR==10)
			print "snp", "sample", "X", "Y", "BAlleleFreq", "LogRRatio", "rawX", "rawY", "quality", "AlleleA", "AlleleB";
		else if ($3=="-" || $4=="-")
			print $1, $2, $5, $6, $7, $8, $9, $10, $11, "-", "-";
		else if($7>=0 && $7<=0.4 && $3==$4)
			print $1, $2, $5, $6, $7, $8, $9, $10, $11, "A", "A";
		else if($7>0.3 && $7<=0.7 && $3!=$4)
			print $1, $2, $5, $6, $7, $8, $9, $10, $11, "A", "B";
		else if($7>0.6 && $7<=1 && $3==$4)
			print $1, $2, $5, $6, $7, $8, $9, $10, $11, "B", "B";
		else if($7=="NaN")
			print $1, $2, $5, $6, $7, $8, $9, $10, $11, "-", "-";
		else
			print $1, $2, $5, $6, $7, $8, $9, $10, $11, "-", "-";
		}' tmp > GEM_Sample_${scan}.txt
		rm tmp temp.txt
done;


### Remove all the sub-files and the header.txt

rm GEM_SubFile_*.txt GEM_header.txt






