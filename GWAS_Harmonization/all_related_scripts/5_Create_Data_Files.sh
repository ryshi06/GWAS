#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J GEM_create_data_files
#SBATCH --mem-per-cpu=64GB                                      # Memory per node
#SBATCH --cpus-per-task=1                                       # number of node
#SBATCH --mail-type=ALL               
#SBATCH --mail-user=rus39@pitt.edu

###--------------------------------------------------------------------------###

module load gcc/12.2.0 r/4.3.0 

cd /zfs1/kfan/Ruyu/harmonization_Sep5/scripts

###--------------------------------------------------------------------------###
### Run R script to generate genotype data files 
###--------------------------------------------------------------------------###

R CMD BATCH --vanilla 5_CreateDataFile.R 




