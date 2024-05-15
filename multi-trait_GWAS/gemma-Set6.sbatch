#!/bin/bash 
#SBATCH --job-name=Set6
#SBATCH --cpus-per-task=1     
#SBATCH --mail-user=rus39@pitt.edu    
#SBATCH --mail-type=ALL            
#SBATCH --mem-per-cpu=128GB   
#SBATCH --time=1-00:00:00 
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err   

module load plink/1.90b6.7 

cd /ix/kfan/Ruyu/multi-trait_GWAS

### run with residuals
/ix/kfan/Ruyu/multi-trait_GWAS/bin/gemma -bfile ./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET_residuals -maf 0.01 -d ./output/ABCDS_multitrait_set6_Tau_plus_PET.eigenD.txt -u ./output/ABCDS_multitrait_set6_Tau_plus_PET.eigenU.txt -lmm 2 -n 1 2 3 4 5 -o ABCDS_multitrait_set6_Tau_plus_PET_residuals

