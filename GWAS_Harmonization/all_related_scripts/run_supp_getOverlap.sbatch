#!/bin/bash                                     
#SBATCH --job-name=GWAS_4_Overlap
#SBATCH --cpus-per-task=1     
#SBATCH --mail-user=rus39@pitt.edu    
#SBATCH --mail-type=ALL            
#SBATCH --time=1-00:00:00   
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

module load gcc/8.2.0 r/4.0.0 

Rscript /zfs1/kfan/Ruyu/harmonization_Sep5/scripts/supp_getOverlap_GWAS4.R