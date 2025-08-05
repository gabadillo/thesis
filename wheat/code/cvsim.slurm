#!/bin/bash
#SBATCH --job-name=wheat  
#SBATCH --output=logs/wheat_%a.out  
#SBATCH --error=logs/wheat_%a.err   
#SBATCH --array=1-30            
#SBATCH --time=02:00:00                  
#SBATCH --mem=2G                        
#SBATCH --cpus-per-task=1                


module load R                   
Rscript 4.3.CV_simulation.R $SLURM_ARRAY_TASK_ID 
