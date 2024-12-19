#!/bin/bash
#SBATCH --job-name=RT1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=7000
#SBATCH --time=5-00:00:00


module load matlab/r2022a
matlab -nodisplay -nosplash -r "run('MAIN_MLP1_RLSSP.m');"
matlab -nodisplay -nosplash -r "run('compute_norms_RLSSP.m'); exit;"

    
squeue