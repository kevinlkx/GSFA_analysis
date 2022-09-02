#!/bin/bash

#SBATCH --job-name=sceptre_LUHMES
#SBATCH --output=sceptre_LUHMES_runtime.out
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --account=pi-xinhe
#SBATCH --mem=10G

module load R/4.0.4

cd /home/kaixuan/projects/GSFA_analysis/code

Rscript run_sceptre_LUHMES_data.R \
  --datadir /project2/xinhe/kevinluo/GSFA/sceptre_analysis/LUHMES_data_runningtime \
  --outdir /project2/xinhe/kevinluo/GSFA/sceptre_analysis/LUHMES_data_runningtime/sceptre_output \
  --regularization 0.1 \
  --B 1000 \
  --seed 4



