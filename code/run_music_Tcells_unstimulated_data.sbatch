#!/bin/bash

#SBATCH --job-name=music
#SBATCH --output=music_Tcells_unstimulated_%J.out
#SBATCH --partition=broadwl
#SBATCH --account=pi-xinhe
#SBATCH --cpus-per-task=5
#SBATCH --mem=50G

module load R/4.0.4

cd /home/kaixuan/projects/GSFA_analysis/code

Rscript run_music_Tcells_unstimulated_data.R


