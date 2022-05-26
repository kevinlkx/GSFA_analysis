#!/bin/bash

#SBATCH --job-name=sceptre
#SBATCH --output=sceptre_%J.out
#SBATCH --error=sceptre_%J.err
#SBATCH --partition=broadwl
#SBATCH --account=pi-xinhe
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

module load R/4.0.4

cd /home/kaixuan/projects/GSFA_analysis/code

DATADIR=$1

Rscript sceptre_analysis.R \
  --datadir ${DATADIR} \
  --outdir ${DATADIR}/sceptre_output


