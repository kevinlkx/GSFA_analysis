#!/bin/bash

#SBATCH --job-name=sceptre_Tcells
#SBATCH --output=sceptre_Tcells_permute_%J.out
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --account=pi-xinhe
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

module load R/4.2.0

cd /home/kaixuan/projects/GSFA_analysis/code

CONDITION=$1
PERMUTE_NUM=$2

OUTDIR=/project2/xinhe/kevinluo/GSFA/sceptre_analysis/Tcells_permuted_data
REGULARIZATION=0.1
B=1000
SEED=4

Rscript run_sceptre_Tcells_permuted_data.R \
  --outdir ${OUTDIR} \
  --condition ${CONDITION} \
  --regularization ${REGULARIZATION} \
  --B ${B} \
  --seed ${SEED} \
  --permute \
  --perm_num ${PERMUTE_NUM}




