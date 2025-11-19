#!/bin/sh

#SBATCH --job-name=methylation_deconvolution
#SBATCH --time=08:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-18
#SBATCH --mem=15g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out


module load Rtidyverse/4.4.0

echo "Running deconvolution with nonnegative constraint only"
Rscript --vanilla methylation_train_deconv_type0.R -s ${SLURM_ARRAY_TASK_ID}
echo "Running deconvolution with sum to one constraint only"
Rscript --vanilla methylation_train_deconv_type1.R -s ${SLURM_ARRAY_TASK_ID}
# echo "Running deconvolution with sum to one constraint and smoothness penalty"
# Rscript --vanilla methylation_train_deconv_type2.R -s ${SLURM_ARRAY_TASK_ID}
# echo "Running deconvolution without sum to one constraint but with sparseness and smoothness penalty"
# Rscript --vanilla methylation_train_deconv_type3.R -s ${SLURM_ARRAY_TASK_ID}

