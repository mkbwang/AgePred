#!/bin/sh

#SBATCH --job-name=methylation_MLE
#SBATCH --time=00:30:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-18
#SBATCH --mem=5g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out


module load Rtidyverse/4.4.0

echo "MLE estimates"
Rscript --vanilla methylation_spline_MLE_prediction.R -s ${SLURM_ARRAY_TASK_ID}

