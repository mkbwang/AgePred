#!/bin/sh

#SBATCH --job-name=methylation_competitors
#SBATCH --time=00:40:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-18
#SBATCH --mem=15g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out


module load Rtidyverse/4.4.0

Rscript --vanilla methylation_train_existing.R -s ${SLURM_ARRAY_TASK_ID}

