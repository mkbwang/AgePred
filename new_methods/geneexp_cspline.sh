#!/bin/sh

#SBATCH --job-name=geneexp_cspline
#SBATCH --time=02:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-4
#SBATCH --mem=10g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out


module load Rtidyverse/4.4.0


Rscript --vanilla geneexp_cspline_fit.R -s ${SLURM_ARRAY_TASK_ID}


