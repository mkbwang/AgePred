#!/bin/sh

#SBATCH --job-name=methylation_deconvolution_downsize
#SBATCH --time=10:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-9
#SBATCH --mem=10g
#SBATCH --account=ligen0
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a-error.out


module load Rtidyverse/4.4.0

# Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 10
# 
# Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 20
# 
# Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 30
# 
# Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 40

# Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 50

# Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 20
# 
# Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 30
# 
# Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 40
# 
# Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 60
# 
Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 80
# 
Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 100
# 
Rscript --vanilla methylation_train_deconv_type1_downsize.R -s ${SLURM_ARRAY_TASK_ID} -u 120


# echo "Running deconvolution with sum to one constraint and smoothness penalty"
# Rscript --vanilla methylation_train_deconv_type2.R -s ${SLURM_ARRAY_TASK_ID}
# echo "Running deconvolution without sum to one constraint but with sparseness and smoothness penalty"
# Rscript --vanilla methylation_train_deconv_type3.R -s ${SLURM_ARRAY_TASK_ID}


