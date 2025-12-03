

library(dplyr)
library(arrow)
library(ggplot2)
library(deconvolution)
rm(list=ls())
# sum-to-one constraint

methylation_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/methylation/"
method_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/new_methods/"

studies_info <- read.csv(file.path(methylation_folder, "studies_info.csv"))
studies_info <- studies_info %>% filter(Total >= 80 & Max_Age - Min_Age >= 30)
studies <- studies_info$Study

library(optparse)
option_list <- list(make_option(c("-s", "--seed"), type="integer", default=1, 
                                help="seed [default=1]"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
id <- opt$seed
study <- studies[id]


train_marker <- read_feather(file.path(methylation_folder, "train",
                                       sprintf("%s.feather", study))) |> as.data.frame()
rownames(train_marker) <- train_marker$Sample
train_marker$Sample <- NULL
train_marker_mat <- t(train_marker) |> as.matrix()

## calculate range of methylation and remove features that have zeros and/or small range
IQR_medians <- robust_scale(train_marker_mat, margin=1)
zero_proportions <- rowMeans(train_marker_mat == 0)
one_proportions <- rowMeans(train_marker_mat == 1)
feature_mask <- (IQR_medians$scale_vals > 0.02) & (zero_proportions == 0) &
  (one_proportions == 0)


train_metadata <- read.csv(file.path(methylation_folder, "train", 
                                     sprintf("%s.csv", study)),
                           row.names=1)

## remove samples in the largest or smallest age range if sample count too small and distant from other ages
train_metadata$age <- round(train_metadata$age)
train_ages <- train_metadata$age

unique_age_ascend <- unique(train_ages) |> sort(decreasing=FALSE)
unique_age_descend <- unique(train_ages) |> sort(decreasing=TRUE)

age_remove <- c()
for (i in 1:length(unique_age_ascend)){
  if(abs(unique_age_ascend[i+1]-unique_age_ascend[i]) >= 3 &
     sum(train_ages == unique_age_ascend[i]) <= 2){
    age_remove <- c(age_remove, unique_age_ascend[i])
  } else{
    break
  }
}

for (i in 1:length(unique_age_descend)){
  if(abs(unique_age_descend[i+1]-unique_age_descend[i]) >= 3 &
     sum(train_ages == unique_age_descend[i]) <= 2){
    age_remove <- c(age_remove, unique_age_descend[i])
  } else{
    break
  }
}

sample_mask <- !(train_ages %in% age_remove)

## take subset of markers and samples, then take logit of the methylation beta values
train_marker_mat <- train_marker_mat[feature_mask, sample_mask]
train_metadata <- train_metadata[sample_mask, ]
train_ages <- train_metadata$age
logit_train_marker_mat <- log(train_marker_mat) - log(1 - train_marker_mat)

reference_ages <- seq(min(train_metadata$age), max(train_metadata$age), 1)

# check that sample names match
stopifnot(all(colnames(train_marker_mat) == rownames(train_metadata)))

## fit splines on logit methylation level
begin <- proc.time()
spline_fit_logit_output <- fit_data(input_pheno=train_ages,
                                    input_marker_value = logit_train_marker_mat,
                                    output_pheno=reference_ages)
end <- proc.time()


feature_logit_summary <- spline_fit_logit_output$spline_summary
feature_logit_summary <- feature_logit_summary %>% arrange(desc(adjustedR2))
feature_logit_summary$padj <- p.adjust(feature_logit_summary$Pval, method="BH")
feature_names <- feature_logit_summary$Feature
spline_knots <- spline_fit_logit_output$spline_knots
predicted_mean_mat_logit <- spline_fit_logit_output$output_mean
predicted_se_mat_logit <- spline_fit_logit_output$output_se


## retain a subset of features
feature_logit_summary$isMarker <- feature_logit_summary$padj < 0.1 & 
  feature_logit_summary$EDF < length(spline_knots)*2/3


spline_prediction <- list(mean_logit=predicted_mean_mat_logit,
                          se_logit=predicted_se_mat_logit,
                          label_ages=reference_ages)


## output
output_file <- sprintf("%s_logit_feature.csv", study)
write.csv(feature_logit_summary, file.path(method_folder, "cspline", output_file),
          row.names=FALSE, quote=FALSE)


output_data <- sprintf("%s_logit_pred_data.rds", study)
saveRDS(spline_prediction, 
        file.path(method_folder, "cspline", output_data))


