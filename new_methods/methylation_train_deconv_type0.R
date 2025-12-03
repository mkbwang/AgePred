

library(dplyr)
library(arrow)
library(deconvolution)
rm(list=ls())
# sum-to-one constraint


# methylation_folder <- "extdata/methylation/" # TODO: folder to change
# method_folder <- "new_methods/"

methylation_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/methylation/" # TODO: folder to change
method_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/new_methods"


studies_info <- read.csv(file.path(methylation_folder, "studies_info.csv"))
studies_info <- studies_info %>% filter(Total >= 80 & Max_Age - Min_Age > 30)
studies <- studies_info$Study

library(optparse)
option_list <- list(make_option(c("-s", "--seed"), type="integer", default=1,
                                help="seed [default=1]"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
id <- opt$seed
study <- studies[id]


# load cspline result
cubic_spline_summary <- read.csv(file.path(method_folder, "cspline", 
                                           sprintf("%s_logit_feature.csv", study)))
marker_features <- cubic_spline_summary$Feature[cubic_spline_summary$isMarker]
if (length(marker_features) > 30){
  marker_features <- marker_features[1:30]
}


train_metadata <- read.csv(file.path(methylation_folder, "train",
                                     sprintf("%s.csv", study)),
                           row.names=1)

# round age
train_metadata$age <- round(train_metadata$age)
train_ages <- train_metadata$age
unique_train_ages <- unique(train_ages) |> sort()

train_marker <- read_feather(file.path(methylation_folder, "train",
                                       sprintf("%s.feather", study)))
train_marker$Sample <- NULL
train_marker <- train_marker[, marker_features]

train_marker$age <- train_ages
train_marker_median <- train_marker %>% group_by(age) %>% summarise_all(median)
age_labels <- train_marker_median$age
train_marker_median$age <- NULL
train_marker_median <- as.matrix(train_marker_median)
train_marker_median <- log(train_marker_median) - log(1-train_marker_median)


test_metadata <- read.csv(file.path(methylation_folder, "test",
                                    sprintf("%s.csv", study)),
                          row.names=1)
test_metadata$age <- round(test_metadata$age)
test_marker <- read_feather(file.path(methylation_folder, "test",
                                      sprintf("%s.feather", study)))
test_marker$Sample <- NULL
test_marker <- test_marker[, marker_features]
test_marker <- as.matrix(test_marker)
test_marker[test_marker == 0] <- min(test_marker[test_marker > 0])/2
test_marker[test_marker == 1] <- (max(test_marker[test_marker < 1]) + 1)/2
test_marker <- log(test_marker) - log(1-test_marker)

# normalize train and test data
scaling_params <- robust_scale(train_marker_median, margin=2)
train_marker_median_normalized <- scale_transform(value_mat=train_marker_median, 
                                           median_vals=scaling_params$median_vals,
                                           scale_vals = scaling_params$scale_vals,
                                           margin=2)
test_marker_normalized <- scale_transform(value_mat=test_marker, 
                                          median_vals=scaling_params$median_vals,
                                          scale_vals = scaling_params$scale_vals,
                                          margin=2)


weights_original <- matrix(0, nrow=nrow(test_marker), ncol=nrow(train_marker_median))
weights_normalized <- matrix(0, nrow=nrow(test_marker), ncol=nrow(train_marker_median))

rownames(weights_original) <- rownames(weights_normalized) <- rownames(test_metadata)
colnames(weights_original) <- colnames(weights_normalized) <- sprintf("Age_%d", age_labels)


test_metadata$Prediction_mean_original <- 0 # features are not normalized before deconvolution
test_metadata$Prediction_mode_original <- 0
test_metadata$Prediction_mean_normalized <- 0 # features are normalized before deconvolution
test_metadata$Prediction_mode_normalized <- 0


XXt_original <- as.matrix(train_marker_median) %*% t(as.matrix(train_marker_median))
trace_original <- diag(XXt_original) |> sum()
XXt_normalized <- as.matrix(train_marker_median_normalized) %*% t(as.matrix(train_marker_median_normalized))
trace_normalized <- diag(XXt_normalized) |> sum()


for (i in 1:nrow(test_metadata)){

  print(i)
  test_data <- as.vector(as.matrix(test_marker[i, ]))
  XY_original <- as.matrix(train_marker_median) %*% test_data
  test_data_normalized <- as.vector(as.matrix(test_marker_normalized[i, ]))
  XY_normalized <- as.matrix(train_marker_median_normalized) %*% test_data_normalized


  # deconvolution with raw methylation values
  fitted_result_original <- deconvolution(XXt=XXt_original, XY=XY_original,
                                          labels=age_labels, log=FALSE, standardize = FALSE,
                                          sum_constraint = F)
  test_metadata$Prediction_mean_original[i] <- fitted_result_original$estim
  weights_original[i, ] <- fitted_result_original$normalized_weights
  test_metadata$Prediction_mode_original[i] <- age_labels[which.max(fitted_result_original$normalized_weights)]


  # deconvolution with normalized methylation values
  fitted_result_normalized <- deconvolution(XXt=XXt_normalized, XY=XY_normalized,
                                          labels=age_labels, log=FALSE, standardize = FALSE,
                                          sum_constraint = F)
  test_metadata$Prediction_mean_normalized[i] <- fitted_result_normalized$estim
  weights_normalized[i, ] <- fitted_result_normalized$normalized_weights
  test_metadata$Prediction_mode_normalized[i] <- age_labels[which.max(fitted_result_normalized$normalized_weights)]


}


# save to excel sheet
weights_original <- round(weights_original, digits = 6) |> as.data.frame()
weights_normalized <- round(weights_normalized, digits = 6) |> as.data.frame()


library(openxlsx)


output <- list(Prediction=test_metadata,
               Weights_original=weights_original,
               Weights_normalized=weights_normalized)

ofile <- sprintf("methylation_%s_type0.xlsx", study)

write.xlsx(output, file = file.path(method_folder, "deconvolution", ofile),
           rowNames=TRUE)

