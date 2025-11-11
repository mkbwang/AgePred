

library(dplyr)
library(arrow)
library(deconvolution)
library(caret)
library(openxlsx)
rm(list=ls())
# sum-to-one constraint and smoothness penalty

geneexp_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/multitissue_transcriptome/"
method_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/new_methods/"


input_files <- list.files(file.path(geneexp_folder, "train"), pattern="*_counts.csv")
tissues <- sapply(input_files, function(fname){
  strsplit(fname, split="[_]")[[1]][1]
})

library(optparse)
option_list <- list(make_option(c("-s", "--seed"), type="integer", default=1,
                                help="seed [default=1]"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
id <- opt$seed
tissue <- tissues[id]


# load cspline result
cspline_result <- read.csv(file.path(method_folder, "cspline",
                                     sprintf("%s_feature.csv", tissue)))
cspline_result <- cspline_result %>% filter(IQR > 0.01)

cspline_result$padj <- p.adjust(cspline_result$Pval, method="BH")

cspline_result <- cspline_result %>% arrange(Pval)

if (sum(cspline_result$padj < 0.05) > 100){
  marker_feature <- cspline_result$Feature[cspline_result$padj < 0.05]
} else{
  marker_feature <- cspline_result$Feature[1:100]
}


train_metadata <- read.csv(file.path(geneexp_folder, "train",
                                     sprintf("metadata_%s.csv", tissue)))

# round age
train_metadata$Age <- round(train_metadata$Age)
train_ages <- train_metadata$Age
unique_train_ages <- unique(train_ages) |> sort()


train_marker <- read.csv(file.path(geneexp_folder, "train",
                                   sprintf("%s_counts.csv", tissue)),
                         row.names=1) |> t() |> as.data.frame()

train_marker <- train_marker[, marker_feature]

train_marker$age <- train_ages
train_marker_median <- train_marker %>% group_by(age) %>% summarise_all(median)
age_labels <- train_marker_median$age
train_marker_median$age <- NULL


test_metadata <- read.csv(file.path(geneexp_folder, "test",
                                    sprintf("metadata_%s.csv", tissue)))
test_metadata$Age <- round(test_metadata$Age)
test_marker <- read.csv(file.path(geneexp_folder, "test",
                                  sprintf("%s_counts.csv", tissue)),
                        row.names=1) |> t() |> as.data.frame()
test_marker <- test_marker[, marker_feature]


# normalize train and test data
scales_median <- robust_scale(train_marker_median)
train_marker_median <- train_marker_median[, scales_median$scale_vals > 0]
test_marker <- test_marker[, scales_median$scale_vals > 0]
normalized_data <- preprocess(X=train_marker_median, Y=test_marker, takelog=FALSE,
                              min_scale=0)

train_marker_median_normalized <- normalized_data$normalized_X
test_marker_normalized <- normalized_data$normalized_Y


weights_original <- matrix(0, nrow=nrow(test_marker), ncol=nrow(train_marker_median))
weights_normalized <- matrix(0, nrow=nrow(test_marker), ncol=nrow(train_marker_median))

rownames(weights_original) <- rownames(weights_normalized) <- test_metadata$SRR.ID
colnames(weights_original) <- colnames(weights_normalized) <- sprintf("Age_%d", age_labels)


test_metadata$Prediction_mean_original <- 0 # features are not normalized before deconvolution
test_metadata$Prediction_mode_original <- 0
test_metadata$lambda1_original <- 0 # optimal smoothness penalty parameter chosen by cross validation
test_metadata$lambda2_original <- 0

test_metadata$Prediction_mean_normalized <- 0 # features are normalized before deconvolution
test_metadata$Prediction_mode_normalized <- 0
test_metadata$lambda1_normalized <- 0 # optimal smoothness penalty parameter chosen by cross validation
test_metadata$lambda2_normalized <- 0
# to make sure that the smoothness penalty values are on the same scale as the

XXt_original <- as.matrix(train_marker_median) %*% t(as.matrix(train_marker_median))
trace_original <- diag(XXt_original) |> sum()
XXt_normalized <- as.matrix(train_marker_median_normalized) %*% t(as.matrix(train_marker_median_normalized))
trace_normalized <- diag(XXt_normalized) |> sum()

lambdas <- expand.grid(2^seq(-4, 4), 2^seq(-4, 4))
set.seed(2025)
for (i in 1:nrow(test_metadata)){

  print(i)
  test_data <- as.vector(as.matrix(test_marker[i, ]))
  XY_original <- as.matrix(train_marker_median) %*% test_data
  mean_XY_original <- abs(XY_original) |> mean()
  test_data_normalized <- as.vector(as.matrix(test_marker_normalized[i, ]))
  XY_normalized <- as.matrix(train_marker_median_normalized) %*% test_data_normalized
  mean_XY_normalized <- abs(XY_normalized) |> mean()

  mse_cv_original <- matrix(0, nrow=nrow(lambdas), ncol=5)
  mse_cv_normalized <- matrix(0, nrow=nrow(lambdas), ncol=5)

  foldids <- createFolds(seq(1, length(test_data)), k=5)
  # deconvolution with raw methylation values
  # five fold cross validation for selecting best lambda
  for (k in 1:5){

    XXt_original_subset <- as.matrix(train_marker_median[, -foldids[[k]]]) %*%
      t(as.matrix(train_marker_median[, -foldids[[k]]]))
    XY_original_subset <- as.matrix(train_marker_median[, -foldids[[k]]]) %*%
      test_data[-foldids[[k]]] |> as.vector()

    for (m in 1:nrow(lambdas)){
      fit_cv <- deconvolution(XXt=XXt_original_subset, XY=XY_original_subset,
                              labels=age_labels, log=FALSE, standardize = FALSE,
                              lambda1=lambdas[m, 1]*trace_original/6/length(age_labels),
                              lambda2=lambdas[m, 2]*mean_XY_original,
                              sum_constraint = F)
      pred_cv <- fit_cv$normalized_weights %*% as.matrix(train_marker_median[, foldids[[k]]])
      mse_cv_original[m, k] <- mean((test_data[foldids[[k]]] - pred_cv)^2)
    }
  }
  lambda_choice <- lambdas[which.min(rowSums(mse_cv_original)), ] |> as.matrix()
  test_metadata$lambda1_original[i] <- lambda_choice[1]
  test_metadata$lambda2_original[i] <- lambda_choice[2]

  fitted_result_original <- deconvolution(XXt=XXt_original, XY=XY_original,
                                          labels=age_labels, log=FALSE, standardize = FALSE,
                                          lambda1=lambda_choice[1]*trace_original/6/length(age_labels),
                                          lambda2=lambda_choice[2]*mean_XY_original,
                                          alpha=1,
                                          sum_constraint = F)
  test_metadata$Prediction_mean_original[i] <- fitted_result_original$estim
  weights_original[i, ] <- fitted_result_original$normalized_weights
  test_metadata$Prediction_mode_original[i] <- age_labels[which.max(fitted_result_original$normalized_weights)]


  # deconvolution with normalized methylation values
  # five fold cross validation for selecting best lambda
  for (k in 1:5){
    XXt_normalized_subset <- as.matrix(train_marker_median_normalized[, -foldids[[k]]]) %*%
      t(as.matrix(train_marker_median_normalized[, -foldids[[k]]]))
    XY_normalized_subset <- as.matrix(train_marker_median_normalized[, -foldids[[k]]]) %*%
      test_data_normalized[-foldids[[k]]] |> as.vector()
    for (m in 1:nrow(lambdas)){
      fit_cv <- deconvolution(XXt=XXt_normalized_subset, XY=XY_normalized_subset,
                              labels=age_labels, log=FALSE, standardize = FALSE,
                              lambda1=lambdas[m, 1]*trace_normalized/6/length(age_labels),
                              lambda2=lambdas[m, 2]*mean_XY_normalized,
                              sum_constraint = F)
      pred_cv <- fit_cv$normalized_weights %*% as.matrix(train_marker_median_normalized[, foldids[[k]]])
      mse_cv_normalized[m, k] <- mean((test_data_normalized[foldids[[k]]] - pred_cv)^2)
    }
  }
  lambda_choice <- lambdas[which.min(rowSums(mse_cv_normalized)), ] |> as.matrix()
  test_metadata$lambda1_normalized[i] <- lambda_choice[1]
  test_metadata$lambda2_normalized[i] <- lambda_choice[2]

  fitted_result_normalized <- deconvolution(XXt=XXt_normalized, XY=XY_normalized,
                                          labels=age_labels, log=FALSE, standardize = FALSE,
                                          lambda1=lambda_choice[1]*trace_normalized/6/length(age_labels),
                                          lambda2=lambda_choice[2]*mean_XY_normalized,
                                          sum_constraint = F)
  test_metadata$Prediction_mean_normalized[i] <- fitted_result_normalized$estim
  weights_normalized[i, ] <- fitted_result_normalized$normalized_weights
  test_metadata$Prediction_mode_normalized[i] <- age_labels[which.max(fitted_result_normalized$normalized_weights)]


}


#save to excel sheet
weights_original <- round(weights_original, digits = 6) |> as.data.frame()
weights_normalized <- round(weights_normalized, digits = 6) |> as.data.frame()

library(openxlsx)


output <- list(Prediction=test_metadata,
               Weights_original=weights_original,
               Weights_normalized=weights_normalized)

ofile <- sprintf("geneexp_%s_type3.xlsx", tissue)

write.xlsx(output, file = file.path(method_folder, "deconvolution", ofile),
           rowNames=TRUE)


