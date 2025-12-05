

library(dplyr)
library(deconvolution)
library(caret)
library(patchwork)
rm(list=ls())
# no sum-to-one constraint


# methylation_folder <- "extdata/methylation/" # TODO: folder to change
# method_folder <- "new_methods/"

planarian_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/planarian_data" #TODO: change on greatlakes
method_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/new_methods" #TODO: change on greatlakes



library(optparse)
option_list <- list(make_option(c("-s", "--seed"), type="integer", default=1,
                                help="seed [default=1]"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
loo_id <- opt$seed


# load data and metadata
gene_tpm <- read.csv(file.path(planarian_folder, "version2.csv"))
metadata <- read.csv(file.path(planarian_folder, "meta_version2.csv"))


rownames(gene_tpm) <- gene_tpm$X
gene_tpm$X <- NULL
sample_names <- colnames(gene_tpm)


# leave one sample out as test
train_ids <- setdiff(seq(1, 40), loo_id)

train_samples <- sample_names[train_ids]
test_sample <- sample_names[loo_id]


gene_tpm_train <- gene_tpm[, train_samples] |> as.matrix()
gene_tpm_test <- gene_tpm[, test_sample] |> as.vector()
names(gene_tpm_test) <- rownames(gene_tpm_train)

train_ages <- as.integer(metadata$chronological.age[train_ids])
test_age <- as.integer(metadata$chronological.age[loo_id])




# load cspline result
cubic_spline_summary <- read.csv(file.path(method_folder, "cspline",
                                           sprintf("planarian_without_s%d.csv", loo_id)))
marker_features <- cubic_spline_summary$Gene[cubic_spline_summary$ismarker]



train_marker <- gene_tpm_train[marker_features, ] |> t() |> as.data.frame()

train_marker$age <- train_ages
train_marker_median <- train_marker %>% group_by(age) %>% summarise_all(median)
age_labels <- train_marker_median$age

train_marker_min <- train_marker %>% summarise_all(min)

train_marker_median$age <- NULL
train_marker_min$age <- NULL

train_marker$age <- NULL

train_marker_median <- as.matrix(train_marker_median)
train_marker_min <- as.matrix(train_marker_min)
train_marker_median <- log(train_marker_median)



test_marker <- gene_tpm_test[marker_features]
if (any(test_marker == 0)){
  test_marker[test_marker == 0] <- train_marker_min[test_marker == 0] / 2
}
test_marker <- log(test_marker)


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


weights_original_0 <- rep(0, nrow(train_marker_median))
weights_normalized_0 <- rep(0, nrow(train_marker_median))
weights_original_1 <- rep(0, nrow(train_marker_median))
weights_normalized_1 <- rep(0, nrow(train_marker_median))
weights_original_2 <- rep(0, nrow(train_marker_median))
weights_normalized_2 <- rep(0, nrow(train_marker_median))
weights_original_3 <- rep(0, nrow(train_marker_median))
weights_normalized_3 <- rep(0, nrow(train_marker_median))



names(weights_original_0) <- names(weights_normalized_0) <-
  names(weights_original_1) <- names(weights_normalized_1) <-
  names(weights_original_2) <- names(weights_normalized_2) <-
  names(weights_original_3) <- names(weights_normalized_3) <-
  sprintf("Age_%d", age_labels)


mean_original_0 <- mode_original_0 <- mean_normalized_0 <- mode_normalized_0 <-
  mean_original_1 <- mode_original_1 <- mean_normalized_1 <- mode_normalized_1 <-
  mean_original_2 <- mode_original_2 <- mean_normalized_2 <- mode_normalized_2 <-
  mean_original_3 <- mode_original_3 <- mean_normalized_3 <- mode_normalized_3 <- 0


# prepare to set up matrices for quadratic programming
XXt_original <- as.matrix(train_marker_median) %*% t(as.matrix(train_marker_median))
trace_original <- diag(XXt_original) |> sum()
XXt_normalized <- as.matrix(train_marker_median_normalized) %*% t(as.matrix(train_marker_median_normalized))
trace_normalized <- diag(XXt_normalized) |> sum()


XY_original <- as.matrix(train_marker_median) %*% test_marker
XY_normalized <- as.matrix(train_marker_median_normalized) %*% test_marker_normalized




##  0: no sum to one constraint, no smoothness or sparseness penalty
# deconvolution with raw values
fitted_result_original <- deconvolution(XXt=XXt_original, XY=XY_original,
                                        labels=age_labels, log=FALSE, standardize = FALSE,
                                        sum_constraint = FALSE)
mean_original_0 <- fitted_result_original$estim
weights_original_0 <- fitted_result_original$normalized_weights
mode_original_0 <- age_labels[which.max(weights_original_0)]
weight_plot_original_0 <- viz_weights(labels=age_labels, weights=weights_original_0,
                                      truth=test_age, title=sprintf("Type0, Raw Input \n(Mean %.2f, Mode %d)",
                                                                    mean_original_0, mode_original_0),
                                      size=7)



# deconvolution with normalized values
fitted_result_normalized <- deconvolution(XXt=XXt_normalized, XY=XY_normalized,
                                        labels=age_labels, log=FALSE, standardize = FALSE,
                                        sum_constraint = F)
mean_normalized_0 <- fitted_result_normalized$estim
weights_normalized_0 <- fitted_result_normalized$normalized_weights
mode_normalized_0 <- age_labels[which.max(weights_normalized_0)]
weight_plot_normalized_0 <- viz_weights(labels=age_labels, weights=weights_normalized_0,
                                      truth=test_age, title=sprintf("Type0, Standardized Input \n(Mean %.2f, Mode %d)",
                                                                    mean_normalized_0, mode_normalized_0),
                                      size=7)



##  1: sum to one constraint, no smoothness or sparseness penalty
# deconvolution with raw values
fitted_result_original <- deconvolution(XXt=XXt_original, XY=XY_original,
                                        labels=age_labels, log=FALSE, standardize = FALSE,
                                        sum_constraint = TRUE)
mean_original_1 <- fitted_result_original$estim
weights_original_1 <- fitted_result_original$normalized_weights
mode_original_1 <- age_labels[which.max(weights_original_1)]
weight_plot_original_1 <- viz_weights(labels=age_labels, weights=weights_original_1,
                                      truth=test_age, title=sprintf("Type1, Raw Input \n(Mean %.2f, Mode %d)",
                                                                    mean_original_1, mode_original_1),
                                      size=7)


# deconvolution with normalized values
fitted_result_normalized <- deconvolution(XXt=XXt_normalized, XY=XY_normalized,
                                          labels=age_labels, log=FALSE, standardize = FALSE,
                                          sum_constraint = TRUE)
mean_normalized_1 <- fitted_result_normalized$estim
weights_normalized_1 <- fitted_result_normalized$normalized_weights
mode_normalized_1 <- age_labels[which.max(weights_normalized_1)]
weight_plot_normalized_1 <- viz_weights(labels=age_labels, weights=weights_normalized_1,
                                      truth=test_age, title=sprintf("Type1, Standardized Input \n(Mean %.2f, Mode %d)",
                                                                    mean_normalized_1, mode_normalized_1),
                                      size=7)



##  2: sum to one constraint, smoothness penalty
lambdas <- 2^seq(-4, 4) # penalty range

mse_cv_original <- matrix(0, nrow=length(lambdas), ncol=5)

foldids <- createFolds(seq(1, length(test_marker)), k=5)
# deconvolution with raw values
# five fold cross validation for selecting best lambda
for (k in 1:5){

  XXt_original_subset <- as.matrix(train_marker_median[, -foldids[[k]]]) %*%
    t(as.matrix(train_marker_median[, -foldids[[k]]]))
  XY_original_subset <- as.matrix(train_marker_median[, -foldids[[k]]]) %*%
    test_marker[-foldids[[k]]] |> as.vector()

  for (m in 1:length(lambdas)){
    fit_cv <- deconvolution(XXt=XXt_original_subset, XY=XY_original_subset,
                            labels=age_labels, log=FALSE, standardize = FALSE,
                            lambda1=lambdas[m]*trace_original/6/length(age_labels),
                            sum_constraint = T)
    pred_cv <- fit_cv$normalized_weights %*% as.matrix(train_marker_median[, foldids[[k]]])
    mse_cv_original[m, k] <- mean((test_marker[foldids[[k]]] - pred_cv)^2)
  }
}
lambda1_2_original <- lambdas[which.min(rowSums(mse_cv_original))]

fitted_result_original <- deconvolution(XXt=XXt_original, XY=XY_original,
                                        labels=age_labels, log=FALSE, standardize = FALSE,
                                        lambda1=lambda1_2_original*trace_original/6/length(age_labels),
                                        sum_constraint = T)

mean_original_2 <- fitted_result_original$estim
weights_original_2 <- fitted_result_original$normalized_weights
mode_original_2 <- age_labels[which.max(weights_original_2)]
weight_plot_original_2 <- viz_weights(labels=age_labels, weights=weights_original_2,
                                      truth=test_age, title=sprintf("Type2, Raw Input \n(Mean %.2f, Mode %d)",
                                                                    mean_original_2, mode_original_2),
                                      size=7)


# deconvolution with normalized values
# five fold cross validation for selecting best lambda
mse_cv_normalized <- matrix(0, nrow=length(lambdas), ncol=5)
for (k in 1:5){
  XXt_normalized_subset <- as.matrix(train_marker_median_normalized[, -foldids[[k]]]) %*%
    t(as.matrix(train_marker_median_normalized[, -foldids[[k]]]))
  XY_normalized_subset <- as.matrix(train_marker_median_normalized[, -foldids[[k]]]) %*%
    test_marker_normalized[-foldids[[k]]] |> as.vector()
  for (m in 1:length(lambdas)){
    fit_cv <- deconvolution(XXt=XXt_normalized_subset, XY=XY_normalized_subset,
                            labels=age_labels, log=FALSE, standardize = FALSE,
                            lambda1=lambdas[m]*trace_normalized/6/length(age_labels),
                            sum_constraint = T)
    pred_cv <- fit_cv$normalized_weights %*% as.matrix(train_marker_median_normalized[, foldids[[k]]])
    mse_cv_normalized[m, k] <- mean((test_marker_normalized[foldids[[k]]] - pred_cv)^2)
  }
}
lambda1_2_normalized <- lambdas[which.min(rowSums(mse_cv_normalized))]

fitted_result_normalized <- deconvolution(XXt=XXt_normalized, XY=XY_normalized,
                                          labels=age_labels, log=FALSE, standardize = FALSE,
                                          lambda1=lambda1_2_normalized*trace_normalized/6/length(age_labels),
                                          sum_constraint = T)
mean_normalized_2 <- fitted_result_normalized$estim
weights_normalized_2 <- fitted_result_normalized$normalized_weights
mode_normalized_2 <- age_labels[which.max(fitted_result_normalized$normalized_weights)]
weight_plot_normalized_2 <- viz_weights(labels=age_labels, weights=weights_normalized_2,
                                      truth=test_age, title=sprintf("Type2, Standardized Input \n(Mean %.2f, Mode %d)",
                                                                    mean_normalized_2, mode_normalized_2),
                                      size=7)


##  3: no sum to one constraint, smoothness and sparseness penalty

lambdas <- expand.grid(2^seq(-4, 4), 2^seq(-4, 4))

mean_XY_original <- abs(XY_original) |> mean()
mean_XY_normalized <- abs(XY_normalized) |> mean()

mse_cv_original <- matrix(0, nrow=nrow(lambdas), ncol=5)
mse_cv_normalized <- matrix(0, nrow=nrow(lambdas), ncol=5)

foldids <- createFolds(seq(1, length(test_marker)), k=5)
# deconvolution with raw methylation values
# five fold cross validation for selecting best lambda
for (k in 1:5){

  XXt_original_subset <- as.matrix(train_marker_median[, -foldids[[k]]]) %*%
    t(as.matrix(train_marker_median[, -foldids[[k]]]))
  XY_original_subset <- as.matrix(train_marker_median[, -foldids[[k]]]) %*%
    test_marker[-foldids[[k]]] |> as.vector()

  for (m in 1:nrow(lambdas)){
    fit_cv <- deconvolution(XXt=XXt_original_subset, XY=XY_original_subset,
                            labels=age_labels, log=FALSE, standardize = FALSE,
                            lambda1=lambdas[m, 1]*trace_original/6/length(age_labels),
                            lambda2=lambdas[m, 2]*mean_XY_original,
                            sum_constraint = F)
    pred_cv <- fit_cv$weights %*% as.matrix(train_marker_median[, foldids[[k]]])
    mse_cv_original[m, k] <- mean((test_marker[foldids[[k]]] - pred_cv)^2)
  }

}
lambdas_3_original <- lambdas[which.min(rowSums(mse_cv_original)), ] |> as.matrix()

fitted_result_original <- deconvolution(XXt=XXt_original, XY=XY_original,
                                        labels=age_labels, log=FALSE, standardize = FALSE,
                                        lambda1=lambdas_3_original[1]*trace_original/6/length(age_labels),
                                        lambda2=lambdas_3_original[2]*mean_XY_original,
                                        alpha=1,
                                        sum_constraint = F)
mean_original_3 <- fitted_result_original$estim
weights_original_3 <- fitted_result_original$normalized_weights
mode_original_3 <- age_labels[which.max(weights_original_3)]
weight_plot_original_3 <- viz_weights(labels=age_labels, weights=weights_original_3,
                                        truth=test_age, title=sprintf("Type3, Raw Input \n(Mean %.2f, Mode %d)",
                                                                      mean_original_3, mode_original_3),
                                      size=7)



# deconvolution with normalized methylation values
# five fold cross validation for selecting best lambda
for (k in 1:5){
  XXt_normalized_subset <- as.matrix(train_marker_median_normalized[, -foldids[[k]]]) %*%
    t(as.matrix(train_marker_median_normalized[, -foldids[[k]]]))
  XY_normalized_subset <- as.matrix(train_marker_median_normalized[, -foldids[[k]]]) %*%
    test_marker_normalized[-foldids[[k]]] |> as.vector()
  for (m in 1:nrow(lambdas)){
    fit_cv <- deconvolution(XXt=XXt_normalized_subset, XY=XY_normalized_subset,
                            labels=age_labels, log=FALSE, standardize = FALSE,
                            lambda1=lambdas[m, 1]*trace_normalized/6/length(age_labels),
                            lambda2=lambdas[m, 2]*mean_XY_normalized,
                            sum_constraint = F)
    pred_cv <- fit_cv$weights %*% as.matrix(train_marker_median_normalized[, foldids[[k]]])
    mse_cv_normalized[m, k] <- mean((test_marker_normalized[foldids[[k]]] - pred_cv)^2)
  }
}

lambdas_3_normalized <- lambdas[which.min(rowSums(mse_cv_normalized)), ] |> as.matrix()

fitted_result_normalized <- deconvolution(XXt=XXt_normalized, XY=XY_normalized,
                                          labels=age_labels, log=FALSE, standardize = FALSE,
                                          lambda1=lambdas_3_normalized[1]*trace_normalized/6/length(age_labels),
                                          lambda2=lambdas_3_normalized[2]*mean_XY_normalized,
                                          sum_constraint = F)

mean_normalized_3 <- fitted_result_normalized$estim
weights_normalized_3 <- fitted_result_normalized$normalized_weights
mode_normalized_3 <- age_labels[which.max(weights_normalized_3)]
weight_plot_normalized_3 <- viz_weights(labels=age_labels, weights=weights_normalized_3,
                                      truth=test_age, title=sprintf("Type3, Standardized Input \n(Mean %.2f, Mode %d)",
                                                                    mean_normalized_3, mode_normalized_3),
                                      size=7)

combined_plots <- wrap_plots(weight_plot_original_0,  weight_plot_original_1, weight_plot_original_2, weight_plot_original_3,
                             weight_plot_normalized_0, weight_plot_normalized_1, weight_plot_normalized_2,  weight_plot_normalized_3,
                             ncol=4)


prediction_df <- data.frame(Truth=test_age,
  Type=rep(c("Type0", "Type1", "Type2", "Type3"), each=2),
                            Input=rep(c("Original", "Normalized"), 4),
                            pred_mean=c(mean_original_0, mean_normalized_0, mean_original_1, mean_normalized_1,
                                        mean_original_2, mean_normalized_2, mean_original_3, mean_normalized_3),
                            pred_mode=c(mode_original_0, mode_normalized_0, mode_original_1, mode_normalized_1,
                                        mode_original_2, mode_normalized_2, mode_original_3, mode_normalized_3),
                            lambda1=c(0,0,0,0, lambda1_2_original, lambda1_2_normalized, lambdas_3_original[1], lambdas_3_normalized[1]),
                            lambda2=c(0,0,0,0,0,0,lambdas_3_original[2], lambdas_3_normalized[2]))



ofile <- sprintf("planarian_%s.csv", test_sample)
ofigure1 <- sprintf("planarian_%s.png", test_sample)
ofigure2 <- sprintf("planarian_%s.svg", test_sample)

write.csv(prediction_df, file = file.path(method_folder, "deconvolution", "planarian", ofile),
          row.names=FALSE, quote=FALSE)

ggsave(file.path(method_folder, "deconvolution", "planarian", ofigure1),
       combined_plots, width=11, height=5)

ggsave(file.path(method_folder, "deconvolution", "planarian", ofigure2),
       combined_plots, width=11, height=5)
