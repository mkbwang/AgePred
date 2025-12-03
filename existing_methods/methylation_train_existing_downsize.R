library(dplyr)
library(arrow)
library(deconvolution)
rm(list=ls())
# sum-to-one constraint


# methylation_folder <- "extdata/methylation/"
# method_folder <- "new_methods/"
# output_folder <- "existing_methods/"

methylation_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/methylation/"
method_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/new_methods/"
output_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/existing_methods/"


studies_info <- read.csv(file.path(methylation_folder, "studies_info.csv"))
studies_info <- studies_info %>% filter(Total >= 215 & Max_Age - Min_Age > 30)
studies <- studies_info$Study

library(optparse)
option_list <- list(make_option(c("-s", "--seed"), type="integer", default=1,
                                help="seed [default=1]"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
id <- opt$seed
study <- studies[id]


# load cspline result
cspline_result <- read.csv(file.path(method_folder, "cspline",
                                     sprintf("%s_logit_feature.csv", study)))

marker_feature <- cspline_result$Feature


train_metadata <- read.csv(file.path(methylation_folder, "train",
                                     sprintf("%s.csv", study)),
                           row.names=1)

# round age
train_metadata$age <- round(train_metadata$age)
train_ages <- train_metadata$age

train_marker <- read_feather(file.path(methylation_folder, "train",
                                       sprintf("%s.feather", study)))
train_marker$Sample <- NULL
train_marker <- as.matrix(train_marker[, marker_feature])
train_marker <- log(train_marker) - log(1-train_marker)


test_metadata <- read.csv(file.path(methylation_folder, "test",
                                    sprintf("%s.csv", study)),
                          row.names=1)
test_metadata$age <- round(test_metadata$age)
test_marker <- read_feather(file.path(methylation_folder, "test",
                                      sprintf("%s.feather", study)))
test_marker$Sample <- NULL
test_marker <- as.matrix(test_marker[, marker_feature])
test_marker[test_marker == 0] <- min(test_marker[test_marker > 0])/2
test_marker[test_marker == 1] <- (max(test_marker[test_marker < 1]) + 1)/2
test_marker <- log(test_marker) - log(1-test_marker)





train_sizes <- c(40, 60, 80, 100, 120, 150)
random_seed <- seq(1, 20)

performance_result_list <- list()
for (i in 1:4){
  
  performance <- data.frame(SampleSize=train_sizes[i],
                            Seed=random_seed,
                            EN=0,
                            RF=0,
                            LGB=0)
  print(sprintf("Training data size %d", train_sizes[i]))
  for (j in 1:5){
    print(j)
    set.seed(j)
    subset_indices <- sample(1:nrow(train_metadata), train_sizes[i])
    train_ages_subset <- train_ages[subset_indices]
    train_marker_subset <- train_marker[subset_indices, ]
    scaling_params <- robust_scale(train_marker_subset, margin=2)
    train_marker_normalized_subset <- scale_transform(value_mat=train_marker_subset, 
                                               median_vals=scaling_params$median_vals,
                                               scale_vals = scaling_params$scale_vals,
                                               margin=2)
    test_marker_normalized <- scale_transform(value_mat=test_marker, 
                                              median_vals=scaling_params$median_vals,
                                              scale_vals = scaling_params$scale_vals,
                                              margin=2)
    
    # elastic net
    library(glmnet)
    
    en_fit <- cv.glmnet(x=train_marker_normalized_subset, y=train_ages_subset, alpha=0.5,
                        type.measure="mae")
    en_coefs <- coef(en_fit, s="lambda.min")
    en_coefs <- as.vector(en_coefs)
    
    en_test_prediction <- predict(en_fit, newx=test_marker_normalized,
                                  s="lambda.min")
    
    performance$EN[j] <- mean(abs(test_metadata$age - en_test_prediction))
    
    
    # random forest
    library(randomForest)
    library(caret)
    
    rf_fit <- randomForest(x=train_marker_subset,
                           y=train_ages_subset,
                           mtry=round(sqrt(ncol(train_marker))))
    rf_test_prediction <- predict(rf_fit, newdata=test_marker)
    
    performance$RF[j] <- mean(abs(test_metadata$age - rf_test_prediction))
    
    # lightGBM
    library(lightgbm)
    dtrain <- lgb.Dataset(data=as.matrix(train_marker_subset),
                          label=train_ages_subset)
    params <- list(
      objective="regression",
      metric="mae",
      learning_rate=0.05,
      max_depth=4,
      num_leaves=30,
      feature_fraction=0.6,
      lambda_l1=0.4,
      lambda_l2=0.4,
      verbose=-1
    )
    
    lgb_fit <-suppressWarnings(suppressMessages(lgb.train(
      params=params,
      data=dtrain,
      nrounds=50
    )))
    
    lgb_test_prediction <- predict(lgb_fit, newdata=test_marker)
    performance$LGB[j] <- mean(abs(test_metadata$age - lgb_test_prediction))
    
  }
  
  performance_result_list[[i]] <- performance
  
}


performance_result_combined <- do.call(rbind, performance_result_list)

#### existing methods


ofile <- sprintf("methylation_%s_train_downsize.csv", study)

write.csv(performance_result_combined, 
          file = file.path(output_folder, "output", ofile),
           row.names=TRUE)

# # support vector regression
# library(e1071)
#
# svr_fit <- svm(x=train_scaled,
#                y=train_ages_transformed,
#                type="eps-regression",
#                kernel="radial")
#
# svr_train_prediction <- predict(svr_fit, newdata=train_scaled)
# svr_train_prediction_raw <- age_anti_transform(svr_train_prediction)
# svr_test_prediction <- predict(svr_fit, newdata=test_scaled)
# svr_test_prediction_raw <- age_anti_transform(svr_test_prediction)
# svr_mae <- mean(abs(test_ages - svr_test_prediction_raw))
# svr_corr <- cor(test_ages, svr_test_prediction_raw)
#
#
# # xgboost
# library(xgboost)
#
# dtrain <- xgb.DMatrix(data=as.matrix(train_scaled),
#                       label=train_ages_transformed)
# dtest <- xgb.DMatrix(data=as.matrix(test_scaled),
#                      label=test_ages_transformed)
#
# params <- list(
#     objective="reg:absoluteerror",
#     eval_metric="mae",
#     eta=0.1,
#     gamma=0.425,
#     max_depth=3,
#     min_child_weight=9.5,
#     colsample_bytree=0.8,
#     alpha=0,
#     lambda=3.5
# )
#
# xgb_fit <- xgb.train(
#     params=params,
#     data=dtrain,
#     nrounds=300,
#     verbose=1
# )
#
# xgboost_train_prediction <- predict(xgb_fit, newdata=as.matrix(train_scaled))
# xgboost_train_prediction_raw <- age_anti_transform(xgboost_train_prediction)
# xgboost_test_prediction <- predict(xgb_fit, newdata=as.matrix(test_scaled))
# xgboost_test_prediction_raw <- age_anti_transform(xgboost_test_prediction)
# xgboost_mae <- mean(abs(test_ages - xgboost_test_prediction_raw))
# xgboost_corr <- cor(test_ages, xgboost_test_prediction_raw)
#
#
#
#
#
#
