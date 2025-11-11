library(dplyr)
library(arrow)
library(deconvolution)
rm(list=ls())
# sum-to-one constraint

geneexp_folder <- "extdata/multitissue_transcriptome/"
method_folder <- "new_methods/"
output_folder <- "existing_methods/"

# geneexp_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/multitissue_transcriptome/"
# method_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/new_methods/"
# output_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/existing_methods/"


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

# if (sum(cspline_result$padj < 0.05) > 100){
#     marker_feature <- cspline_result$Feature[cspline_result$padj < 0.05]
# } else{
#     marker_feature <- cspline_result$Feature[1:100]
# }

# marker_feature <- cspline_result$Feature[cspline_result$Pval < 0.05]

marker_feature <- cspline_result$Feature# [cspline_result$Pval < 0.05]


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




test_metadata <- read.csv(file.path(geneexp_folder, "test",
                                    sprintf("metadata_%s.csv", tissue)))
test_metadata$Age <- round(test_metadata$Age)
test_marker <- read.csv(file.path(geneexp_folder, "test",
                                  sprintf("%s_counts.csv", tissue)),
                        row.names=1) |> t() |> as.data.frame()
test_marker <- test_marker[, marker_feature]




normalized_data <- preprocess(X=train_marker, Y=test_marker, takelog=FALSE,
                              min_scale=0)

train_marker_normalized <- normalized_data$normalized_X
test_marker_normalized <- normalized_data$normalized_Y



#### existing methods


# elastic net
library(glmnet)

en_fit <- cv.glmnet(x=train_marker_normalized, y=train_ages, alpha=0.5,
                    type.measure="mae")
en_coefs <- coef(en_fit, s="lambda.min")
en_coefs <- as.vector(en_coefs)


en_test_prediction <- predict(en_fit, newx=test_marker_normalized,
                              s="lambda.min")

test_metadata$EN <- as.vector(en_test_prediction)

# random forest
library(randomForest)
library(caret)

rf_fit <- randomForest(x=train_marker_normalized,
                       y=train_ages,
                       mtry=round(sqrt(ncol(train_marker_normalized))))
rf_test_prediction <- predict(rf_fit, newdata=test_marker_normalized)


test_metadata$RF <- rf_test_prediction

# lightGBM
library(lightgbm)

dtrain <- lgb.Dataset(data=as.matrix(train_marker_normalized),
                      label=train_ages)
params <- list(
    objective="regression",
    metric="mae",
    learning_rate=0.05,
    max_depth=4,
    num_leaves=30,
    feature_fraction=0.6,
    lambda_l1=0.4,
    lambda_l2=0.4
)

lgb_fit <- lgb.train(
    params=params,
    data=dtrain,
    nrounds=50
)


lgb_test_prediction <- predict(lgb_fit, newdata=test_marker_normalized)

test_metadata$LGB <- lgb_test_prediction

ofile <- sprintf("geneexp_%s.csv", tissue)

write.csv(test_metadata, file = file.path(output_folder, "output", ofile),
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
