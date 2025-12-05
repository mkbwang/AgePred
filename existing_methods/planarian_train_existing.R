
library(dplyr)
library(deconvolution)
rm(list=ls())
library(glmnet)
library(randomForest)
library(lightgbm)


# load data and metadata
gene_tpm <- read.csv(file.path("planarian_data", "version2.csv"))
metadata <- read.csv(file.path("planarian_data", "meta_version2.csv"))


rownames(gene_tpm) <- gene_tpm$X
gene_tpm$X <- NULL
sample_names <- colnames(gene_tpm)
all_ages <- as.integer(metadata$chronological.age)[1:40]

test_metadata <- data.frame(Sample=sample_names[1:40],
                            Age=all_ages,
                            EN=0,
                            RF=0,
                            LGB=0)



# leave one sample out as test
for (loo_id in 1:40){

    train_ids <- setdiff(seq(1, 40), loo_id)

    train_samples <- sample_names[train_ids]
    test_sample <- sample_names[loo_id]


    gene_tpm_train <- gene_tpm[, train_samples] |> as.matrix()
    gene_tpm_test <- gene_tpm[, test_sample] |> as.vector()

    train_ages <- all_ages[train_ids]
    test_age <- all_ages[loo_id]

    # retain features that are present in all training samples
    presence_rate <- rowMeans(gene_tpm_train > 0)
    log_gene_tpm_train <- log(gene_tpm_train[presence_rate== 1, ])
    log_gene_tpm_test <- log(gene_tpm_test[presence_rate == 1])
    if (any(log_gene_tpm_test == -Inf)){
        log_gene_tpm_test[log_gene_tpm_test==-Inf] <-
            apply(log_gene_tpm_train[log_gene_tpm_test==-Inf, ,drop=FALSE], 1, function(myvec){
                min(myvec) - log(2)
            })
    }


    train_marker <- t(log_gene_tpm_train)
    test_marker <- log_gene_tpm_test

    # standardize all features
    scaling_params <- robust_scale(train_marker, margin=2)
    train_marker_normalized <- scale_transform(value_mat=train_marker,
                                               median_vals=scaling_params$median_vals,
                                               scale_vals = scaling_params$scale_vals,
                                               margin=2)
    test_marker_normalized <- scale_transform(value_mat=test_marker,
                                              median_vals=scaling_params$median_vals,
                                              scale_vals = scaling_params$scale_vals,
                                              margin=2)

    ## elastic net fit
    en_fit <- cv.glmnet(x=train_marker_normalized, y=train_ages, alpha=0.5,
                        type.measure="mae")
    en_coefs <- coef(en_fit, s="lambda.min")

    en_test_prediction <- predict(en_fit, newx=test_marker_normalized,
                                  s="lambda.min")
    test_metadata$EN[loo_id] <- en_test_prediction


    ## random forest fit
    rf_fit <- randomForest(x=train_marker_normalized,
                           y=train_ages,
                           mtry=round(sqrt(ncol(train_marker_normalized))))
    rf_test_prediction <- predict(rf_fit, newdata=test_marker_normalized)
    test_metadata$RF[loo_id] <- rf_test_prediction

    ## lightGBM
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


    lgb_test_prediction <- predict(lgb_fit, newdata=as.matrix(test_marker_normalized) |> t())
    test_metadata$LGB[loo_id] <- lgb_test_prediction

}


write.csv(test_metadata, "existing_methods/output/planarian_test_prediction.csv",
          row.names=FALSE, quote=FALSE)


