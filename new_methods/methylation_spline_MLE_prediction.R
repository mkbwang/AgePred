

library(dplyr)
library(arrow)
library(deconvolution)
library(openxlsx)
rm(list=ls())

methylation_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/methylation/"
method_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/new_methods/"

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


## retain features with high adjustedR2, small EDF and small p value based on cubic spline fit
cubic_spline_summary <- read.csv(file.path(method_folder, "cspline", 
                                           sprintf("%s_logit_feature.csv", study)))
marker_features <- cubic_spline_summary$Feature[cubic_spline_summary$isMarker]
feature_significance <- cubic_spline_summary$adjustedR2[cubic_spline_summary$isMarker]
if (length(marker_features) > 30){
  marker_features <- marker_features[1:30]
  feature_significance <- feature_significance[1:30]
}


## load predicted values from cspline
cspline_prediction <- readRDS(file.path(method_folder, "cspline",
                                        sprintf("%s_logit_pred_data.rds", study)))
predicted_mean_mat_logit <- cspline_prediction$mean_logit
predicted_se_mat_logit <- cspline_prediction$se_logit
reference_ages <- cspline_prediction$label_ages

logit_reference_mean <- predicted_mean_mat_logit[marker_features, ]
logit_reference_se <- predicted_se_mat_logit[marker_features, ]



# load test data
test_marker <- read_feather(file.path(methylation_folder, "test", 
                                       sprintf("%s.feather", study))) |> as.data.frame()
test_metadata <- read.csv(file.path(methylation_folder, "test", 
                                    sprintf("%s.csv", study)), row.names=1)

rownames(test_marker) <- test_marker$Sample
test_marker$Sample <- NULL
test_marker_mat <- t(test_marker) |> as.matrix()
test_marker_mat[test_marker_mat == 0] <- min(test_marker_mat[test_marker_mat > 0])/2
test_marker_mat[test_marker_mat == 1] <- (max(test_marker_mat[test_marker_mat < 1]) + 1)/2
logit_test_marker_mat <- log(test_marker_mat) - log(1-test_marker_mat)

logit_test_marker_mat <- logit_test_marker_mat[marker_features, ]


## aggregate MLE curves from different features
predicted_age_mle <- rep(0, nrow(test_metadata))
consensus_lik_mat <- matrix(0, nrow=nrow(test_metadata),
                            ncol=length(reference_ages))

for (j in 1:nrow(test_metadata)){

    values_vec <- logit_test_marker_mat[, j]

    # with a flat prior
    lik_result <- likelihood(mean_mat=logit_reference_mean,
                             se_mat=logit_reference_se,
                             y=values_vec,
                             exclude_factor=1)

    lik_mat <- lik_result$lik
    # plot(reference_ages, lik_mat[28, ])
    validity_mask <- lik_result$mask

    mle_result <- mle(lik=lik_mat,
                      feature_significance=feature_significance,
                      mask=validity_mask,
                      labels=reference_ages)


    consensus_lik <- mle_result$lik
    consensus_lik_mat[j, ] <- consensus_lik$lik
    predicted_age <- mle_result$estimate


    # plot(reference_ages, consensus_lik$lik)
    predicted_age_mle[j] <- predicted_age

}

test_metadata$mle_pred <- predicted_age_mle
rownames(consensus_lik_mat) <- rownames(test_metadata)
colnames(consensus_lik_mat) <- sprintf("Age%d", reference_ages)


output <- list(Prediction=test_metadata,
               likelihood=as.data.frame(consensus_lik_mat))

ofile <- sprintf("methylation_%s.xlsx", study)

write.xlsx(output, file = file.path(method_folder, "MLE", ofile),
           rowNames=TRUE)


