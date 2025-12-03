

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
studies_info <- studies_info %>% filter(Total >= 215 & Max_Age - Min_Age > 30)
studies <- studies_info$Study



library(optparse)
option_list <- list(make_option(c("-s", "--seed"), type="integer", default=1,
                                help="seed [default=1]"),
                    make_option(c("-u", "--subset", type="integer", default=40,
                                  help="subset [default=40]")))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
id <- as.integer(opt$seed)
subsample <- as.integer(opt$subset) # subsample size

study <- studies[id]

# load cspline results from fitting to all the training data
cubic_spline_summary <- read.csv(file.path(method_folder, "cspline", 
                                           sprintf("%s_logit_feature.csv", study)))
cubic_spline_summary <- cubic_spline_summary %>% arrange(Pval)
master_marker_features <- cubic_spline_summary$Feature[cubic_spline_summary$adjustedR2 < 0.05]

if (length(master_marker_features) > 50){
  master_marker_features <- master_marker_features[1:50]
}

# load training data
train_metadata_total <- read.csv(file.path(methylation_folder, "train",
                                     sprintf("%s.csv", study)),
                           row.names=1)

# round age
train_metadata_total$age <- round(train_metadata_total$age)
train_ages_total <- train_metadata_total$age
unique_train_ages <- unique(train_ages_total) |> sort()

train_marker_total <- read_feather(file.path(methylation_folder, "train",
                                       sprintf("%s.feather", study)))
train_marker_total$Sample <- NULL
train_marker_total <- as.matrix(train_marker_total)



test_metadata_total <- read.csv(file.path(methylation_folder, "test",
                                    sprintf("%s.csv", study)),
                          row.names=1)
test_metadata_total$age <- round(test_metadata_total$age)
test_ages_total <- test_metadata_total$age
test_marker_total <- read_feather(file.path(methylation_folder, "test",
                                      sprintf("%s.feather", study)))
test_marker_total$Sample <- NULL
test_marker_total <- as.matrix(test_marker_total)
test_marker_total[test_marker_total == 0] <- min(test_marker_total[test_marker_total > 0])/2
test_marker_total[test_marker_total == 1] <- (max(test_marker_total[test_marker_total < 1]) + 1)/2
test_marker_total <- log(test_marker_total) - log(1-test_marker_total)


test_maes <- data.frame(TrainSize=subsample,
                        MAE=rep(0, 20),
                        MAE_subset=rep(0, 20),
                        Success=TRUE)


for (j in 1:20){
  
  print(j)
  set.seed(j)
  subset_indices <- sample(1:nrow(train_marker_total), subsample)
  train_marker <- train_marker_total[subset_indices, ]
  train_ages <- train_ages_total[subset_indices]

  # fit cubic splines to decide top features
  IQR_medians <- robust_scale(train_marker, margin=2)
  zero_proportions <- colMeans(train_marker == 0)
  one_proportions <- colMeans(train_marker == 1)
  feature_mask <- (IQR_medians$scale_vals > 0.02) & (zero_proportions == 0) &
    (one_proportions == 0)
  
  
  train_marker <- train_marker[, feature_mask]
  logit_train_marker <- log(train_marker) - log(1 - train_marker)
  reference_ages <- seq(min(train_ages), max(train_ages), 1)
  spline_fit_logit_output <- tryCatch(
    expr = {
    fit_data(input_pheno=train_ages,
             input_marker_value = t(logit_train_marker),
             output_pheno=reference_ages)
    },
    error = function(e) {
      return(NULL)
    }
  )
  if (is.null(spline_fit_logit_output)){
    test_maes$Success[j] <- FALSE
    next
  }
  feature_logit_summary <- spline_fit_logit_output$spline_summary
  feature_logit_summary <- feature_logit_summary %>% arrange(Pval)
  feature_logit_summary$padj <- p.adjust(feature_logit_summary$Pval, method="BH")
  marker_features <- feature_logit_summary$Feature[feature_logit_summary$padj < 0.05]
  
  if (length(marker_features) > 50){
    marker_features <- marker_features[1:50]
  } else if (length(marker_features) < 30){
    marker_features <- feature_logit_summary$Feature[1:30]
  }
  
  ## take median of feature values at each age
  train_marker_subset <- as.data.frame(train_marker_total[subset_indices, marker_features])
  train_marker_subset$age <- train_ages
  train_marker_subset_median <- train_marker_subset %>% group_by(age) %>% 
    summarise_all(median) %>% arrange(age)
  age_labels <- train_marker_subset_median$age
  train_marker_subset_median$age <- NULL
  train_marker_subset_median <- as.matrix(train_marker_subset_median)
  train_marker_subset_median <- log(train_marker_subset_median) - log(1-train_marker_subset_median)
  
  test_marker_subset <- test_marker_total[, marker_features]
  
  prediction_mean <- rep(0, nrow(test_metadata_total))
  
  XXt_original <- as.matrix(train_marker_subset_median) %*% 
    t(as.matrix(train_marker_subset_median))
  trace_original <- diag(XXt_original) |> sum()
  
  
  for (i in 1:nrow(test_metadata_total)){
    
    test_data <- as.vector(as.matrix(test_marker_subset[i, ]))
    XY_original <- as.matrix(train_marker_subset_median) %*% test_data
    
    # deconvolution with raw methylation values
    fitted_result_original <- deconvolution(XXt=XXt_original, XY=XY_original,
                                            labels=age_labels, log=FALSE, standardize = FALSE,
                                            sum_constraint = T)
    prediction_mean[i] <- fitted_result_original$estim
    
  }
  
  test_maes$MAE[j] <- mean(abs(prediction_mean - test_ages_total))
  range_min <- min(train_ages)
  range_max <- max(train_ages)
  range_mask <- test_ages_total >= range_min & test_ages_total <= range_max
  test_maes$MAE_subset[j] <- mean(abs(prediction_mean[range_mask] - 
                                        test_ages_total[range_mask]))
  
}


ofile <- sprintf("methylation_%s_downsize_%d.csv", study, subsample)
write.csv(test_maes, file = file.path(method_folder, "deconvolution_downsize", ofile),
          row.names=FALSE, quote=FALSE)

