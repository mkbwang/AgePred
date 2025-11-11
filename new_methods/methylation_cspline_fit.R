

library(dplyr)
library(arrow)
library(deconvolution)
rm(list=ls())
# sum-to-one constraint

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


train_metadata <- read.csv(file.path(methylation_folder, "train", 
                                     sprintf("%s.csv", study)))

# round age
train_metadata$age <- round(train_metadata$age)
train_ages <- train_metadata$age
unique_train_ages <- unique(train_ages) |> sort()
x <- c(unique_train_ages[seq(1, length(unique_train_ages), 2)],
       min(unique_train_ages), max(unique_train_ages)) |> unique() # knots for fitting cubic splines


train_marker <- read_feather(file.path(methylation_folder, "train",
                                       sprintf("%s.feather", study)))
train_marker$Sample <- NULL


# fit cubic spline to all the marker values against time and apply F test
feature_cspline_summary <- data.frame(Feature=colnames(train_marker),
                                      IQR=0,
                                      EDF=0,
                                      Fstat=0,
                                      Pval=0)


for (i in 1:ncol(train_marker)){
  
  if (i %% 100 == 0){
    print(i)
  }  
  values <- train_marker[, i] |> as.matrix() |> as.vector()
  
  
  cspline_fit <- fit_cspline(t=train_ages, y=values, x=x)
  feature_cspline_summary$EDF[i] <- cspline_fit$EDF
  feature_cspline_summary$Fstat[i] <- cspline_fit$Fstat
  feature_cspline_summary$Pval[i] <- cspline_fit$pval
  
}

## calculate IQR
scale_estimates <- robust_scale(as.matrix(train_marker))
feature_cspline_summary$IQR <- scale_estimates$scale_vals


output_file <- sprintf("%s_feature.csv", study)
write.csv(feature_cspline_summary, file.path(method_folder, "cspline", output_file),
          row.names=FALSE, quote=FALSE)


