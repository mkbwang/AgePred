

library(dplyr)
library(deconvolution)
rm(list=ls())
# sum-to-one constraint


geneexp_folder <- "extdata/multitissue_transcriptome/"
method_folder <- "new_methods/"

# geneexp_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/multitissue_transcriptome/"
# method_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/new_methods/"


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


train_metadata <- read.csv(file.path(geneexp_folder, "train",
                                     sprintf("metadata_%s.csv", tissue)))

# round age
train_metadata$Age <- round(train_metadata$Age)
train_ages <- train_metadata$Age
unique_train_ages <- unique(train_ages) |> sort()
x <- c(unique_train_ages[seq(1, length(unique_train_ages), 2)],
       min(unique_train_ages), max(unique_train_ages)) |> unique() # knots for fitting cubic splines


train_marker <- read.csv(file.path(geneexp_folder, "train",
                                       sprintf("%s_counts.csv", tissue)),
                         row.names=1) |> t()


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


output_file <- sprintf("%s_feature.csv", tissue)
write.csv(feature_cspline_summary, file.path(method_folder, "cspline", output_file),
          row.names=FALSE, quote=FALSE)


