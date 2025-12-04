
rm(list=ls())
library(deconvolution)
library(mgcv)
library(dplyr)
library(ggplot2)


planarian_folder <- "planarian_data" #TODO: change on greatlakes
method_folder <- "new_methods" #TODO: change on greatlakes





library(optparse)
option_list <- list(make_option(c("-s", "--seed"), type="integer", default=1,
                                help="seed [default=1]"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
loo_id  <- opt$seed


# load data and metadata
gene_tpm <- read.csv(file.path("planarian_data", "version2.csv"))
metadata <- read.csv(file.path("planarian_data", "meta_version2.csv"))


rownames(gene_tpm) <- gene_tpm$X
gene_tpm$X <- NULL
sample_names <- colnames(gene_tpm)


# leave one sample out as test
loo_id <- 25 # TODO: remove when on great lakes
train_ids <- setdiff(seq(1, 40), loo_id)

train_samples <- sample_names[train_ids]
test_sample <- sample_names[loo_id]


gene_tpm_train <- gene_tpm[, train_samples] |> as.matrix()
gene_tpm_test <- gene_tpm[, test_sample] |> as.vector()


train_ages <- as.integer(metadata$chronological.age[train_ids])
test_age <- as.integer(metadata$chronological.age[loo_id])


# only keep genes that are nonzero in all the training samples
presence_rate <- rowMeans(gene_tpm_train > 0)
log_gene_tpm_train <- log(gene_tpm_train[presence_rate== 1, ])
log_gene_tpm_test <- log(gene_tpm_test[presence_rate == 1])
log_gene_tpm_test[log_gene_tpm_test==-Inf] <-
    apply(log_gene_tpm_train[log_gene_tpm_test==-Inf, ], 1, function(myvec){
        min(myvec) - log(2)
    })

gene_names <- rownames(log_gene_tpm_train)



# load utility functions
source(file.path(planarian_folder, "planarian_inference_utils.R"))


# fit between age 6-23
spline_6_23 <- spline_fit(ages=train_ages, expression_mat=log_gene_tpm_train,
                          min_age=6, max_age=23)
summary_6_23 <- spline_6_23$fit_summary
summary_6_23$ismarker <- summary_6_23$adjustedR2 > 0.6


write.csv(summary_6_23, file.path(method_folder, "cspline",
                                  sprintf("planarian_without_s%d.csv", loo_id)))



