
library(dplyr)
library(tidyr)
library(ggplot2)
rm(list=ls())
train_sizes <- c(10, 20, 30, 40, 60, 80, 100, 120)

studies_info <- read.csv("extdata/methylation/studies_info_subset.csv")
studies_info <- studies_info %>% filter(Total > 215)

studies <- studies_info$Study

deconv_results <- list()

for (j in 1:length(studies)){

    current_study <- studies[j]
    results <- list()
    for (k in 1:length(train_sizes)){
        current_result <- read.csv(file.path("new_methods", "deconvolution_downsize",
                                             sprintf("methylation_%s_downsize_%d.csv", current_study, train_sizes[k])))
        if ("Success" %in% colnames(current_result)){
            current_result <- current_result %>% filter(Success)
            current_result$Success <- NULL
        }
        results[[k]] <- current_result

    }
    combined_results <- do.call(rbind, results)
    deconv_results[[current_study]] <- combined_results

}


# load existing method performances


competitor_results <- list()

for (j in 1:length(studies)){
    current_study <- studies[j]
    current_result <- read.csv(file.path("existing_methods",
                                         "output",
                                         sprintf("methylation_%s_train_downsize.csv", current_study)),
                               row.names=1)
    competitor_results[[current_study]] <- current_result
}

boxplots_list <- list()

for (j in 1:length(studies)){

    current_study <- studies[j]
    competitor_result <- competitor_results[[current_study]]
    deconv_result <- deconv_results[[current_study]]

    # mae for all test samples
    deconv_result_full <- deconv_result %>% select(TrainSize, MAE) %>%
        rename(SampleSize=TrainSize)
    deconv_result_full$Method <- "Deconv"
    deconv_result_full$Test <- "All Test Samples"
    competitor_result_full <- competitor_result %>% select(SampleSize, EN, RF, LGB)
    competitor_result_full_long <- pivot_longer(competitor_result_full, cols=c("EN", "RF", "LGB"),
                                           names_to="Method",
                                           values_to="MAE")
    competitor_result_full_long$Test <- "All Test Samples"
    combined_results_full <- rbind(deconv_result_full, competitor_result_full_long)


    # mae for test samples whose ages are within the train sample range
    deconv_result_subset <- deconv_result %>% select(TrainSize, MAE_subset) %>%
        rename(SampleSize=TrainSize, MAE=MAE_subset)
    deconv_result_subset$Method <- "Deconv"
    deconv_result_subset$Test <- "Subset Test Samples"
    competitor_result_subset <- competitor_result %>% select(SampleSize, EN_subset, RF_subset, LGB_subset) %>%
        rename(EN=EN_subset, RF=RF_subset, LGB=LGB_subset)
    competitor_result_subset_long <- pivot_longer(competitor_result_subset,
                                                  cols=c("EN", "RF", "LGB"),
                                                names_to="Method",
                                                values_to="MAE")
    competitor_result_subset_long$Test <- "Subset Test Samples"
    combined_results_subset <- rbind(deconv_result_subset, competitor_result_subset_long)


    combined_results <- rbind(combined_results_full, combined_results_subset)
    combined_results$SampleSize <- factor(combined_results$SampleSize, levels=train_sizes)

    min_mae <- round(min(combined_results$MAE))
    max_mae <- round(max(combined_results$MAE))

    boxplots_list[[current_study]] <- ggplot(combined_results, aes(x=SampleSize, y=MAE, color=Method)) +
        geom_boxplot(alpha=0.8, outlier.shape=NA) +
        theme_bw() +
        facet_wrap(~Test)+
        scale_y_continuous(breaks=seq(min_mae, max_mae, 1))+
        xlab("Train Data Size") + ylab("Test MAE") +
        ggtitle(current_study)


}



for (j in 1:length(studies)){

    current_study <- studies[j]
    ofile1 <- sprintf("%s.png", current_study)
    ofile2 <- sprintf("%s.svg", current_study)

    ggsave(file.path("comparison", "boxplots", "png", ofile1),
           plot=boxplots_list[[current_study]], width=8, height=4)

    ggsave(file.path("comparison", "boxplots", "svg", ofile1),
           boxplots_list[[current_study]], width=8, height=4)

}

