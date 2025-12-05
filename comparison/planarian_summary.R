

rm(list=ls())
library(ggplot2)
library(dplyr)
library(deconvolution)
library(stringr)

# load existing method result
competitors <- read.csv("existing_methods/output/planarian_test_prediction.csv")
sample_names <- competitors$Sample
truth <- str_extract_all(sample_names, "\\d+\\.?\\d*") |> unlist() |> as.integer()

mae_en <- mean(abs(competitors$EN - competitors$Age))
en_plot <- viz_predict(truth=truth, predicted=competitors$EN,
                       title=sprintf("Elastic Net (MAE %.2f)", mae_en),
                       size=7)
mae_rf <- mean(abs(competitors$RF - competitors$Age))
rf_plot <- viz_predict(truth=truth, predicted=competitors$RF,
                       title=sprintf("Random Forest (MAE %.2f)", mae_rf),
                       size=7)
mae_lgb <- mean(abs(competitors$LGB - competitors$Age))
lgb_plot <- viz_predict(truth=truth, predicted=competitors$RF,
                       title=sprintf("LightGBM (MAE %.2f)", mae_lgb),
                       size=7)


# load deconvolution results

deconv_meanresults <- data.frame(Sample=competitors$Sample,
                             Deconv0_original=0,
                             Deconv0_normalized=0,
                             Deconv1_original=0,
                             Deconv1_normalized=0,
                             Deconv2_original=0,
                             Deconv2_normalized=0,
                             Deconv3_original=0,
                             Deconv3_normalized=0)


deconv_moderesults <- data.frame(Sample=competitors$Sample,
                                 Deconv0_original=0,
                                 Deconv0_normalized=0,
                                 Deconv1_original=0,
                                 Deconv1_normalized=0,
                                 Deconv2_original=0,
                                 Deconv2_normalized=0,
                                 Deconv3_original=0,
                                 Deconv3_normalized=0)


for (i in 1:nrow(deconv_meanresults)){

    current_sample <- deconv_meanresults$Sample[i]
    fname <- sprintf("planarian_%s.csv", current_sample)
    current_results <- read.csv(file.path("new_methods", "deconvolution", "planarian",
                                          sprintf("planarian_%s.csv", current_sample)))

    deconv_meanresults[i, seq(2, 9)] <- current_results$pred_mean
    deconv_moderesults[i, seq(2, 9)] <- current_results$pred_mode

}

combined_performance <- competitors %>% inner_join(deconv_meanresults, by="Sample")
write.csv(combined_performance, "comparison/planarian_results.csv",
          row.names=FALSE, quote=FALSE)

mae_mean_original_0 <- mean(abs(deconv_meanresults$Deconv0_original - truth))
deconv_original_0_plot <- viz_predict(truth=truth, predicted=deconv_meanresults$Deconv0_original,
                        title=sprintf("Type 0 (Raw Input)\n(MAE %.2f)", mae_mean_original_0),
                        size=7)


mae_mean_normalized_0 <- mean(abs(deconv_meanresults$Deconv0_normalized - truth))
deconv_normalized_0_plot <- viz_predict(truth=truth, predicted=deconv_meanresults$Deconv0_normalized,
                                      title=sprintf("Type 0 (Standardized Input)\n(MAE %.2f)",
                                                    mae_mean_normalized_0),
                                      size=7)


mae_mean_original_1 <- mean(abs(deconv_meanresults$Deconv1_original - truth))
deconv_original_1_plot <- viz_predict(truth=truth, predicted=deconv_meanresults$Deconv1_original,
                                        title=sprintf("Type 1 (Raw Input)\n(MAE %.2f)",
                                                      mae_mean_original_1),
                                      size=7)


mae_mean_normalized_1 <- mean(abs(deconv_meanresults$Deconv1_normalized - truth))
deconv_normalized_1_plot <- viz_predict(truth=truth, predicted=deconv_meanresults$Deconv1_normalized,
                                      title=sprintf("Type 1 (Standardized Input)\n(MAE %.2f)",
                                                    mae_mean_normalized_1),
                                      size=7)

mae_mean_original_2 <- mean(abs(deconv_meanresults$Deconv2_original - truth))
deconv_original_2_plot <- viz_predict(truth=truth, predicted=deconv_meanresults$Deconv2_original,
                                      title=sprintf("Type 2 (Raw Input)\n(MAE %.2f)",
                                                    mae_mean_original_2),
                                      size=7)

mae_mean_normalized_2 <- mean(abs(deconv_meanresults$Deconv2_normalized - truth))
deconv_normalized_2_plot <- viz_predict(truth=truth, predicted=deconv_meanresults$Deconv2_normalized,
                                      title=sprintf("Type 2 (Standardized Input)\n(MAE %.2f)",
                                                    mae_mean_normalized_2),
                                      size=7)


mae_mean_original_3 <- mean(abs(deconv_meanresults$Deconv3_original - truth))
deconv_original_3_plot <- viz_predict(truth=truth, predicted=deconv_meanresults$Deconv3_original,
                                      title=sprintf("Type 3 (Raw Input)\n(MAE %.2f)",
                                                    mae_mean_original_3),
                                      size=7)

mae_mean_normalized_3 <- mean(abs(deconv_meanresults$Deconv3_normalized - truth))
deconv_normalized_3_plot <- viz_predict(truth=truth, predicted=deconv_meanresults$Deconv3_normalized,
                                      title=sprintf("Type 3 (Standardized Input)\n(MAE %.2f)",
                                                    mae_mean_normalized_3),
                                      size=7)



library(patchwork)

combined_plots <- wrap_plots(deconv_original_0_plot, deconv_original_1_plot, deconv_original_2_plot, deconv_original_3_plot,
                             deconv_normalized_0_plot, deconv_normalized_1_plot, deconv_normalized_2_plot, deconv_normalized_3_plot,
                             en_plot, rf_plot, lgb_plot, ncol=4)

ggsave("comparison/planarian_comparison.png", combined_plots,
       width=11, height=7)
