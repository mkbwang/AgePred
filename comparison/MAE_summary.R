
library(dplyr)
library(ggplot2)
library(patchwork)
library(openxlsx)

rm(list=ls())
# metadata_folder <- "extdata"
# deconv_folder <- "new_methods/deconvolution"
# competitor_folder <- "existing_methods/output"

metadata_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata"
deconv_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/new_methods/deconvolution"
mle_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/new_methods/MLE"
competitor_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/existing_methods/output"


# first look at methylation
studies_info <- read.csv(file.path(metadata_folder, "methylation", "studies_info.csv"))
studies_info <- studies_info %>% filter(Total >= 80 & Max_Age - Min_Age > 30)
studies <- studies_info$Study

methylation_result_df <- data.frame(Study=studies,
                                    SampleSize=studies_info$Total,
                                    MLE=0,
                                    Deconv0_original=0,
                                    Deconv0_normalized=0,
                                    Deconv1_original=0,
                                    Deconv1_normalized=0,
                                    EN=0,
                                    RF=0,
                                    LGB=0)

methylation_plots <- list()

## look at study by study
for (j in 1:length(studies)){

    study <- studies[j]
  
    ## deconvolution without sum to one constraint
    deconv_type0 <- read.xlsx(file.path(deconv_folder,
                                        sprintf("methylation_%s_type0.xlsx", study)),
                              sheet="Prediction",
                              rowNames=TRUE)
    methylation_result_df$Deconv0_original[j] <- mean(abs(deconv_type0$age - deconv_type0$Prediction_mean_original))
    methylation_result_df$Deconv0_normalized[j] <- mean(abs(deconv_type0$age - deconv_type0$Prediction_mean_normalized))
    deconv0_original <- data.frame(Type=sprintf("Deconv0_original (MAE %.2f)", mean(abs(deconv_type0$age - deconv_type0$Prediction_mean_original))),
                                   Age=deconv_type0$age,
                                   Prediction=deconv_type0$Prediction_mean_original)
    deconv0_normalized <- data.frame(Type=sprintf("Deconv0_normalized (MAE %.2f)", mean(abs(deconv_type0$age - deconv_type0$Prediction_mean_normalized))),
                                     Age=deconv_type0$age,
                                     Prediction=deconv_type0$Prediction_mean_normalized)
    
    ## deconvolution without sum to one constraint
    deconv_type1 <- read.xlsx(file.path(deconv_folder,
                                        sprintf("methylation_%s_type1.xlsx", study)),
                              sheet="Prediction",
                              rowNames=TRUE)
    methylation_result_df$Deconv1_original[j] <- mean(abs(deconv_type1$age - deconv_type1$Prediction_mean_original))
    methylation_result_df$Deconv1_normalized[j] <- mean(abs(deconv_type1$age - deconv_type1$Prediction_mean_normalized))
    deconv1_original <- data.frame(Type=sprintf("Deconv1_original (MAE %.2f)", mean(abs(deconv_type1$age - deconv_type1$Prediction_mean_original))),
                                   Age=deconv_type1$age,
                                   Prediction=deconv_type1$Prediction_mean_original)
    deconv1_normalized <- data.frame(Type=sprintf("Deconv1_normalized (MAE %.2f)", mean(abs(deconv_type1$age - deconv_type1$Prediction_mean_normalized))),
                                   Age=deconv_type1$age,
                                   Prediction=deconv_type1$Prediction_mean_normalized)

    # MLE 
    mle <- read.xlsx(file.path(mle_folder, sprintf("methylation_%s.xlsx", study)),
                     sheet="Prediction", rowNames=TRUE)
    methylation_result_df$MLE[j] <- mean(abs(mle$age - mle$mle_pred))
    mle_output <- data.frame(Type=sprintf("MLE (MAE %.2f)", mean(abs(mle$age - mle$mle_pred))),
                             Age=mle$age,
                             Prediction=mle$mle_pred)
    
    
    competitor_result <- read.csv(file.path(competitor_folder,
                                            sprintf("methylation_%s.csv", study)),
                                  row.names=1)

    methylation_result_df$EN[j] <- mean(abs(competitor_result$age - competitor_result$EN))
    EN_output <- data.frame(Type=sprintf("Elastic Net (MAE %.2f)", mean(abs(competitor_result$age - competitor_result$EN))),
                            Age=competitor_result$age,
                            Prediction=competitor_result$EN)

    methylation_result_df$RF[j] <- mean(abs(competitor_result$age - competitor_result$RF))
    RF_output <- data.frame(Type=sprintf("Random Forest (MAE %.2f)", mean(abs(competitor_result$age - competitor_result$RF))),
                            Age=competitor_result$age,
                            Prediction=competitor_result$RF)

    methylation_result_df$LGB[j] <- mean(abs(competitor_result$age - competitor_result$LGB))
    LGB_output <- data.frame(Type=sprintf("LightGBM (MAE %.2f)", mean(abs(competitor_result$age - competitor_result$LGB))),
                            Age=competitor_result$age,
                            Prediction=competitor_result$LGB)


    combined_predictions_long <- rbind(deconv0_original, deconv1_original, mle_output,
                                       EN_output, RF_output, LGB_output)
    combined_predictions_long$Type <- factor(combined_predictions_long$Type,
                                             levels=unique(combined_predictions_long$Type))

    min_age <- studies_info$Min_Age[j] |> round()
    max_age <- studies_info$Max_Age[j] |> round()


    prediction_plot <- ggplot(combined_predictions_long, aes(x=Age, y=Prediction)) +
        geom_point(alpha=0.8) +
        geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
        facet_wrap(~Type, ncol=3) +
        scale_x_continuous(limits=c(min_age, max_age), breaks=seq(min_age, max_age, 5))+
        scale_y_continuous(limits=c(min_age, max_age), breaks=seq(min_age, max_age, 5))+
        xlab("True Age") + ylab("Prediction") + ggtitle(study)
  
    methylation_plots[[j]] <- prediction_plot
    png_name <- sprintf("methylation_%s_logit.png", study)
    svg_name <- sprintf("methylation_%s_logit.svg", study)

    ggsave(file.path("comparison", "scatter", "svg", svg_name),
           plot = prediction_plot, width = 8, height = 5, units = "in")

    ggsave(file.path("comparison", "scatter", "png", png_name),
           plot = prediction_plot, width = 8, height = 5, units = "in")

}


methylation_result_df_sorted <- methylation_result_df %>% arrange(SampleSize)
methylation_result_df_sorted$Deconv0_normalized <- NULL
methylation_result_df_sorted$Deconv1_normalized <- NULL
methylation_result_df_sorted$ID <- seq(1, nrow(methylation_result_df_sorted))
methylation_result_df_sorted$sname <- sprintf("%s (%d)", methylation_result_df_sorted$Study,
                                              methylation_result_df_sorted$SampleSize)

methylation_result_df_sorted$Study <- NULL
methylation_result_df_sorted$SampleSize <- NULL
library(tidyr)
methylation_result_df_sorted_long <- methylation_result_df_sorted %>% pivot_longer(
  cols = c("MLE", "Deconv0_original", "Deconv1_original", "EN", "RF", "LGB"), # Selects the columns to pivot: Q1_Sales, Q2_Sales, Q1_Profit, Q2_Profit
  names_to = "Method", # Names for the new key columns
  values_to = "MAE" # Name for the new value column
)

mae_plot <- ggplot(methylation_result_df_sorted_long, aes(x=ID, y=MAE, color=Method)) +
  geom_point(alpha=0.7) + geom_line(alpha=0.7) + 
  scale_x_continuous(breaks=seq(1, 18), labels=methylation_result_df_sorted$sname) +
  scale_color_manual(values=c( "#A020F0", "#FF00FF", "#00725A", "#B2AC29", "#03396B", "#555555"))+
  xlab("Study") + ylab("MAE") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


write.csv(methylation_result_df, "comparison/methylation_result.csv", row.names=FALSE, quote=FALSE)

# write.csv(methylation_result_df, "/nfs/turbo/sph-ligen/wangmk/AgePred/comparison/methylation_result.csv",
#           row.names=FALSE, quote=FALSE)


# then look at gene expression


gene_meta_files <- list.files(file.path(metadata_folder, "multitissue_transcriptome",
                                        "test"),
                              pattern="metadata*")

tissues <- sapply(gene_meta_files, function(fname){
    strsplit(gsub(".csv", "", fname),
             split="_")[[1]][2]
})


geneexp_result_df <- data.frame(Tissue=tissues,
                                Deconv1_original=0,
                                Deconv1_normalized=0,
                                Deconv2_original=0,
                                Deconv2_normalized=0,
                                Deconv3_original=0,
                                Deconv3_normalized=0,
                                EN=0,
                                RF=0,
                                LGB=0)

geneexp_plots <- list()

## look at tissue by tissue
for (j in 1:length(tissues)){

    tissue <- tissues[j]

    deconv_type1 <- read.xlsx(file.path(deconv_folder,
                                        sprintf("geneexp_%s_type1.xlsx", tissue)),
                              sheet="Prediction",
                              rowNames=TRUE)

    geneexp_result_df$Deconv1_original[j] <- mean(abs(deconv_type1$Age - deconv_type1$Prediction_mean_original))
    geneexp_result_df$Deconv1_normalized[j] <- mean(abs(deconv_type1$Age - deconv_type1$Prediction_mean_normalized))
    deconv1_original <- data.frame(Type=sprintf("Deconv1_original (MAE %.2f)", mean(abs(deconv_type1$Age - deconv_type1$Prediction_mean_original))),
                                   Age=deconv_type1$Age,
                                   Prediction=deconv_type1$Prediction_mean_original)
    deconv1_normalized <- data.frame(Type=sprintf("Deconv1_normalized (MAE %.2f)", mean(abs(deconv_type1$Age - deconv_type1$Prediction_mean_normalized))),
                                     Age=deconv_type1$Age,
                                     Prediction=deconv_type1$Prediction_mean_normalized)


    deconv_type2 <- read.xlsx(file.path(deconv_folder,
                                        sprintf("geneexp_%s_type2.xlsx", tissue)),
                              sheet="Prediction",
                              rowNames=TRUE)
    geneexp_result_df$Deconv2_original[j] <- mean(abs(deconv_type2$Age - deconv_type2$Prediction_mean_original))
    geneexp_result_df$Deconv2_normalized[j] <- mean(abs(deconv_type2$Age - deconv_type2$Prediction_mean_normalized))
    deconv2_original <- data.frame(Type=sprintf("Deconv2_original (MAE %.2f)", mean(abs(deconv_type2$Age - deconv_type2$Prediction_mean_original))),
                                   Age=deconv_type2$Age,
                                   Prediction=deconv_type2$Prediction_mean_original)
    deconv2_normalized <- data.frame(Type=sprintf("Deconv2_normalized (MAE %.2f)", mean(abs(deconv_type2$Age - deconv_type2$Prediction_mean_normalized))),
                                     Age=deconv_type2$Age,
                                     Prediction=deconv_type2$Prediction_mean_normalized)


    deconv_type3 <- read.xlsx(file.path(deconv_folder,
                                        sprintf("geneexp_%s_type3.xlsx", tissue)),
                              sheet="Prediction",
                              rowNames=TRUE)

    geneexp_result_df$Deconv3_original[j] <- mean(abs(deconv_type3$Age - deconv_type3$Prediction_mean_original))
    geneexp_result_df$Deconv3_normalized[j] <- mean(abs(deconv_type3$Age - deconv_type3$Prediction_mean_normalized))
    deconv3_original <- data.frame(Type=sprintf("Deconv3_original (MAE %.2f)", mean(abs(deconv_type3$Age - deconv_type3$Prediction_mean_original))),
                                   Age=deconv_type3$Age,
                                   Prediction=deconv_type3$Prediction_mean_original)

    deconv3_normalized <- data.frame(Type=sprintf("Deconv3_normalized (MAE %.2f)", mean(abs(deconv_type3$Age - deconv_type3$Prediction_mean_normalized))),
                                     Age=deconv_type3$Age,
                                     Prediction=deconv_type3$Prediction_mean_normalized)


    competitor_result <- read.csv(file.path(competitor_folder,
                                            sprintf("geneexp_%s.csv", tissue)),
                                  row.names=1)

    geneexp_result_df$EN[j] <- mean(abs(competitor_result$Age - competitor_result$EN))
    EN_output <- data.frame(Type=sprintf("Elastic Net (MAE %.2f)", mean(abs(competitor_result$Age - competitor_result$EN))),
                            Age=competitor_result$Age,
                            Prediction=competitor_result$EN)

    geneexp_result_df$RF[j] <- mean(abs(competitor_result$Age - competitor_result$RF))
    RF_output <- data.frame(Type=sprintf("Random Forest (MAE %.2f)", mean(abs(competitor_result$Age - competitor_result$RF))),
                            Age=competitor_result$Age,
                            Prediction=competitor_result$RF)

    geneexp_result_df$LGB[j] <- mean(abs(competitor_result$Age - competitor_result$LGB))
    LGB_output <- data.frame(Type=sprintf("LightGBM (MAE %.2f)", mean(abs(competitor_result$Age - competitor_result$LGB))),
                             Age=competitor_result$Age,
                             Prediction=competitor_result$LGB)


    combined_predictions_long <- rbind(deconv1_original, deconv1_normalized,
                                       deconv2_original, deconv2_normalized,
                                       deconv3_original, deconv3_normalized,
                                       EN_output, RF_output, LGB_output)
    combined_predictions_long$Type <- factor(combined_predictions_long$Type,
                                             levels=unique(combined_predictions_long$Type))

    min_Age <- min(deconv_type1$Age) |> round()
    max_Age <- max(deconv_type1$Age) |> round()


    prediction_plot <- ggplot(combined_predictions_long, aes(x=Age, y=Prediction)) +
        geom_point(alpha=0.8) +
        geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
        facet_wrap(~Type, ncol=3) +
        scale_x_continuous(limits=c(min_Age, max_Age), breaks=seq(min_Age, max_Age, 5))+
        scale_y_continuous(breaks=seq(min_Age, max_Age, 5))+
        xlab("True Age") + ylab("Prediction") + ggtitle(tissue)

    png_name <- sprintf("geneexp_%s.png", tissue)
    svg_name <- sprintf("geneexp_%s.svg", tissue)

    ggsave(file.path("comparison", "scatter", "svg", svg_name),
           plot = prediction_plot, width = 8, height = 7, units = "in")

    ggsave(file.path("comparison", "scatter", "png", png_name),
           plot = prediction_plot, width = 8, height = 7, units = "in")

}



write.csv(geneexp_result_df, "comparison/geneexp_result.csv", row.names=FALSE, quote=FALSE)

# write.csv(geneexp_result_df, "/nfs/turbo/sph-ligen/wangmk/AgePred/comparison/geneexp_result.csv",
#           row.names=FALSE, quote=FALSE)


