

library(openxlsx)

metadata_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata"
deconv_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/new_methods/deconvolution"
competitor_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/existing_methods/output"

# first look at methylation
studies_info <- read.csv(file.path(metadata_folder, "methylation", "studies_info.csv"))
studies_info <- studies_info %>% filter(Total >= 80 & Max_Age - Min_Age > 30)
studies <- studies_info$Study


for (j in 1:length(studies)){
    print(j)
    study <- studies[j]

    prediction_df_1 <- read.xlsx(file.path(deconv_folder,
                                         sprintf("methylation_%s_type1.xlsx", study)),
                               sheet="Prediction", rowNames=TRUE)

    weights_original_1 <- read.xlsx(file.path(deconv_folder,
                                            sprintf("methylation_%s_type1.xlsx", study)),
                                  sheet="Weights_original", rowNames=TRUE) |> as.matrix()

    unique_ages <- sapply(colnames(weights_original_1),
                          function(cname){
                              strsplit(cname, split="_")[[1]][2]
                          }) |> as.integer()

    weights_normalized_1 <- read.xlsx(file.path(deconv_folder,
                                              sprintf("methylation_%s_type1.xlsx", study)),
                                    sheet="Weights_normalized", rowNames=TRUE) |> as.matrix()
    prediction_df_2 <- read.xlsx(file.path(deconv_folder,
                                           sprintf("methylation_%s_type2.xlsx", study)),
                                 sheet="Prediction", rowNames=TRUE) |> as.matrix()
    weights_original_2 <- read.xlsx(file.path(deconv_folder,
                                              sprintf("methylation_%s_type2.xlsx", study)),
                                    sheet="Weights_original", rowNames=TRUE) |> as.matrix()
    weights_normalized_2 <- read.xlsx(file.path(deconv_folder,
                                                sprintf("methylation_%s_type2.xlsx", study)),
                                      sheet="Weights_normalized", rowNames=TRUE) |> as.matrix()
    prediction_df_3 <- read.xlsx(file.path(deconv_folder,
                                           sprintf("methylation_%s_type3.xlsx", study)),
                                 sheet="Prediction", rowNames=TRUE)
    weights_original_3 <- read.xlsx(file.path(deconv_folder,
                                              sprintf("methylation_%s_type3.xlsx", study)),
                                    sheet="Weights_original", rowNames=TRUE) |> as.matrix()
    weights_normalized_3 <- read.xlsx(file.path(deconv_folder,
                                                sprintf("methylation_%s_type3.xlsx", study)),
                                      sheet="Weights_normalized", rowNames=TRUE) |> as.matrix()


    weight_plot_list <- list()

    for (k in 1:nrow(prediction_df_1)){

        sample_name <- rownames(prediction_df_1)[k]
        truth <- prediction_df_1$age[k]
        weight_original_df <- data.frame(Age=rep(unique_ages, 3),
                                         Weight=c(weights_original_1[k, ],
                                                  weights_original_2[k, ],
                                                  weights_original_3[k, ]))
        weight_original_df$Type <- rep(c("Type1", "Type2", "Type3"),
                                       each=length(unique_ages))
        weight_original_df$Feature <- "Original Feature Value"

        weight_normalized_df <- data.frame(Age=rep(unique_ages, 3),
                                         Weight=c(weights_normalized_1[k, ],
                                                  weights_normalized_2[k, ],
                                                  weights_normalized_3[k, ]))
        weight_normalized_df$Type <- rep(c("Type1", "Type2", "Type3"),
                                       each=length(unique_ages))
        weight_normalized_df$Feature <- "Standardized Feature Value"

        combined_weights_df <- rbind(weight_original_df, weight_normalized_df)
        combined_weights_df$Type <- factor(combined_weights_df$Type,
                                           levels=c("Type1", "Type2", "Type3"))
        combined_weights_df$Feature <- factor(combined_weights_df$Feature,
                                              levels=c("Original Feature Value",
                                                       "Standardized Feature Value"))

        max_weight <- max(combined_weights_df$Weight)
        min_label <- min(combined_weights_df$Age)
        max_label <- max(combined_weights_df$Age)

        weight_plot_list[[k]] <- ggplot(combined_weights_df, aes(x=Age, y=Weight, color=Type)) + geom_point(alpha=0.6) +
            geom_line(alpha=0.6) + xlab("Age") + ylab("Weight") +
            scale_x_continuous(limits=c(min_label, max_label),
                               breaks=seq(min_label, max_label, 5)) +
            scale_y_continuous(limits=c(0, max_weight),
                               breaks=seq(0, max_weight, 0.05))+
            geom_vline(xintercept=truth, linetype="dashed")+
            scale_color_manual(values=c("#A93226", "#196F3D", "#154360"))+
            facet_wrap(~Feature, ncol=2) + xlab("Age") + ylab("Weight")+
            ggtitle(sample_name)

    }

    ofile <- sprintf("Methylation_%s.pdf", study)
    pdf(file.path("comparison", "weights", ofile), width = 7, height = 4)
    for(i in 1:length(weight_plot_list)){
        print(weight_plot_list[[i]])
    }
    dev.off()

}

write.csv(studies_info, file.path(metadata_folder, "methylation", "studies_info_subset.csv"),
          row.names=FALSE)

# look at gene expression

gene_meta_files <- list.files(file.path(metadata_folder, "multitissue_transcriptome",
                                        "test"),
                              pattern="metadata*")

tissues <- sapply(gene_meta_files, function(fname){
    strsplit(gsub(".csv", "", fname),
             split="_")[[1]][2]
})

tissues_info <- data.frame(Tissue=tissues,
                           Total=0,
                           Train=0,
                           Test=0,
                           Min_Age=0,
                           Max_Age=0)

## check metadata information
for (j in 1:length(tissues)){

    tissue <- tissues[j]
    train_metadata <- read.csv(file.path(metadata_folder, "multitissue_transcriptome", "train",
                                         sprintf("metadata_%s.csv", tissue)))

    test_metadata <- read.csv(file.path(metadata_folder, "multitissue_transcriptome", "test",
                                         sprintf("metadata_%s.csv", tissue)))
    total_metadata <- rbind(train_metadata, test_metadata)

    tissues_info$Total[j] <- nrow(train_metadata) + nrow(test_metadata)
    tissues_info$Train[j] <- nrow(train_metadata)
    tissues_info$Test[j] <- nrow(test_metadata)
    tissues_info$Min_Age[j] <- min(total_metadata$Age)
    tissues_info$Max_Age[j] <- max(total_metadata$Age)

}
write.csv(tissues_info, file.path(metadata_folder, "multitissue_transcriptome", "tissues_info.csv"),
          row.names=FALSE)

for (j in 1:length(tissues)){
    print(j)
    tissue <- tissues[j]

    prediction_df_1 <- read.xlsx(file.path(deconv_folder,
                                           sprintf("geneexp_%s_type1.xlsx", tissue)),
                                 sheet="Prediction", rowNames=TRUE)

    weights_original_1 <- read.xlsx(file.path(deconv_folder,
                                              sprintf("geneexp_%s_type1.xlsx", tissue)),
                                    sheet="Weights_original", rowNames=TRUE) |> as.matrix()

    unique_ages <- sapply(colnames(weights_original_1),
                          function(cname){
                              strsplit(cname, split="_")[[1]][2]
                          }) |> as.integer()

    weights_normalized_1 <- read.xlsx(file.path(deconv_folder,
                                                sprintf("geneexp_%s_type1.xlsx", tissue)),
                                      sheet="Weights_normalized", rowNames=TRUE) |> as.matrix()
    prediction_df_2 <- read.xlsx(file.path(deconv_folder,
                                           sprintf("geneexp_%s_type2.xlsx", tissue)),
                                 sheet="Prediction", rowNames=TRUE) |> as.matrix()
    weights_original_2 <- read.xlsx(file.path(deconv_folder,
                                              sprintf("geneexp_%s_type2.xlsx", tissue)),
                                    sheet="Weights_original", rowNames=TRUE) |> as.matrix()
    weights_normalized_2 <- read.xlsx(file.path(deconv_folder,
                                                sprintf("geneexp_%s_type2.xlsx", tissue)),
                                      sheet="Weights_normalized", rowNames=TRUE) |> as.matrix()
    prediction_df_3 <- read.xlsx(file.path(deconv_folder,
                                           sprintf("geneexp_%s_type3.xlsx", tissue)),
                                 sheet="Prediction", rowNames=TRUE)
    weights_original_3 <- read.xlsx(file.path(deconv_folder,
                                              sprintf("geneexp_%s_type3.xlsx", tissue)),
                                    sheet="Weights_original", rowNames=TRUE) |> as.matrix()
    weights_normalized_3 <- read.xlsx(file.path(deconv_folder,
                                                sprintf("geneexp_%s_type3.xlsx", tissue)),
                                      sheet="Weights_normalized", rowNames=TRUE) |> as.matrix()


    weight_plot_list <- list()

    for (k in 1:nrow(prediction_df_1)){

        sample_name <- rownames(prediction_df_1)[k]
        truth <- prediction_df_1$age[k]
        weight_original_df <- data.frame(Age=rep(unique_ages, 3),
                                         Weight=c(weights_original_1[k, ],
                                                  weights_original_2[k, ],
                                                  weights_original_3[k, ]))
        weight_original_df$Type <- rep(c("Type1", "Type2", "Type3"),
                                       each=length(unique_ages))
        weight_original_df$Feature <- "Original Feature Value"

        weight_normalized_df <- data.frame(Age=rep(unique_ages, 3),
                                           Weight=c(weights_normalized_1[k, ],
                                                    weights_normalized_2[k, ],
                                                    weights_normalized_3[k, ]))
        weight_normalized_df$Type <- rep(c("Type1", "Type2", "Type3"),
                                         each=length(unique_ages))
        weight_normalized_df$Feature <- "Standardized Feature Value"

        combined_weights_df <- rbind(weight_original_df, weight_normalized_df)
        combined_weights_df$Type <- factor(combined_weights_df$Type,
                                           levels=c("Type1", "Type2", "Type3"))
        combined_weights_df$Feature <- factor(combined_weights_df$Feature,
                                              levels=c("Original Feature Value",
                                                       "Standardized Feature Value"))

        max_weight <- max(combined_weights_df$Weight)
        min_label <- min(combined_weights_df$Age)
        max_label <- max(combined_weights_df$Age)

        weight_plot_list[[k]] <- ggplot(combined_weights_df, aes(x=Age, y=Weight, color=Type)) + geom_point(alpha=0.6) +
            geom_line(alpha=0.6) + xlab("Age") + ylab("Weight") +
            scale_x_continuous(limits=c(min_label, max_label),
                               breaks=seq(min_label, max_label, 5)) +
            scale_y_continuous(limits=c(0, max_weight),
                               breaks=seq(0, max_weight, 0.05))+
            geom_vline(xintercept=truth, linetype="dashed")+
            scale_color_manual(values=c("#A93226", "#196F3D", "#154360"))+
            facet_wrap(~Feature, ncol=2) + xlab("Age") + ylab("Weight")+
            ggtitle(sample_name)

    }

    ofile <- sprintf("geneexp_%s.pdf", tissue)
    pdf(file.path("comparison", "weights", ofile), width = 7, height = 4)
    for(i in 1:length(weight_plot_list)){
        print(weight_plot_list[[i]])
    }
    dev.off()

}



