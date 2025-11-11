
library(arrow)

old_train_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/methylation/train_old"
old_test_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/methylation/test_old"

combined_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/methylation/combined"
new_train_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/methylation/train"
new_test_folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/methylation/test"

files <- list.files(path=old_train_folder, pattern="*.feather")


study_names <- sapply(files, function(longname){
  strsplit(longname, split="[.]")[[1]][1]
})


study_df <- data.frame(Study=study_names, Tissue="", 
                       Total=0, Train=0, Test=0, Min_Age=0, Max_Age=0)
rownames(study_df) <- NULL

set.seed(2025)
for (j in 1:length(study_names)){
  
  print(study_names[j])  
  ## load old train/test split and combine
  biomarker_train_old <- read_feather(file.path(old_train_folder, 
                                                sprintf("%s.feather", study_names[j])))
  
  metadata_train_old <- read.csv(file.path(old_train_folder, 
                                           sprintf("%s.csv", study_names[j])), row.names=1)
  
  biomarker_test_old <- read_feather(file.path(old_test_folder, 
                                               sprintf("%s.feather", study_names[j]))) 
  
  metadata_test_old <- read.csv(file.path(old_test_folder, 
                                           sprintf("%s.csv", study_names[j])), row.names=1)
  
  metadata_combined <- rbind(metadata_train_old, metadata_test_old)
  
  age_filter <- metadata_combined$age > 20 & metadata_combined$age < 90
  metadata_combined <- metadata_combined[age_filter, ]
  write.csv(metadata_combined, 
            file.path(combined_folder, sprintf("%s.csv", study_names[j])),
            quote=FALSE)
  biomarker_combined <- rbind(biomarker_train_old, biomarker_test_old)
  biomarker_combined <- biomarker_combined[age_filter, ]
  
  write_feather(biomarker_combined, 
                file.path(combined_folder, sprintf("%s.feather", study_names[j])))
  
  study_df$Tissue[j] <- paste(unique(metadata_combined$tissue_type), collapse="|")
  study_df$Total[j] <- nrow(metadata_combined)
  study_df$Min_Age[j] <- min(metadata_combined$age)
  study_df$Max_Age[j] <- max(metadata_combined$age)
  
  # resplit into new train/test samples
  min_max_indices <- c(which.min(metadata_combined$age), which.max(metadata_combined$age))
  
  train_indices <- sample(seq(1, nrow(metadata_combined)), 0.7*nrow(metadata_combined))
  test_indices <- setdiff(seq(1, nrow(metadata_combined)), train_indices)
  train_indices <- unique(c(train_indices, min_max_indices))
  test_indices <- setdiff(test_indices, min_max_indices)
  
  study_df$Train[j] <- length(train_indices)
  study_df$Test[j] <- length(test_indices)
  
  biomarker_train <- biomarker_combined[train_indices, ]
  write_feather(biomarker_train, 
                file.path(new_train_folder, sprintf("%s.feather", study_names[j])))
  biomarker_test <- biomarker_combined[test_indices, ]
  write_feather(biomarker_test, 
                file.path(new_test_folder, sprintf("%s.feather", study_names[j])))
  
  metadata_train <- metadata_combined[train_indices, ]
  write.csv(metadata_train, 
            file.path(new_train_folder, sprintf("%s.csv", study_names[j])),
            quote=FALSE)
  metadata_test <- metadata_combined[test_indices, ]
  write.csv(metadata_test, 
            file.path(new_test_folder, sprintf("%s.csv", study_names[j])),
            quote=FALSE)
      
}


write.csv(study_df, "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/methylation/studies_info.csv",
          row.names=FALSE, quote=FALSE)


