
library(DESeq2)
library(dplyr)
library(arrow)
rm(list=ls())

folder <- "extdata/multitissue_transcriptome"

# folder <- "/nfs/turbo/sph-ligen/wangmk/AgePred/extdata/multitissue_transcriptome"


metadata <- read.table(file.path(folder, "meta_filtered.txt"), sep='\t',
                       quote="")
rownames(metadata) <- metadata$SRR.ID

raw_counts <- read.table(file.path(folder, "raw_filtered.txt"), sep='\t',
                         quote="", header=TRUE, row.names=1)

raw_counts_mat <- as.matrix(raw_counts)
# normalize and log transform
dds <- DESeqDataSetFromMatrix(countData=raw_counts_mat, colData=metadata, design = ~1)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind=TRUE)

normalized_counts <- assay(vsd)

# batch correction with ComBat
library(sva)

bc_counts <- ComBat(dat=normalized_counts,
                    batch=as.factor(metadata$Batch))


write_feather(as.data.frame(bc_counts),
              file.path(folder, "combined", "all_samples.feather"))
write.csv(metadata,
          file.path(folder, "combined", "metadata.csv"), quote=FALSE)



# retain healthy samples
healthy_filter <- metadata$Healthy
metadata_healthy <- metadata %>% filter(Healthy)

counts_healthy <- bc_counts[, healthy_filter]


metadata_heart <- metadata_healthy %>% filter(grepl("Heart", Tissue))
metadata_blood <- metadata_healthy %>% filter(Tissue == "Blood;PBMC")
metadata_retina <- metadata_healthy %>% filter(Tissue == "Retina")
metadata_adipose <- metadata_healthy %>% filter(Tissue == "Adipose;")


metadata_list <- list(metadata_heart, metadata_blood, metadata_retina, metadata_adipose)


set.seed(2025)
for (j in 1:length(metadata_list)){
  print(j)

  metadata <- metadata_list[[j]]
  tissue_name <- strsplit(unique(metadata$Tissue), split=";")[[1]][1]
  counts_combined <- counts_healthy[, metadata$SRR.ID]
  write.csv(as.data.frame(counts_combined),
            file.path(folder, "combined", sprintf("%s_counts.csv", tissue_name)),
            quote=FALSE)
  write.csv(metadata, file.path(folder, "combined", sprintf("metadata_%s.csv", tissue_name)),
            quote=FALSE)

  # split into train and test
  min_max_ids <- c(which.min(metadata$Age), which.max(metadata$Age))
  train_ids <- sample(1:nrow(metadata), nrow(metadata)*0.7)
  test_ids <- setdiff(1:nrow(metadata), train_ids)

  train_ids <- unique(c(train_ids, min_max_ids))
  test_ids <- setdiff(test_ids, min_max_ids)

  metadata_train <- metadata[train_ids, ]
  counts_train <- counts_combined[, train_ids]


  write.csv(as.data.frame(counts_train),
            file.path(folder, "train", sprintf("%s_counts.csv", tissue_name)),
            quote=FALSE)
  write.csv(metadata_train, file.path(folder, "train", sprintf("metadata_%s.csv", tissue_name)),
            quote=FALSE)

  metadata_test <- metadata[test_ids, ]
  counts_test <- counts_combined[, test_ids]

  write.csv(as.data.frame(counts_test),
            file.path(folder, "test", sprintf("%s_counts.csv", tissue_name)),
            quote=FALSE)
  write.csv(metadata_test, file.path(folder, "test", sprintf("metadata_%s.csv", tissue_name)),
            quote=FALSE)

}

