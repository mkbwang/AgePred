
rm(list=ls())
library(dplyr)
library(DWLS)



#load single-cell data and labels
load("extdata/singlecell/dataSC.RData")
load("extdata/singlecell/dataBulk.RData")
load("extdata/singlecell/trueLabels.RData")

labels<-trueLabels

#change to real labels
newcat<-c("NonCycISC","CycISC","TA","Ent","PreEnt","Goblet","Paneth","Tuft","EE")
for (i in 1:length(newcat)){
  labels[which(labels==(i-1))]<-newcat[i]
}

labels_df <- data.frame(Label=labels, ID=seq(1, length(labels)))
labels_df <- labels_df %>% group_by(Label) %>% slice_sample(prop=0.1)

labels_subset <- labels_df$Label
dataSC_subset <- dataSC[, labels_df$ID]
prevalence <- rowMeans(dataSC_subset > 0)
dataSC_subset <- dataSC_subset[prevalence > 0.5, ]

Signature <- buildSignatureMatrixMAST(scdata=dataSC_subset,
                                      id=labels_subset,
                                      path="competitors/DWLS",
                                      diff.cutoff=0.5,pval.cutoff=0.01)

tr <- trimData(Signature, dataBulk)
View(dataSC_subset[rownames(Signature), ])

solDWLS <- solveDampenedWLS(tr$sig, tr$bulk)
solSVR <- solveSVR(tr$sig, tr$bulk)



