
rm(list=ls())
library(MuSiC)
library(SingleCellExperiment)

GSE50244.bulk.eset = readRDS('extdata/singlecell/GSE50244bulkeset.rds')
# GSE50244.bulk.eset
bulk.mtx <- exprs(GSE50244.bulk.eset)
bulk.meta <- pData(GSE50244.bulk.eset)


# sce template
EMTAB.sce <- readRDS("extdata/singlecell/EMTABsce_healthy.rds")
EMTAB_meta <- colData(EMTAB.sce) |> as.data.frame()
EMTAB_meta[EMTAB_meta$sampleID==1,"cellType"] |> table()

# run music for deconvolution
est.prop.GSE50244 <- music_prop(bulk.mtx=bulk.mtx,
                                sc.sce=EMTAB.sce,
                                clusters="cellType",
                                samples="sampleID",
                                select.ct=c("alpha", "beta", "delta",
                                            "gamma", "acinar", "ductal")) 

View(est.prop.GSE50244$Est.prop.weighted)
View(est.prop.GSE50244$Est.prop.allgene)



