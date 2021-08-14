library(tidyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(limma)
library(batchelor)
library(scran)
library(scater)
library(bluster)
library(rdist)
# library(harmony)
# library(Biocmanager)
# library(SingleCellExperiment)

setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/integrated/neighbor/")

gc.combined.seurat <- readRDS("../limma.regressed.mnn.integrated.seurat.rds")

gc.scaledata <- gc.combined.seurat@assays$RNA@scale.data[,
                                        gc.combined.seurat$orig.ident=="gcdata"]
egc.scaledata <- gc.combined.seurat@assays$RNA@scale.data[,
                                        gc.combined.seurat$orig.ident=="egc.data"]
gc.scaledata <- t(gc.scaledata)
egc.scaledata <- t(egc.scaledata)

distance <- cdist(gc.scaledata, egc.scaledata)
write.csv(distance, "distance.csv")

# #test
# A <- matrix(c(0,0,2,2), byrow=T, nrow=2)
# B <- matrix(c(4,4), nrow=1)
# cdist(A,B)
