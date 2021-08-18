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



#the above is too slow, try something else
setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/integrated/neighbor/")

gc.combined.seurat <- readRDS("../limma.regressed.mnn.integrated.seurat.rds")

gc.scaledata.hvg <- gc.combined.seurat@assays$RNA@scale.data[
                                        VariableFeatures(gc.combined.seurat)[1:1000],
                                        gc.combined.seurat$orig.ident=="gcdata"]
egc.scaledata.hvg <- gc.combined.seurat@assays$RNA@scale.data[
                                        VariableFeatures(gc.combined.seurat)[1:1000],
                                        gc.combined.seurat$orig.ident=="egc.data"]
gc.scaledata.hvg <- t(gc.scaledata.hvg)
egc.scaledata.hvg <- t(egc.scaledata.hvg)

# system.time(cdist(gc.scaledata.hvg[1:2,], egc.scaledata.hvg))

distance.hvg <- cdist(gc.scaledata.hvg, egc.scaledata.hvg)
rownames(distance.hvg) <- rownames(gc.scaledata.hvg)
colnames(distance.hvg) <- rownames(egc.scaledata.hvg)



#save rds
saveRDS(distance.hvg, "distance.hvg.rds")
distance.hvg <- readRDS("distance.hvg.rds")

# system.time(distance.hvg.list <- apply(distance.hvg[1:2,], 1,
#                 function(x) head(colnames(distance.hvg)[order(x)], 10)))

distance.hvg.list <- apply(distance.hvg, 1,
                function(x) head(order(x), 10))
#avoid string matrices
#store the order/index

distance.hvg.list <- t(distance.hvg.list)

distance.hvg.list.celltype <- distance.hvg.list
# gc.combined.seurat@meta.data[distance.hvg.list[1,], "cell.type"]
for (i in 1:nrow(distance.hvg.list)){
  distance.hvg.list.celltype[i,] = gc.combined.seurat@meta.data[
                                                    distance.hvg.list[i,], "cell.type"]
}
