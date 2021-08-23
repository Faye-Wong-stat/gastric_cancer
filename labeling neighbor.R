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
saveRDS(distance, "distance.rds")

# #test
# A <- matrix(c(0,0,2,2), byrow=T, nrow=2)
# B <- matrix(c(4,4), nrow=1)
# cdist(A,B)



#the above is too slow, try something else
setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/integrated/neighbor/")

gc.combined.seurat <- readRDS("../limma.regressed.mnn.integrated.seurat.rds")

table(gc.combined.seurat[,gc.combined.seurat$orig.ident=="egc.data"]@meta.data$cell.type)
# B cell     Cancer cell      Chief cell              EC      Enterocyte
#   2158             774              44             994            3364
# Enteroendocrine      Fibroblast             GMC     Goblet cell      Macrophage
#   2814            1381            1886             552             397
# Mast cell             MSC  Neck-like cell              PC             PMC
#    227            1302             383            1133           13026
# SM cell          T cell
#    279            1618


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

#try 100 neighbors
# system.time(distance.hvg.list <- apply(distance.hvg[1:2,], 1,
#                 function(x) head(order(x), 100)))

distance.hvg.list <- apply(distance.hvg, 1,
                function(x) head(order(x), 100))
#avoid string matrices
#store the order/index
distance.hvg.list <- t(distance.hvg.list)
# rownames(distance.hvg.list) <- rownames(gc.scaledata.hvg)

saveRDS(distance.hvg.list, "distance.hvg.list.100.rds")
distance.hvg.list <- readRDS("distance.hvg.list.100.rds")

distance.hvg.list.celltype <- as.data.frame(matrix(NA, nrow=nrow(distance.hvg.list), ncol=1))
rownames(distance.hvg.list.celltype) <- rownames(distance.hvg.list)
# rownames(distance.hvg.list.celltype) <- rownames(gc.scaledata.hvg)
colnames(distance.hvg.list.celltype) <- "80.100"

gc.cell.number = nrow(distance.hvg.list.celltype)
for (i in 1:nrow(distance.hvg.list.celltype)){
  neighbor.cell.types = table(
          gc.combined.seurat$cell.type[(distance.hvg.list[i,]+gc.cell.number)])
  if (identical(names(neighbor.cell.types)[neighbor.cell.types>=80], character(0))){
    distance.hvg.list.celltype[i,] = NA
  } else {
    distance.hvg.list.celltype[i,] = names(neighbor.cell.types)[neighbor.cell.types>=80]
  }
  # print(i)
}

gc.combined.seurat$neighbor_cell.type_hvg100 <- NA
gc.combined.seurat$neighbor_cell.type_hvg100[1:nrow(distance.hvg.list.celltype)] <-
                        distance.hvg.list.celltype[,1]

# saveRDS(gc.combined.seurat, "integrated.seurat.neighbor100.rds")

pdf("umap query neighbor hvg 80.100.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="neighbor_cell.type_hvg100", label=T, repel=T)
dev.off()

#try 50 neighbors
distance.hvg.list.celltype["40.50"] <- NA

for (i in 1:nrow(distance.hvg.list.celltype)){
  neighbor.cell.types = table(
          gc.combined.seurat$cell.type[(distance.hvg.list[i,1:50]+gc.cell.number)])
  if (identical(names(neighbor.cell.types)[neighbor.cell.types>=40], character(0))){
    distance.hvg.list.celltype[i,2] = NA
  } else {
    distance.hvg.list.celltype[i,2] = names(neighbor.cell.types)[neighbor.cell.types>=40]
  }
  # print(i)
}

gc.combined.seurat$neighbor_cell.type_hvg50 <- NA
gc.combined.seurat$neighbor_cell.type_hvg50[1:nrow(distance.hvg.list.celltype)] <-
                        distance.hvg.list.celltype[,2]

pdf("umap query neighbor hvg 40.50.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="neighbor_cell.type_hvg50", label=T, repel=T)
dev.off()

#try 20 neighbors
distance.hvg.list.celltype["16.20"] <- NA

for (i in 1:nrow(distance.hvg.list.celltype)){
  neighbor.cell.types = table(
          gc.combined.seurat$cell.type[(distance.hvg.list[i,1:20]+gc.cell.number)])
  if (identical(names(neighbor.cell.types)[neighbor.cell.types>=16], character(0))){
    distance.hvg.list.celltype[i,3] = NA
  } else {
    distance.hvg.list.celltype[i,3] = names(neighbor.cell.types)[neighbor.cell.types>=16]
  }
  # print(i)
}

gc.combined.seurat$neighbor_cell.type_hvg20 <- NA
gc.combined.seurat$neighbor_cell.type_hvg20[1:nrow(distance.hvg.list.celltype)] <-
                        distance.hvg.list.celltype[,3]

pdf("umap query neighbor hvg 16.20.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="neighbor_cell.type_hvg20", label=T, repel=T)
dev.off()

#try 10 neighbors
distance.hvg.list.celltype["8.10"] <- NA

for (i in 1:nrow(distance.hvg.list.celltype)){
  neighbor.cell.types = table(
          gc.combined.seurat$cell.type[(distance.hvg.list[i,1:10]+gc.cell.number)])
  if (identical(names(neighbor.cell.types)[neighbor.cell.types>=8], character(0))){
    distance.hvg.list.celltype[i,4] = NA
  } else {
    distance.hvg.list.celltype[i,4] = names(neighbor.cell.types)[neighbor.cell.types>=8]
  }
  # print(i)
}

gc.combined.seurat$neighbor_cell.type_hvg10 <- NA
gc.combined.seurat$neighbor_cell.type_hvg10[1:nrow(distance.hvg.list.celltype)] <-
                        distance.hvg.list.celltype[,4]

pdf("umap query neighbor hvg 8.10.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="neighbor_cell.type_hvg10", label=T, repel=T)
dev.off()

#try 5 neighbors
distance.hvg.list.celltype["4.5"] <- NA

for (i in 1:nrow(distance.hvg.list.celltype)){
  neighbor.cell.types = table(
          gc.combined.seurat$cell.type[(distance.hvg.list[i,1:5]+gc.cell.number)])
  if (identical(names(neighbor.cell.types)[neighbor.cell.types>=4], character(0))){
    distance.hvg.list.celltype[i,5] = NA
  } else {
    distance.hvg.list.celltype[i,5] = names(neighbor.cell.types)[neighbor.cell.types>=4]
  }
  # print(i)
}

gc.combined.seurat$neighbor_cell.type_hvg5 <- NA
gc.combined.seurat$neighbor_cell.type_hvg5[1:nrow(distance.hvg.list.celltype)] <-
                        distance.hvg.list.celltype[,5]

pdf("umap query neighbor hvg 4.5.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="neighbor_cell.type_hvg5", label=T, repel=T)
dev.off()



#consistent labels
distance.consistent <- distance.hvg.list.celltype[complete.cases(distance.hvg.list.celltype),]
distance.consistent <-
    distance.consistent[apply(distance.consistent, 1,function(x){length(unique(x))==1}),]
consistent.labels <- rownames(distance.consistent)
length(consistent.labels)
#19547
pdf("umap query neighbor hvg consistent.pdf")
DimPlot(gc.combined.seurat[,consistent.labels], reduction="umap",
        group.by="neighbor_cell.type_hvg5", label=T, repel=T)
dev.off()



#cells get the same labels by two methods
two.methods <- gc.combined.seurat@meta.data[c("cluster_cell.type", "neighbor_cell.type_hvg5")]
two.methods <- two.methods[complete.cases(two.methods),]
two.methods <- two.methods[apply(two.methods, 1,function(x){length(unique(x))==1}),]

cells.same.label.two.methods <- intersect(intersect(consistent.labels,
                                                    consistent.labels.cluster),
                                  rownames(two.methods))

length(consistent.labels.cluster)
#13823
length(consistent.labels)
#19547
length(consistent.labels.cluster)
#13339
nrow(two.methods)
#16242

pdf("umap query neighbor hvg consistent cluster.neighbor.pdf")
DimPlot(gc.combined.seurat[,cells.same.label.two.methods], reduction="umap",
        group.by="neighbor_cell.type_hvg5", label=T, repel=T)
dev.off()





#only look at 10 neighbors and classify using 8 of them
setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/integrated/neighbor/")

gc.combined.seurat <- readRDS("../limma.regressed.mnn.integrated.seurat.rds")
distance.hvg.list <- readRDS("distance.hvg.list.100.rds")

distance.hvg.list.celltype <- as.data.frame(matrix(NA, nrow=nrow(distance.hvg.list), ncol=1))
rownames(distance.hvg.list.celltype) <- rownames(distance.hvg.list)
colnames(distance.hvg.list.celltype) <- "8.10"

gc.cell.number = nrow(distance.hvg.list.celltype)
for (i in 1:nrow(distance.hvg.list.celltype)){
  neighbor.cell.types = table(
          gc.combined.seurat$cell.type[(distance.hvg.list[i,1:10]+gc.cell.number)])
  if (identical(names(neighbor.cell.types)[neighbor.cell.types>=8], character(0))){
    distance.hvg.list.celltype[i,] = NA
  } else {
    distance.hvg.list.celltype[i,] = names(neighbor.cell.types)[neighbor.cell.types>=8]
  }
  # print(i)
}

gc.combined.seurat$neighbor_cell.type_hvg10 <- NA
gc.combined.seurat$neighbor_cell.type_hvg10[1:nrow(distance.hvg.list.celltype)] <-
                        distance.hvg.list.celltype[,1]
pdf("umap query neighbor hvg 8.10.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="neighbor_cell.type_hvg10", label=T, repel=T)
dev.off()



#proportion plot bar plot goes 100%

#h.plori
#careful k-fold cross validation, leave one individual out, do it for each cell type.
#accuracy per cell type, persistency of infection
#if accuracy high, fastmnn, corrected our data into reference data space, predict status in our data, show the distribution of predicted infection per individual.











# distance.hvg.list.celltype$six <- NA
#
# for (i in 1:nrow(distance.hvg.list.celltype)){
#   neighbor.cell.types = table(
#           gc.combined.seurat$cell.type[(distance.hvg.list[i,]+gc.cell.number)])
#   if (identical(names(neighbor.cell.types)[neighbor.cell.types>=6], character(0))){
#     distance.hvg.list.celltype[i,2] = NA
#   } else {
#     distance.hvg.list.celltype[i,2] = names(neighbor.cell.types)[neighbor.cell.types>=6]
#   }
#   # print(i)
# }
#
# gc.combined.seurat$neighbor_cell.type_hvg6 <- NA
# gc.combined.seurat$neighbor_cell.type_hvg6[1:nrow(distance.hvg.list.celltype)] <-
#                         distance.hvg.list.celltype[,2]
#
# pdf("umap query neighbor hvg 6.pdf")
# DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
#   group.by="neighbor_cell.type_hvg6", label=T, repel=T)
# dev.off()
