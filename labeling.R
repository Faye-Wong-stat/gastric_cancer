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
library(harmony)
# library(Biocmanager)
# library(SingleCellExperiment)

setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/integrated/")
# gc.combined <- readRDS("limma.regressed.mnn.integrated.rds")
# # gc.new <- readRDS("gc.new.rds")
# # egc.annt.new <- readRDS("egc.annt.new.rds")
#
# #use clusters
# g <- buildSNNGraph(gc.combined, k=20, use.dimred="corrected")
# clust <- igraph::cluster_walktrap(g)$membership
#
# colLabels(gc.combined) <- factor(clust)
#
# pdf("umap by cluster.pdf")
# plotUMAP(gc.combined, colour_by="label")
# dev.off()




#look at what cells get assigned what clusters
gc.combined.seurat <- readRDS("limma.regressed.mnn.integrated.seurat.rds")

# #finding a better number of neighbors
# gc.combined.seurat <- FindNeighbors(gc.combined.seurat, k.param=15)
# gc.combined.seurat <- FindClusters(gc.combined.seurat)
# pdf("umap seurat by cluster.pdf")
# DimPlot(gc.combined.seurat, reduction="umap", split.by="orig.ident",
#   label=T, repel=T)
# dev.off()
# #15 is orginal


gc.combined.seurat <- FindNeighbors(gc.combined.seurat, k.param=15)
gc.combined.seurat <- FindClusters(gc.combined.seurat)

gc.reference.cluter.celltype <- subset(gc.combined.seurat, subset=orig.ident=="egc.data")@
                                      meta.data[,c("seurat_clusters","cell.type")]
# gc.query.cluter.celltype <- subset(gc.combined.seurat, subset=orig.ident=="gcdata")@
#                                   meta.data[,c("seurat_clusters","cell.type")]

gc.reference.table <- as.data.frame(unclass(table(gc.reference.cluter.celltype$seurat_clusters,
  gc.reference.cluter.celltype$cell.type)))
# gc.query.table <- as.data.frame(unclass(table(gc.query.cluter.celltype$seurat_clusters,
#   gc.query.cluter.celltype$cell.type)))

gc.reference.table <- apply(gc.reference.table, 1, function(x){x/sum(x)})
labels <- names(unlist(apply(gc.reference.table, 2, function(x){which(x>0.8)})))
labels <- data.frame(cluster_cell.type=gsub(".*\\.", "", labels),
                    clusters=gsub("\\..*", "", labels))
rownames(labels) <- labels[,2]

gc.combined.seurat[["cluster_cell.type"]] <- NA
gc.combined.seurat@meta.data[,"cluster_cell.type"] <-
            labels[as.character(gc.combined.seurat@meta.data$seurat_clusters),1]
gc.combined.seurat@meta.data[gc.combined.seurat$orig.ident=="egc.data",
                            "cluster_cell.type"] <- NA

pdf("cluster/umap reference.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="egc.data"), reduction="umap",
  group.by="cell.type", label=T, repel=T)
dev.off()
pdf("cluster/umap query cluster.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="cluster_cell.type", label=T, repel=T)
dev.off()
pdf("cluster/umap cluster.pdf", width=14)
DimPlot(gc.combined.seurat, reduction="umap", split.by="orig.ident",
  label=T, repel=T)
dev.off()

#how many cells didn't get labeled
ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata" &
                  is.na(gc.combined.seurat$cluster_cell.type)]) /
  ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata"])
#38% unlabeled cells

labeled <- as.data.frame(gc.combined.seurat@meta.data[
                              gc.combined.seurat$orig.ident=="gcdata",
                            "cluster_cell.type"])
rownames(labeled) <- colnames(gc.combined.seurat[,gc.combined.seurat$
                              orig.ident=="gcdata"])
names(labeled) <- "15k"

saveRDS(gc.combined.seurat, "limma.regressed.mnn.integrated.seurat.rds")



#tried k.param=15 previously
#try different parameters
gc.combined.seurat <- FindNeighbors(gc.combined.seurat, k.param=10)
gc.combined.seurat <- FindClusters(gc.combined.seurat)
pdf("cluster/umap cluster 10k.pdf", width=14)
DimPlot(gc.combined.seurat, reduction="umap", split.by="orig.ident",
  label=T, repel=T)
dev.off()

gc.reference.cluter.celltype <- subset(gc.combined.seurat, subset=orig.ident=="egc.data")@
                                      meta.data[,c("seurat_clusters","cell.type")]

gc.reference.table <- as.data.frame(unclass(table(gc.reference.cluter.celltype$seurat_clusters,
  gc.reference.cluter.celltype$cell.type)))

gc.reference.table <- apply(gc.reference.table, 1, function(x){x/sum(x)})
labels <- names(unlist(apply(gc.reference.table, 2, function(x){which(x>0.8)})))
labels <- data.frame(cluster_cell.type=gsub(".*\\.", "", labels),
                    clusters=gsub("\\..*", "", labels))
rownames(labels) <- labels[,2]

gc.combined.seurat[["cluster_cell.type"]] <- NA
gc.combined.seurat@meta.data[,"cluster_cell.type"] <-
            labels[as.character(gc.combined.seurat@meta.data$seurat_clusters),1]
gc.combined.seurat@meta.data[gc.combined.seurat$orig.ident=="egc.data",
                            "cluster_cell.type"] <- NA

pdf("cluster/umap reference 10k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="egc.data"), reduction="umap",
  group.by="cell.type", label=T, repel=T)
dev.off()
pdf("cluster/umap query cluster 10k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="cluster_cell.type", label=T, repel=T)
dev.off()

ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata" &
                  is.na(gc.combined.seurat$cluster_cell.type)]) /
  ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata"])
#33% unlabeled cells

labeled[,"10k"] <- gc.combined.seurat@meta.data[
                              gc.combined.seurat$orig.ident=="gcdata",
                            "cluster_cell.type"]



#try k.param=20
gc.combined.seurat <- FindNeighbors(gc.combined.seurat, k.param=20)
gc.combined.seurat <- FindClusters(gc.combined.seurat)
pdf("cluster/umap cluster 20k.pdf", width=14)
DimPlot(gc.combined.seurat, reduction="umap", split.by="orig.ident",
  label=T, repel=T)
dev.off()

gc.reference.cluter.celltype <- subset(gc.combined.seurat, subset=orig.ident=="egc.data")@
                                      meta.data[,c("seurat_clusters","cell.type")]

gc.reference.table <- as.data.frame(unclass(table(gc.reference.cluter.celltype$seurat_clusters,
  gc.reference.cluter.celltype$cell.type)))

gc.reference.table <- apply(gc.reference.table, 1, function(x){x/sum(x)})
labels <- names(unlist(apply(gc.reference.table, 2, function(x){which(x>0.8)})))
labels <- data.frame(cluster_cell.type=gsub(".*\\.", "", labels),
                    clusters=gsub("\\..*", "", labels))
rownames(labels) <- labels[,2]

gc.combined.seurat[["cluster_cell.type"]] <- NA
gc.combined.seurat@meta.data[,"cluster_cell.type"] <-
            labels[as.character(gc.combined.seurat@meta.data$seurat_clusters),1]
gc.combined.seurat@meta.data[gc.combined.seurat$orig.ident=="egc.data",
                            "cluster_cell.type"] <- NA

pdf("cluster/umap reference 20k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="egc.data"), reduction="umap",
  group.by="cell.type", label=T, repel=T)
dev.off()
pdf("cluster/umap query cluster 20k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="cluster_cell.type", label=T, repel=T)
dev.off()

ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata" &
                  is.na(gc.combined.seurat$cluster_cell.type)]) /
  ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata"])
#25% unlabeled cells

labeled[,"20k"] <- gc.combined.seurat@meta.data[
                              gc.combined.seurat$orig.ident=="gcdata",
                            "cluster_cell.type"]



#try k.param=30
gc.combined.seurat <- FindNeighbors(gc.combined.seurat, k.param=30)
gc.combined.seurat <- FindClusters(gc.combined.seurat)
pdf("cluster/umap cluster 30k.pdf", width=14)
DimPlot(gc.combined.seurat, reduction="umap", split.by="orig.ident",
  label=T, repel=T)
dev.off()

gc.reference.cluter.celltype <- subset(gc.combined.seurat, subset=orig.ident=="egc.data")@
                                      meta.data[,c("seurat_clusters","cell.type")]

gc.reference.table <- as.data.frame(unclass(table(gc.reference.cluter.celltype$seurat_clusters,
  gc.reference.cluter.celltype$cell.type)))

gc.reference.table <- apply(gc.reference.table, 1, function(x){x/sum(x)})
labels <- names(unlist(apply(gc.reference.table, 2, function(x){which(x>0.8)})))
labels <- data.frame(cluster_cell.type=gsub(".*\\.", "", labels),
                    clusters=gsub("\\..*", "", labels))
rownames(labels) <- labels[,2]

gc.combined.seurat[["cluster_cell.type"]] <- NA
gc.combined.seurat@meta.data[,"cluster_cell.type"] <-
            labels[as.character(gc.combined.seurat@meta.data$seurat_clusters),1]
gc.combined.seurat@meta.data[gc.combined.seurat$orig.ident=="egc.data",
                            "cluster_cell.type"] <- NA

pdf("cluster/umap reference 30k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="egc.data"), reduction="umap",
  group.by="cell.type", label=T, repel=T)
dev.off()
pdf("cluster/umap query cluster 30k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="cluster_cell.type", label=T, repel=T)
dev.off()

ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata" &
                  is.na(gc.combined.seurat$cluster_cell.type)]) /
  ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata"])
#50% unlabeled cells

labeled[,"30k"] <- gc.combined.seurat@meta.data[
                              gc.combined.seurat$orig.ident=="gcdata",
                            "cluster_cell.type"]



#try k.param=25
gc.combined.seurat <- FindNeighbors(gc.combined.seurat, k.param=25)
gc.combined.seurat <- FindClusters(gc.combined.seurat)
pdf("cluster/umap cluster 25k.pdf", width=14)
DimPlot(gc.combined.seurat, reduction="umap", split.by="orig.ident",
  label=T, repel=T)
dev.off()

gc.reference.cluter.celltype <- subset(gc.combined.seurat, subset=orig.ident=="egc.data")@
                                      meta.data[,c("seurat_clusters","cell.type")]

gc.reference.table <- as.data.frame(unclass(table(gc.reference.cluter.celltype$seurat_clusters,
  gc.reference.cluter.celltype$cell.type)))

gc.reference.table <- apply(gc.reference.table, 1, function(x){x/sum(x)})
labels <- names(unlist(apply(gc.reference.table, 2, function(x){which(x>0.8)})))
labels <- data.frame(cluster_cell.type=gsub(".*\\.", "", labels),
                    clusters=gsub("\\..*", "", labels))
rownames(labels) <- labels[,2]

gc.combined.seurat[["cluster_cell.type"]] <- NA
gc.combined.seurat@meta.data[,"cluster_cell.type"] <-
            labels[as.character(gc.combined.seurat@meta.data$seurat_clusters),1]
gc.combined.seurat@meta.data[gc.combined.seurat$orig.ident=="egc.data",
                            "cluster_cell.type"] <- NA

pdf("cluster/umap reference 25k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="egc.data"), reduction="umap",
  group.by="cell.type", label=T, repel=T)
dev.off()
pdf("cluster/umap query cluster 25k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="cluster_cell.type", label=T, repel=T)
dev.off()

ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata" &
                  is.na(gc.combined.seurat$cluster_cell.type)]) /
  ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata"])
#48% unlabeled cells

labeled[,"25k"] <- gc.combined.seurat@meta.data[
                              gc.combined.seurat$orig.ident=="gcdata",
                            "cluster_cell.type"]



#try k.param=22
gc.combined.seurat <- FindNeighbors(gc.combined.seurat, k.param=22)
gc.combined.seurat <- FindClusters(gc.combined.seurat)
pdf("cluster/umap cluster 22k.pdf", width=14)
DimPlot(gc.combined.seurat, reduction="umap", split.by="orig.ident",
  label=T, repel=T)
dev.off()

gc.reference.cluter.celltype <- subset(gc.combined.seurat, subset=orig.ident=="egc.data")@
                                      meta.data[,c("seurat_clusters","cell.type")]

gc.reference.table <- as.data.frame(unclass(table(gc.reference.cluter.celltype$seurat_clusters,
  gc.reference.cluter.celltype$cell.type)))

gc.reference.table <- apply(gc.reference.table, 1, function(x){x/sum(x)})
labels <- names(unlist(apply(gc.reference.table, 2, function(x){which(x>0.8)})))
labels <- data.frame(cluster_cell.type=gsub(".*\\.", "", labels),
                    clusters=gsub("\\..*", "", labels))
rownames(labels) <- labels[,2]

gc.combined.seurat[["cluster_cell.type"]] <- NA
gc.combined.seurat@meta.data[,"cluster_cell.type"] <-
            labels[as.character(gc.combined.seurat@meta.data$seurat_clusters),1]
gc.combined.seurat@meta.data[gc.combined.seurat$orig.ident=="egc.data",
                            "cluster_cell.type"] <- NA

pdf("cluster/umap reference 22k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="egc.data"), reduction="umap",
  group.by="cell.type", label=T, repel=T)
dev.off()
pdf("cluster/umap query cluster 22k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="cluster_cell.type", label=T, repel=T)
dev.off()

ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata" &
                  is.na(gc.combined.seurat$cluster_cell.type)]) /
  ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata"])
#49% unlabeled cells

labeled[,"22k"] <- gc.combined.seurat@meta.data[
                              gc.combined.seurat$orig.ident=="gcdata",
                            "cluster_cell.type"]



#try k.param=18
gc.combined.seurat <- FindNeighbors(gc.combined.seurat, k.param=18)
gc.combined.seurat <- FindClusters(gc.combined.seurat)
pdf("cluster/umap cluster 18k.pdf", width=14)
DimPlot(gc.combined.seurat, reduction="umap", split.by="orig.ident",
  label=T, repel=T)
dev.off()

gc.reference.cluter.celltype <- subset(gc.combined.seurat, subset=orig.ident=="egc.data")@
                                      meta.data[,c("seurat_clusters","cell.type")]

gc.reference.table <- as.data.frame(unclass(table(gc.reference.cluter.celltype$seurat_clusters,
  gc.reference.cluter.celltype$cell.type)))

gc.reference.table <- apply(gc.reference.table, 1, function(x){x/sum(x)})
labels <- names(unlist(apply(gc.reference.table, 2, function(x){which(x>0.8)})))
labels <- data.frame(cluster_cell.type=gsub(".*\\.", "", labels),
                    clusters=gsub("\\..*", "", labels))
rownames(labels) <- labels[,2]

gc.combined.seurat[["cluster_cell.type"]] <- NA
gc.combined.seurat@meta.data[,"cluster_cell.type"] <-
            labels[as.character(gc.combined.seurat@meta.data$seurat_clusters),1]
gc.combined.seurat@meta.data[gc.combined.seurat$orig.ident=="egc.data",
                            "cluster_cell.type"] <- NA

pdf("cluster/umap reference 18k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="egc.data"), reduction="umap",
  group.by="cell.type", label=T, repel=T)
dev.off()
pdf("cluster/umap query cluster 18k.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="cluster_cell.type", label=T, repel=T)
dev.off()

ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata" &
                  is.na(gc.combined.seurat$cluster_cell.type)]) /
  ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata"])
#43% unlabeled cells

labeled[,"18k"] <- gc.combined.seurat@meta.data[
                              gc.combined.seurat$orig.ident=="gcdata",
                            "cluster_cell.type"]



labeled <- labeled[complete.cases(labeled),]
consistent.labeled <- labeled[apply(labeled, 1, function(x){length(unique(x))==1}),]
consistent.labels <- rownames(consistent.labeled)

#check the umap of cell which got consistent labels
pdf("cluster/umap query consistent.pdf")
DimPlot(gc.combined.seurat[,consistent.labels], reduction="umap",
        group.by="cluster_cell.type", label=T, repel=T)
dev.off()



#k=15, threshold=0.6
gc.combined.seurat <- FindNeighbors(gc.combined.seurat, k.param=15)
gc.combined.seurat <- FindClusters(gc.combined.seurat)

gc.reference.cluter.celltype <- subset(gc.combined.seurat, subset=orig.ident=="egc.data")@
                                      meta.data[,c("seurat_clusters","cell.type")]

gc.reference.table <- as.data.frame(unclass(table(gc.reference.cluter.celltype$seurat_clusters,
  gc.reference.cluter.celltype$cell.type)))

gc.reference.table <- apply(gc.reference.table, 1, function(x){x/sum(x)})
labels <- names(unlist(apply(gc.reference.table, 2, function(x){which(x>0.6)})))
labels <- data.frame(cluster_cell.type=gsub(".*\\.", "", labels),
                    clusters=gsub("\\..*", "", labels))
rownames(labels) <- labels[,2]

gc.combined.seurat[["cluster_cell.type"]] <- NA
gc.combined.seurat@meta.data[,"cluster_cell.type"] <-
            labels[as.character(gc.combined.seurat@meta.data$seurat_clusters),1]
gc.combined.seurat@meta.data[gc.combined.seurat$orig.ident=="egc.data",
                            "cluster_cell.type"] <- NA

pdf("cluster/umap reference 60per.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="egc.data"), reduction="umap",
  group.by="cell.type", label=T, repel=T)
dev.off()
pdf("cluster/umap query cluster 60per.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="cluster_cell.type", label=T, repel=T)
dev.off()
# pdf("cluster/umap cluster.pdf", width=14)
# DimPlot(gc.combined.seurat, reduction="umap", split.by="orig.ident",
#   label=T, repel=T)
# dev.off()

#how many cells didn't get labeled
ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata" &
                  is.na(gc.combined.seurat$cluster_cell.type)]) /
  ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata"])
#7.7% unlabeled cells

labeled15k <- as.data.frame(gc.combined.seurat@meta.data[
                              gc.combined.seurat$orig.ident=="gcdata",
                            "cluster_cell.type"])
rownames(labeled15k) <- colnames(gc.combined.seurat[,gc.combined.seurat$
                              orig.ident=="gcdata"])
names(labeled15k) <- "15k60per"


#threshold=0.9
labels <- names(unlist(apply(gc.reference.table, 2, function(x){which(x>0.9)})))
labels <- data.frame(cluster_cell.type=gsub(".*\\.", "", labels),
                    clusters=gsub("\\..*", "", labels))
rownames(labels) <- labels[,2]

gc.combined.seurat[["cluster_cell.type"]] <- NA
gc.combined.seurat@meta.data[,"cluster_cell.type"] <-
            labels[as.character(gc.combined.seurat@meta.data$seurat_clusters),1]
gc.combined.seurat@meta.data[gc.combined.seurat$orig.ident=="egc.data",
                            "cluster_cell.type"] <- NA

pdf("cluster/umap reference 90per.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="egc.data"), reduction="umap",
  group.by="cell.type", label=T, repel=T)
dev.off()
pdf("cluster/umap query cluster 90per.pdf")
DimPlot(subset(gc.combined.seurat, subset=orig.ident=="gcdata"), reduction="umap",
  group.by="cluster_cell.type", label=T, repel=T)
dev.off()
# pdf("cluster/umap cluster.pdf", width=14)
# DimPlot(gc.combined.seurat, reduction="umap", split.by="orig.ident",
#   label=T, repel=T)
# dev.off()

#how many cells didn't get labeled
ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata" &
                  is.na(gc.combined.seurat$cluster_cell.type)]) /
  ncol(gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata"])
#7.7% unlabeled cells

labeled15k[,"15k90per"] <- gc.combined.seurat@meta.data[
                              gc.combined.seurat$orig.ident=="gcdata",
                            "cluster_cell.type"]










#15k and 30k seemed more accurate, 15k with more labeled cells








#dist order



#look at the mixed cells




# #look at neighbors
# gc.combined.neighbors <- FindNeighbors(gc.combined.seurat,return.neighbor=T,
#                                       graph.name=c("nn", "knn"))
#
# neighbors <- matrix(NA,
#               nrow=ncol(gc.combined.seurat[,gc.combined.seurat[["orig.ident"]]=="gcdata"]),
#               ncol=2)
# rownames(neighbors) <-
#             colnames(gc.combined.seurat[,gc.combined.seurat[["orig.ident"]]=="gcdata"])
#
# for (i in 1:nrow(neighbors)){
#   neighbors[i,] = TopNeighbors
# }


#harmony integration
