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

write.csv(as.data.frame(unclass(table(gc.combined.seurat[,gc.combined.seurat$orig.ident=="egc.data"]@meta.data$cell.type))), "/share/quonlab/workspaces/fangyiwang/gastric cancer/early gastric cancer/cell type number.csv")
# B cell     Cancer cell      Chief cell              EC      Enterocyte
#   2119             522              44             961            3193
# Enteroendocrine      Fibroblast             GMC     Goblet cell      Macrophage
#   2629            1362            1771             542             395
# Mast cell             MSC  Neck-like cell              PC             PMC
#    227            1243             375             966           12808
# SM cell          T cell
#    272            1617


gc.scaledata.hvg <- gc.combined.seurat@assays$RNA@scale.data[
                                        VariableFeatures(gc.combined.seurat)[1:1000],
                                        gc.combined.seurat$orig.ident=="gcdata"]
egc.scaledata.hvg <- gc.combined.seurat@assays$RNA@scale.data[
                                        VariableFeatures(gc.combined.seurat)[1:1000],
                                        gc.combined.seurat$orig.ident=="egc.data"]
gc.scaledata.hvg <- t(gc.scaledata.hvg)
egc.scaledata.hvg <- t(egc.scaledata.hvg)

# system.time(cdist(gc.scaledata.hvg[1:2,], egc.scaledata.hvg))
# dim(gc.scaledata.hvg)

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
#almost 80% of cells are labeled
saveRDS(gc.combined.seurat,"../limma.regressed.mnn.integrated.neighbor.labeled.seurat.rds")


#this is for presentation
gc.labeled <- gc.combined.seurat[,
              gc.combined.seurat$orig.ident=="gcdata" &
              !is.na(gc.combined.seurat$neighbor_cell.type_hvg10)]
pdf("umap by cell type.pdf")
DimPlot(gc.labeled, reduction="umap",
  group.by="neighbor_cell.type_hvg10", label=T, repel=T)
dev.off()
pdf("umap by sample.pdf")
DimPlot(gc.labeled, reduction="umap",
  group.by="sample")
dev.off()
pdf("umap by ethnicity.pdf")
DimPlot(gc.labeled, reduction="umap",
  group.by="ethnicity")
dev.off()
pdf("umap by biopsy location.pdf")
DimPlot(gc.labeled, reduction="umap",
  group.by="biopsy")
dev.off()

#heatmap celltype by sample
sample.info <- gc.labeled@meta.data[,c("sample","neighbor_cell.type_hvg10")]
sample.info <- as.matrix(data.frame(unclass(table(sample.info))))
sample.info <- apply(sample.info, 1, FUN=function(x){
  x/sum(x)
})
sample.info <- as.matrix(sample.info)
sample.info <- t(sample.info)
sample.info <- round(sample.info, 3)

pdf("heatmap sample.pdf", width=25, height=15)
image(1:ncol(sample.info), 1:nrow(sample.info), t(sample.info), axes=F);
axis(1, 1:ncol(sample.info), colnames(sample.info));
axis(2, 1:nrow(sample.info), rownames(sample.info));
for(x in 1:ncol(sample.info)){
  for(y in 1:nrow(sample.info)){
    text(x, y, sample.info[y,x])
  }
}
dev.off()

Idents(gc.labeled) <- gc.labeled$neighbor_cell.type_hvg10
gc.labeled.markers <- FindAllMarkers(gc.labeled, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25)
gc.labeled.markers.5 <- as.data.frame(gc.labeled.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC))
gc.labeled.markers.2 <- as.data.frame(gc.labeled.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC))
gc.labeled.markers.1 <- as.data.frame(gc.labeled.markers %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = avg_log2FC))

pdf("violin of DEGs.pdf",width=50,height=50)
VlnPlot(gc.labeled, features = gc.labeled.markers.1$gene)
dev.off()

Idents(gc.labeled) <- gc.labeled$ethnicity
gc.labeled.markers.ethnicity <- FindAllMarkers(gc.labeled, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25)
gc.labeled.markers.ethnicity.5 <- as.data.frame(gc.labeled.markers.ethnicity %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC))

pdf("violin of DEGs by ethnicity.pdf",width=50,height=50)
VlnPlot(gc.labeled, features = gc.labeled.markers.ethnicity.5$gene)
dev.off()

Idents(gc.labeled) <- gc.labeled$biopsy
gc.labeled.markers.biopsy <- FindAllMarkers(gc.labeled, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25)
gc.labeled.markers.biopsy.5 <- as.data.frame(gc.labeled.markers.biopsy %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC))

pdf("violin of DEGs by biopsy.pdf",width=50,height=50)
VlnPlot(gc.labeled, features = gc.labeled.markers.biopsy.5$gene)
dev.off()

gc.labeled$disease.status[gc.labeled$disease.status=="MCAGOG"] <- "MCAG"
Idents(gc.labeled) <- gc.labeled$disease.status
gc.labeled.markers.disease <- FindAllMarkers(gc.labeled, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25)
gc.labeled.markers.disease.5 <- as.data.frame(gc.labeled.markers.disease %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC))

pdf("violin of DEGs by disease.pdf",width=50,height=50)
VlnPlot(gc.labeled, features = gc.labeled.markers.disease.5$gene)
dev.off()




#subset on cell types
#stratify on different levels and see how genes differently express




saveRDS(gc.labeled, "gc.labeled.rds")
saveRDS(gc.labeled.markers, "gc.labeled.markers.rds")
gc.labeled <- readRDS("gc.labeled.rds")
gc.labeled.markers <- readRDS("gc.labeled.markers.rds")

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
