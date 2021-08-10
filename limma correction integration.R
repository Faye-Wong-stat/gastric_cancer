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
library(Biocmanager)



setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/early gastric cancer")
gc <- readRDS("../gc.rds")
egc.annt <- readRDS("egc.annt.rds")



#limma
#regress out disease status, linear
egc.annt@assays$RNA@scale.data <- removeBatchEffect(egc.annt@assays$RNA@data,
                                                batch=egc.annt@meta.data$disease.status)

gc@assays$RNA@scale.data <- removeBatchEffect(gc@assays$RNA@data,
                                                batch=gc@meta.data$disease.status)

#combine reference data with our data
shared.gene <- intersect(VariableFeatures(gc), VariableFeatures(egc.annt))
gc.new <- gc[shared.gene,]
egc.annt.new <- egc.annt[shared.gene,]

gc.new@assays$RNA@scale.data <- scale(gc.new@assays$RNA@scale.data)
egc.annt.new@assays$RNA@scale.data <- scale(egc.annt.new@assays$RNA@scale.data)

gc.annt.combined.limma <- merge(gc.new, egc.annt.new)
gc.annt.combined.limma@assays$RNA@scale.data <- cbind(gc.new@assays$RNA@scale.data,
                                  egc.annt.new@assays$RNA@scale.data[rownames(gc.new@assays$RNA@scale.data),])


#regress out dataset difference
gc.annt.combined.limma@assays$RNA@scale.data <- removeBatchEffect(gc.annt.combined.limma@assays$RNA@scale.data,
                                                batch=gc.annt.combined.limma@meta.data$orig.ident)

# gc.annt.combined.limma <- FindVariableFeatures(gc.annt.combined.limma, nfeatures=2000)
gc.annt.combined.limma <- RunPCA(gc.annt.combined.limma, npcs=40, verbose=F,
                                feature=rownames(gc.annt.combined.limma))
gc.annt.combined.limma <- RunUMAP(gc.annt.combined.limma, reduction="pca", dims=1:40)

# pdf("../integration/umap integrated limma.pdf", width=25,height=15)
# p1 <- DimPlot(subset(gc.annt.combined.limma, subset=orig.ident=="egc.data"), reduction="umap",
#               group.by="cell.type", label=T, repel=T) +
#   ggtitle("reference annotation")
# p2 <- DimPlot(subset(gc.annt.combined.limma, subset=orig.ident=="gcdata"), reduction="umap",
#               group.by="predicted.celltype", label=T, repel=T) +
#   ggtitle("transferred labels")
# p1+p2
# dev.off()

pdf("../integration/umap integrated limma test.pdf", width=15,height=15)
DimPlot(gc.annt.combined.limma, reduction="umap", group.by="orig.ident")
dev.off()
pdf("../integration/umap integrated limma.pdf", width=15,height=15)
DimPlot(gc.annt.combined.limma, reduction="umap", group.by="cell.type", split.by="orig.ident", label=T, repel=T)
dev.off()


pdf("../integration/umap integrated limma norm data test.pdf", width=15,height=15)
DimPlot(gc.annt.combined.limma, reduction="umap", group.by="orig.ident")
dev.off()
pdf("../integration/umap integrated limma norm data.pdf", width=15,height=15)
DimPlot(gc.annt.combined.limma, reduction="umap", group.by="cell.type", split.by="orig.ident", label=T, repel=T)
dev.off()
