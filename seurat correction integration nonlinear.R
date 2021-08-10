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



#regress out disease status, non-linear
egc.annt.list <- SplitObject(egc.annt, split.by="disease.status")
egc.annt.list <- lapply(X=egc.annt.list, FUN=function(x){
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
})
features <- SelectIntegrationFeatures(object.list=egc.annt.list)

egc.annt.anchors <- FindIntegrationAnchors(object.list=egc.annt.list, anchor.features=features)
egc.annt.combined <- IntegrateData(anchorset=egc.annt.anchors)

DefaultAssay(egc.annt.combined) <- "integrated"
egc.annt.combined <- ScaleData(egc.annt.combined, verbose=F)
egc.annt.combined <- RunPCA(egc.annt.combined, npcs=30, verbose=F)
egc.annt.combined <- RunUMAP(egc.annt.combined, reduction="pca", dims=1:30)
pdf("umap regress disase status non-l.pdf")
p1 <- DimPlot(egc.annt.combined, reduction = "umap", group.by = "disease status")
p2 <- DimPlot(egc.annt.combined, reduction = "umap", group.by = "cell.type", label = TRUE, repel = TRUE)
p1 + p2
dev.off()

#transfer cell.type data from reference to our data, 2 ways
egc.annt.anchors <- FindTransferAnchors(reference=egc.annt.combined, query=gc, dims=1:30,
                                      reference.reduction="pca")
predictions <- TransferData(anchorset=egc.annt.anchors, refdata=egc.annt.combined$cell.type, dims=1:30)
gc <- AddMetaData(gc, metadata=predictions)

egc.annt.combined <- RunUMAP(egc.annt.combined, reduction="pca", dims=1:30, return.model=T)
gc <- MapQuery(anchorset=egc.annt.anchors, reference=egc.annt.combined, query=gc,
                refdata=list(celltype="cell.type"), reference.reduction="pca", reduction.model="umap")

pdf("../integration/umap integrated.pdf", width=25,height=15)
p1 <- DimPlot(egc.annt.combined, reduction="umap", group.by="cell.type", label=T, repel=T) +
  ggtitle("reference annotation")
p2 <- DimPlot(gc, reduction="ref.umap", group.by="predicted.celltype", label=T, repel=T) +
  ggtitle("transferred labels")
p1+p2
dev.off()

#integrate reference data with ours
shared.gene <- intersect(rownames(gc), rownames(egc.annt.combined))
gc.new <- gc[shared.gene,]
egc.annt.combined.new <- egc.annt.combined[shared.gene,]

gc.combined.nonl.list <- list(egc.annt.combined.new, gc.new)
features <- SelectIntegrationFeatures(object.list=gc.combined.nonl.list)

gc.combined.nonl.anchors <- FindIntegrationAnchors(object.list=gc.combined.nonl.list, anchor.features=features)
gc.combined.nonl <- IntegrateData(anchorset=gc.combined.nonl.anchors)

DefaultAssay(gc.combined.nonl) <- "integrated"
gc.combined.nonl <- ScaleData(gc.combined.nonl, verbose=F)
gc.combined.nonl <- RunPCA(gc.combined.nonl, npcs=30, verbose=F)
gc.combined.nonl <- RunUMAP(gc.combined.nonl, reduction="pca", dims=1:30)
pdf("../integration/umap regress disase status non-l.pdf", width=25,height=15)
p1 <- DimPlot(gc.combined.nonl, reduction = "umap", group.by = "disease.status")
p2 <- DimPlot(gc.combined.nonl, reduction = "umap", group.by = "cell.type", label = TRUE, repel = TRUE)
p1 + p2
dev.off()
pdf("../integration/umap integrated non-l.pdf")
DimPlot(gc.combined.nonl, reduction = "umap", group.by = "orig.ident")
dev.off()
