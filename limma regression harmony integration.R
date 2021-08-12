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

setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/integration/")
gc.new <- readRDS("../integrated/gc.new.rds")
egc.annt.new <- readRDS("../integrated/egc.annt.new.rds")



#merge the two data sets as a seurat object
gc.combined.seurat.har <- merge(gc.new, egc.annt.new)
gc.combined.seurat.har@assays$RNA@scale.data <- cbind(gc.new@assays$RNA@scale.data,
                                                  egc.annt.new@assays$RNA@scale.data)
gc.combined.seurat.har <- FindVariableFeatures(gc.combined.seurat.har)
gc.combined.seurat.har <- RunPCA(gc.combined.seurat.har,
                                features=VariableFeatures(gc.combined.seurat.har))
gc.combined.seurat.har <- RunUMAP(gc.combined.seurat.har, dims=1:50)
pdf("harmony integration/umap.pdf")
DimPlot(gc.combined.seurat.har, reduction="umap", group.by="orig.ident",
        label=T, repel=T)
dev.off()

#harmony integration
pdf("harmony integration/harmony plot.pdf")
gc.combined.seurat.har <- RunHarmony(gc.combined.seurat.har, "orig.ident",
                                    plot_convergence=T)
dev.off()
pdf("harmony integration/harmony.pdf")
DimPlot(gc.combined.seurat.har, reduction="harmony", group.by="orig.ident")
dev.off()

#plot UMAP by harmony instead of pcs
gc.combined.seurat.har <- RunUMAP(gc.combined.seurat.har, reduction="harmony", dims=1:50)
pdf("harmony integration/umap harmony.pdf")
DimPlot(gc.combined.seurat.har, reduction="umap", group.by="orig.ident")
dev.off()
#looks bad

saveRDS(gc.combined.seurat.har, "harmony integration/harmony.integrated.seurat.rds")
