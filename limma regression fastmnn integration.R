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
# library(Biocmanager)


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
shared.gene <- intersect(rownames(gc), rownames(egc.annt)) #select genes with var>3000
gc.new <- gc[shared.gene,]
egc.annt.new <- egc.annt[shared.gene,]

gc.new@assays$RNA@scale.data <- scale(gc.new@assays$RNA@scale.data)
egc.annt.new@assays$RNA@scale.data <- scale(egc.annt.new@assays$RNA@scale.data)

#do pca and umap
gc.new <- FindVariableFeatures(gc.new)
gc.new <- RunPCA(gc.new, features=VariableFeatures(gc.new))
egc.annt.new <- FindVariableFeatures(egc.annt.new)
egc.annt.new <- RunPCA(egc.annt.new, features=VariableFeatures(egc.annt.new))
gc.new <- RunUMAP(gc.new, dims=1:50)
egc.annt.new <- RunUMAP(egc.annt.new, dims=1:50)
pdf("../integration/regression and integration/umap of gc.pdf")
DimPlot(gc.new, reduction="umap")
dev.off()
pdf("../integration/regression and integration/umap of egc.pdf")
DimPlot(egc.annt.new, reduction="umap", group.by="cell.type", label=T, repel=T)
dev.off()

gc.combined <- fastMNN(gc=gc.new@assays$RNA@scale.data,
                    egc=egc.annt.new@assays$RNA@scale.data[rownames(gc.new@assays$RNA@scale.data),])
gc.combined <- runUMAP(gc.combined, dimred="corrected")

gc.new[["cell.type"]] <- NA
gc.combined$cell.type <- c(gc.new@meta.data[,"cell.type"],
      egc.annt.new@meta.data[,"cell.type"])

#checking if cell types are transferred correctly
# sum(rownames(colData(gc.combined, "cell.type"))==rownames(rbind(gc.new@meta.data[,c("orig.ident","cell.type")],
#       egc.annt.new@meta.data[,c("orig.ident","cell.type")])))
# sum(colData(gc.combined,"cell.type")[,2]==c(gc.new@meta.data[,"cell.type"],
#       egc.annt.new@meta.data[,"cell.type"]), na.rm=T)

#merge two seurat objects
gc.combined.seurat <- merge(gc.new, egc.annt.new)
gc.combined.seurat@assays$RNA@scale.data <- as.matrix(assay(gc.combined,"reconstructed"))
gc.combined.seurat@assays$RNA@data <- cbind(gc.new@assays$RNA@data, egc.annt.new@assays$RNA@data)

#run pca and umap on the merged seurat objects
gc.combined.seurat <- FindVariableFeatures(gc.combined.seurat)
gc.combined.seurat <- RunPCA(gc.combined.seurat,
                            features=VariableFeatures(gc.combined.seurat))
gc.combined.seurat <- RunUMAP(gc.combined.seurat, dims=1:50)
pdf("../integration/regression and integration/umap seurat.pdf")
DimPlot(gc.combined.seurat, reduction="umap", group.by="cell.type", split.by="orig.ident",
  label=T, repel=T)
dev.off()

#plot umap by cluster
gc.combined.seurat <- FindNeighbors(gc.combined.seurat, k.param=15)
gc.combined.seurat <- FindClusters(gc.combined.seurat)
pdf("../integration/regression and integration/umap seurat by cluster.pdf")
DimPlot(gc.combined.seurat, reduction="umap", split.by="orig.ident",
  label=T, repel=T)
dev.off()
#undifferentiated clusters 0, 3, 5, 8, 9, 21

#plot umap with ggplot on the merged seurat objects
gc.combined.seurat.umap <- gc.combined.seurat[["umap"]]@cell.embeddings
gc.combined.seurat.umap <- cbind(gc.combined.seurat.umap,
                          rbind(gc.new@meta.data[,c("orig.ident","cell.type")],
                                egc.annt.new@meta.data[,c("orig.ident","cell.type")]))
colnames(gc.combined.seurat.umap)[1:2] <- c("axisa", "axisb")

pdf("../integration/regression and integration/umap seurat by cell type.pdf",
    width=30, height=30)
ggplot(data=gc.combined.seurat.umap, aes(x=axisa, y=axisb)) +
  geom_point(data=select(gc.combined.seurat.umap, -cell.type), color="grey", size=0.05) +
  geom_point(aes(color=cell.type), size=0.05) +
  facet_wrap("cell.type") +
  guides(colour = guide_legend(override.aes = list(size=1)))
dev.off()



saveRDS(gc.new, "../integrated/gc.new.rds")
saveRDS(egc.annt.new, "../integrated/egc.annt.new.rds")
saveRDS(gc.combined, "../integrated/limma.regressed.mnn.integrated.rds")
saveRDS(gc.combined.seurat, "../integrated/limma.regressed.mnn.integrated.seurat.rds")

gc.new <- readRDS("../integrated/gc.new.rds")
egc.annt.new <- readRDS("../integrated/egc.annt.new.rds")
gc.combined <- readRDS("../integrated/limma.regressed.mnn.integrated.rds")
gc.combined.seurat <- readRDS("../integrated/limma.regressed.mnn.integrated.seurat.rds")



# #plot umap
# pdf("../integration/regression and integration/umap.pdf", width=25, height=10)
# gridExtra::grid.arrange(
#   plotUMAP(gc.combined[,gc.combined$batch=="gc"]),
#   plotUMAP(gc.combined[,gc.combined$batch=="egc"],
#             colour_by=I(egc.annt.new@meta.data$cell.type),
#             text_by=I(egc.annt.new@meta.data$cell.type)),
#   ncol=2
# )
# dev.off()
#
# plot umap with ggplot
# gc.combined.umap <- reducedDim(gc.combined, "UMAP")
# gc.combined.umap <- cbind(gc.combined.umap,
#                           rbind(gc.new@meta.data[,c("orig.ident","cell.type")],
#                                 egc.annt.new@meta.data[,c("orig.ident","cell.type")]))
# colnames(gc.combined.umap)[1:2] <- c("axisa", "axisb")
# # gc.combined.umap[unique(gc.combined.umap$cell.type)] <- NA
# # for (i in 6:dim(gc.combined.umap)[2]){
# #   gc.combined.umap[,i] <- ifelse(gc.combined.umap$cell.type==colnames(gc.combined.umap)[i],1,0)
# # }
#
# pdf("../integration/regression and integration/umap by cell type.pdf", width=30, height=30)
# ggplot(data=gc.combined.umap, aes(x=axisa, y=axisb)) +
#   geom_point(data=select(gc.combined.umap, -cell.type), color="grey", size=0.05) +
#   geom_point(aes(color=cell.type), size=0.05) +
#   facet_wrap("cell.type") +
#   guides(colour = guide_legend(override.aes = list(size=1)))
# dev.off()
# #cancer cell, enterocyte, goblet cell, GMC, MSC, Neck-like cell, PC, PMC,







#check markers and highly variable genes intersection
#use new highly variable genes
gc.selected.combined.seurat <- subset(gc.combined.seurat,
          idents=c("0","3","5","8","9","21","2","4","6","18","19"))

gc.selected.combined.seurat.list <- SplitObject(gc.selected.combined.seurat,
  split.by="orig.ident")
saveRDS(gc.selected.combined.seurat.list, "../integrated/limma.regressed.mnn.integrated.seurat.selected.list.rds")

#find HVGs from the selected cells from both datasets
for (i in 1:2){
  gc.selected.combined.seurat.list[[i]] = FindVariableFeatures(gc.selected.combined.seurat.list[[i]])
  # gc.selected.combined.seurat[[i]] = RunPCA(gc.selected.combined.seurat[[i]],
  #   features=VariableFeatures(gc.selected.combined.seurat[[i]]))
  # gc.selected.combined.seurat[[i]] = RunUMAP(gc.selected.combined.seurat[[i]], dims=1:50)
}

shared.HVGs <- intersect(VariableFeatures(gc.selected.combined.seurat.list[[1]]),
                        VariableFeatures(gc.selected.combined.seurat.list[[2]]))

gc.selected.combined.seurat <- RunPCA(gc.selected.combined.seurat, features=shared.HVGs)
gc.selected.combined.seurat <- RunUMAP(gc.selected.combined.seurat, dims=1:50)

saveRDS(gc.selected.combined.seurat, "../integrated/limma.regressed.mnn.integrated.seurat.selected.rds")

pdf("../integration/regression and integration/umap seurat select.pdf")
DimPlot(gc.selected.combined.seurat, reduction="umap", group.by="cell.type",
  split.by="orig.ident",label=T, repel=T)
dev.off()

gc.selected.combined.seurat.umap <- gc.selected.combined.seurat[["umap"]]@cell.embeddings
gc.selected.combined.seurat.umap <- cbind(gc.selected.combined.seurat.umap,
                        gc.selected.combined.seurat@meta.data[,c("orig.ident","cell.type")])
colnames(gc.selected.combined.seurat.umap)[1:2] <- c("axisa", "axisb")

pdf("../integration/regression and integration/umap seurat select by cell type.pdf",
    width=30, height=30)
ggplot(data=gc.selected.combined.seurat.umap, aes(x=axisa, y=axisb)) +
  geom_point(data=select(gc.selected.combined.seurat.umap, -cell.type), color="grey",
            size=0.05) +
  geom_point(aes(color=cell.type), size=0.05) +
  facet_wrap("cell.type") +
  guides(colour = guide_legend(override.aes = list(size=1)))
dev.off()





# gc.selected.combined <- fastMNN(gc=gc.new.selected@assays$RNA@scale.data,
#   egc=egc.annt.new.selected@assays$RNA@scale.data[rownames(gc.new.selected@assays$RNA@scale.data),])
#
# gc.selected.combined <- runUMAP(gc.selected.combined, dimred="corrected",
#                                 subset_row=top.HVGs.gc.selected.combined)
#
# #plot umap
# pdf("../integration/regression and integration/umap selected.pdf", width=25, height=10)
# gridExtra::grid.arrange(
#   plotUMAP(gc.selected.combined[,gc.selected.combined$batch=="gc"]),
#   plotUMAP(gc.selected.combined[,gc.selected.combined$batch=="egc"],
#             colour_by=I(egc.annt.new.selected@meta.data$cell.type),
#             text_by=I(egc.annt.new.selected@meta.data$cell.type)),
#   ncol=2
# )
# dev.off()
#
# gc.selected.combined.umap <- reducedDim(gc.selected.combined, "UMAP")
#
# gc.new.selected[["cell.type"]] <- "NA"
# gc.selected.combined.umap <- cbind(gc.selected.combined.umap,
#                           rbind(gc.new.selected@meta.data[,c("orig.ident","cell.type")],
#                                 egc.annt.new.selected@meta.data[,c("orig.ident","cell.type")]))
# colnames(gc.selected.combined.umap)[1:2] <- c("axisa", "axisb")
#
# pdf("../integration/regression and integration/umap selected by cell type.pdf", width=30, height=30)
# ggplot(data=gc.selected.combined.umap, aes(x=axisa, y=axisb)) +
#   geom_point(data=select(gc.selected.combined.umap, -cell.type), color="grey", size=0.05) +
#   geom_point(aes(color=cell.type), size=0.05) +
#   facet_wrap("cell.type") +
#   guides(colour = guide_legend(override.aes = list(size=1)))
#
# dev.off()

#transfer labels from reference to ours
#software, or basic stragey/principle
#do it on seurat and on fastmnn


#g to p, phenotyping is hard
#gwas
#image analysis,
#
