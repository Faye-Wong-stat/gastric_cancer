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



#fastmann
shared.gene <- intersect(VariableFeatures(gc), VariableFeatures(egc.annt))
gc.new <- gc[shared.gene,]
egc.annt.new <- egc.annt[shared.gene,]



#remove disease/heath status effect from both data
gc.new.fast <- fastMNN(gc.new@assays$RNA@scale.data, batch=gc.new@meta.data$disease.status)
egc.annt.new.fast <- fastMNN(egc.annt.new@assays$RNA@scale.data,
  batch=egc.annt.new@meta.data$disease.status)

gc.new.fast <- runUMAP(gc.new.fast, dimred="corrected")
egc.annt.new.fast <- runUMAP(egc.annt.new.fast, dimred="corrected")

gc.new.fast.umap <- cbind(reducedDim(gc.new.fast, "UMAP"), gc.new@meta.data$disease.status)
egc.annt.new.fast.umap <- cbind(reducedDim(egc.annt.new.fast, "UMAP"),
  egc.annt.new@meta.data$disease.status)

colnames(gc.new.fast.umap) <- c("axisa", "axisb", "disease_status")
colnames(egc.annt.new.fast.umap) <- c("axisa", "axisb", "disease_status")

gc.new.fast.umap <- as.data.frame(gc.new.fast.umap)
egc.annt.new.fast.umap <- as.data.frame(egc.annt.new.fast.umap)

gc.new.fast.umap[,1:2] <- sapply(gc.new.fast.umap[,1:2], as.numeric)
egc.annt.new.fast.umap[,1:2] <- sapply(egc.annt.new.fast.umap[,1:2], as.numeric)

pdf("../integration/A.pdf")
ggplot(gc.new.fast.umap, aes(x=axisa, y=axisb, col=disease_status)) + geom_point(size=0.05)
dev.off()
pdf("../integration/B.pdf")
ggplot(egc.annt.new.fast.umap, aes(x=axisa, y=axisb, col=disease_status)) + geom_point(size=0.05)
dev.off()

#put two datasets together
gc.combined.fast.corr <- fastMNN(gc=assay(gc.new.fast, "reconstructed"),
  egc=assay(egc.annt.new.fast, "reconstructed")[rownames(assay(gc.new.fast, "reconstructed")),])
gc.combined.fast.corr <- runUMAP(gc.combined.fast.corr, dimred="corrected")

gc.new[["cell.type"]] <- "na"
gc.combined.fast.corr.umap <- cbind(reducedDim(gc.combined.fast.corr, "UMAP"),
  rbind(gc.new@meta.data[,c("orig.ident", "cell.type")],
        egc.annt.new@meta.data[,c("orig.ident", "cell.type")]))
colnames(gc.combined.fast.corr.umap)[1:2] <- c("axisa", "axisb")

pdf("../integration/umap integrated regress batch disease fast.pdf", width=20, height=10)
ggplot(gc.combined.fast.corr.umap, aes(x=axisa, y=axisb, col=cell.type)) +
  geom_point(size=0.05, alpha=0.5) + facet_wrap("orig.ident") +
  guides(colour = guide_legend(override.aes = list(size=10)))
dev.off()

pdf("../integration/umap integrated regress batch disease fast labeled.pdf", width=20, height=10)
gridExtra::grid.arrange(
  plotUMAP(gc.combined.fast.corr[,gc.combined.fast.corr$batch=="gc"]),
  plotUMAP(gc.combined.fast.corr[,gc.combined.fast.corr$batch=="egc"],
            colour_by=I(egc.annt.new@meta.data$cell.type),
            text_by=I(egc.annt.new@meta.data$cell.type)),
  ncol=2
)
dev.off()

unique.cell.types <- unique(egc.annt.new@meta.data$cell.type)
for (i in 1:length(unique.cell.types)){
pdf(paste("../integration/umap integrated regress batch disease fast", unique.cell.types[i],
          "sep.pdf", sep=" "))
print(
  ggplot() +
    geom_point(gc.combined.fast.corr.umap[gc.combined.fast.corr.umap$orig.ident=="egc.data",],
            mapping=aes(x=axisa,y=axisb), color="grey", size=0.1) +
    geom_point(gc.combined.fast.corr.umap[gc.combined.fast.corr.umap$orig.ident=="egc.data" &
                                      gc.combined.fast.corr.umap$cell.type==unique.cell.types[i],],
            mapping=aes(x=axisa,y=axisb), color="red", size=0.1) +
    ggtitle(unique.cell.types[i])
)
print(i)
dev.off()
}

pdf("../integration/umap integrated regress batch disease fast gc.pdf")
ggplot(gc.combined.fast.corr.umap[gc.combined.fast.corr.umap$orig.ident=="gcdata",],
      aes(x=axisa, y=axisb)) +
  geom_point(size=0.5, color="grey")
dev.off()



#regress out disease effect instead of integrate it out
#use limma
#then integrate data with fast mnn

#top cluster cells only run umap on them




gc.combined.fast <- fastMNN(gc=gc.new@assays$RNA@scale.data,
                            egc=egc.annt.new@assays$RNA@scale.data[rownames(gc.new@assays$RNA@scale.data),])

# gc.combined.fast <- runTSNE(gc.combined.fast, dimred="corrected")
# pdf("../integration/tsne integrated fast test.pdf")
# plotTSNE(gc.combined.fast, colour_by="batch")
# dev.off()
gc.combined.fast <- runUMAP(gc.combined.fast, dimred="corrected") #number of pcs
# pdf("../integration/umap integrated fast test.pdf")
# plotUMAP(gc.combined.fast, colour_by="batch")
# dev.off()
pdf("../integration/umap integrated fast test plot.pdf")
plot(reducedDim(gc.combined.fast,"UMAP")[1:dim(gc.new)[2],1],
      reducedDim(gc.combined.fast,"UMAP")[1:dim(gc.new)[2],2],
      col="red")
points(reducedDim(gc.combined.fast,"UMAP")[(dim(gc.new)[2]+1):dim(gc.combined.fast)[2],1],
      reducedDim(gc.combined.fast,"UMAP")[(dim(gc.new)[2]+1):dim(gc.combined.fast)[2],2],
      col="blue")
legend("topleft", c("gc", "egc"), fill=c("red", "blue"))
dev.off()
pdf("../integration/umap integrated fast test plot reverse.pdf")
plot(reducedDim(gc.combined.fast,"UMAP")[(dim(gc.new)[2]+1):dim(gc.combined.fast)[2],1],
      reducedDim(gc.combined.fast,"UMAP")[(dim(gc.new)[2]+1):dim(gc.combined.fast)[2],2],
      col="blue")
points(reducedDim(gc.combined.fast,"UMAP")[1:dim(gc.new)[2],1],
      reducedDim(gc.combined.fast,"UMAP")[1:dim(gc.new)[2],2],
      col="red")
legend("topleft", c("gc", "egc"), fill=c("red", "blue"))
dev.off()



# #reverse the order
# gc.combined.fast2 <- fastMNN(egc=egc.annt.new@assays$RNA@scale.data,
#                             gc=gc.new@assays$RNA@scale.data[rownames(egc.annt.new@assays$RNA@scale.data),],
#                             subset.row=shared.gene)
#
# # gc.combined.fast <- runTSNE(gc.combined.fast, dimred="corrected")
# # pdf("../integration/tsne integrated fast test.pdf")
# # plotTSNE(gc.combined.fast, colour_by="batch")
# # dev.off()
# gc.combined.fast2 <- runUMAP(gc.combined.fast2, dimred="corrected") #number of pcs
# pdf("../integration/umap integrated fast test reversed.pdf")
# plotUMAP(gc.combined.fast2, colour_by="batch")
# dev.off()



#plot umap using seurat
gc.combined.mnn <- merge(gc.new, egc.annt.new)
gc.combined.mnn@assays$RNA@scale.data <-
  as.matrix(assay(gc.combined.fast, "reconstructed")[rownames(gc.combined.mnn@assays$RNA@scale.data),])

# gc.combined.mnn <- FindVariableFeatures(gc.combined.mnn, selection.method="vst", nfeatures=2000,
#                                         features=shared.genes)
gc.combined.mnn <- RunPCA(gc.combined.mnn, features=rownames(gc.combined.mnn))
gc.combined.mnn <- RunUMAP(gc.combined.mnn, dims=1:30)
# pdf("../integration/umpa integrated mnn.pdf",width=25, height=25)
# DimPlot(gc.combined.mnn, reduction="umap", group.by="cell.type", split.by="orig.ident", label=T, repel=T)
# dev.off()
# pdf("../integration/umpa integrated mnn test.pdf",width=25, height=25)
# DimPlot(gc.combined.mnn, reduction="umap", group.by="orig.ident", label=T, repel=T)
# dev.off()

#change the seurat object into SingleCellExperiment object
gc.combined.mnn.sce <- as.SingleCellExperiment(gc.combined.mnn)

reducedDim(gc.combined.mnn.sce, "corrected") <- reducedDim(gc.combined.fast, "corrected")
gc.combined.mnn.sce <- runUMAP(gc.combined.mnn.sce, dimred="corrected")
# pdf("../integration/umap integrated fast test2.pdf")
# plotUMAP(gc.combined.mnn.sce, colour_by="orig.ident")
# dev.off()
pdf("../integration/umap integrated mnn test plot.pdf")
plot(reducedDim(gc.combined.mnn.sce,"UMAP")[1:dim(gc.new)[2],1],
      reducedDim(gc.combined.mnn.sce,"UMAP")[1:dim(gc.new)[2],2],
      col="red")
points(reducedDim(gc.combined.mnn.sce,"UMAP")[(dim(gc.new)[2]+1):dim(gc.combined.mnn.sce)[2],1],
      reducedDim(gc.combined.mnn.sce,"UMAP")[(dim(gc.new)[2]+1):dim(gc.combined.mnn.sce)[2],2],
      col="blue")
legend("topleft", c("gc", "egc"), fill=c("red", "blue"))
dev.off()
pdf("../integration/umap integrated mnn test plot reverse.pdf")
plot(reducedDim(gc.combined.mnn.sce,"UMAP")[(dim(gc.new)[2]+1):dim(gc.combined.mnn.sce)[2],1],
      reducedDim(gc.combined.mnn.sce,"UMAP")[(dim(gc.new)[2]+1):dim(gc.combined.mnn.sce)[2],2],
      col="blue")
points(reducedDim(gc.combined.mnn.sce,"UMAP")[1:dim(gc.new)[2],1],
      reducedDim(gc.combined.mnn.sce,"UMAP")[1:dim(gc.new)[2],2],
      col="red")
legend("topleft", c("gc", "egc"), fill=c("red", "blue"))
dev.off()





# #reverse the order of combining two datasets
# gc.combined.mnn2 <- merge(egc.annt.new, gc.new)
# gc.combined.mnn2@assays$RNA@scale.data <-
#   as.matrix(assay(gc.combined.fast, "reconstructed")[rownames(gc.combined.mnn2@assays$RNA@scale.data),
#                                                       colnames(gc.combined.mnn2@assays$RNA@scale.data)])
#
# gc.combined.mnn2 <- RunPCA(gc.combined.mnn2, features=rownames(gc.combined.mnn2))
# gc.combined.mnn2 <- RunUMAP(gc.combined.mnn2, dims=1:30)
#
# gc.combined.mnn2.sce <- as.SingleCellExperiment(gc.combined.mnn2)
#
# # gc.combined.mnn2.sce <- runUMAP(gc.combined.mnn2.sce, dimred="pca")
# pdf("../integration/umap integrated fast test2 reversed.pdf")
# plotUMAP(gc.combined.mnn2.sce, colour_by="orig.ident")
# dev.off()
