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

setwd("./gastric cancer")

#load our data
gc.data <- Read10X(data.dir="/share/quonlab/workspaces/fangyiwang/gastric cancer/filtered_feature_bc_matrix/")
gc <- CreateSeuratObject(counts=gc.data, project="gcdata", min.cells=3, min.features=200)

#qc
gc[["percent.mt"]] <- PercentageFeatureSet(gc, pattern="^MT-")
gc[["sample"]] <- gsub('[ACGT]+-','',rownames(gc@meta.data))
pdf("gc_plots/violin plot.pdf")
VlnPlot(gc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3) #, group.by = "orig.ident"
dev.off()
# plot1 <- FeatureScatter(gc, feature1="nCount_RNA", feature2="percent.mt")
# plot2 <- FeatureScatter(gc, feature1="nCount_RNA", feature2="nFeature_RNA")
# plot1 + plot2

gc <- subset(gc, subset= nFeature_RNA<6000 & nCount_RNA<10000) # why didn't subset on percent mt
# gc.25 <- subset(gc, subset=nFeature_RNA<6000 & percent.mt<25 & nCount_RNA<10000)
# dim(gc.25)
# dim(subset(gc, subset=nFeature_RNA<6000 & nCount_RNA<5000))

#add more meta data
gc[["disease.status"]] <- NA
gc[["ethnicity"]] <- NA
gc[["biopsy"]] <- NA
gc.patient <- data.frame(sample=1:7,
                        disease.status=c("NG", "NG", "NG", "NG", "MCAG", "MCAG", "MCAGOG"),
                        ethnicity=c(rep("Hispanic",3), rep("White",3), "Hispanic"),
                        biopsy=c("antrum", "antrum", "cardia", "antrum", "antrum", "cardia", "antrum")
                        )
gc@meta.data[,c("disease.status", "ethnicity", "biopsy")] = gc.patient[gc@meta.data$sample,2:4]


#Normalization
gc <- NormalizeData(gc, normalization.method="LogNormalize", scale.factor=1e6)

#Feature Selection
gc <- FindVariableFeatures(gc, selection.method="vst", nfeatures=2000)
# top10 <- head(VariableFeatures(gc), 10)
# plot1 <- VariableFeaturePlot(gc)
# plot2 <- LabelPoints(plot=plot1, points=top10, repel=T, xnudge=0, ynudge=0)
# plot1 + plot2

#Scaling
gc <- ScaleData(gc, features=rownames(gc))

#Linear Dimension Reduction
gc <- RunPCA(gc, features=VariableFeatures(object=gc))
# print(gc[["pca"]], dims=1:5, nfeatures=5)
# VizDimLoadings(gc, dims=1:2, reduction="pca")
pdf("gc_plots/pca plot.pdf")
DimPlot(gc, reduction="pca") # pc2 looks werid
dev.off()

# DimHeatmap(gc, dims=1, cells=500, balanced=T)

# gc <- JackStraw(gc, num.replicate=100)
# gc <- ScoreJackStraw(gc, dims=1:20)
# JackStrawPlot(gc, dims=1:15)
pdf("gc_plots/elbow plot.pdf")
ElbowPlot(gc)
dev.off()

#Clustering
gc <- FindNeighbors(gc, dims=1:10)
gc <- FindClusters(gc, resolution=0.5)
# head(Idents(gc),5)

#Nonlinear Dimension Reduction
gc <- RunUMAP(gc, dims=1:10)
pdf("gc_plots/umap.pdf")
DimPlot(gc, reduction="umap", label=T)
dev.off()
# gc <- RunTSNE(gc, dims=1:10)
# pdf("gc_plots/tsne.pdf")
# DimPlot(gc, reduction="tsne", label=T)
# dev.off()



#load reference data
setwd("./early gastric cancer")
temp = list.files(pattern="*.gz")
data <- list()
for (i in 1:length(temp)){
  data[[i]] = read.table(temp[i])
}

project.names <- c()
for (i in 1:13){
  project.names[i] = gsub(".*_", "", temp[i])
  project.names[i] = gsub("\\..*", "", project.names[i])
}

egc.object <- list()
for (i in 1:13){
  egc.object[[i]] = CreateSeuratObject(counts=data[[i]], project=project.names[i], min.cells=0, min.features=0)
}

egc <- merge(x=egc.object[[1]],y=c(egc.object[[2]],
                                    egc.object[[3]],
                                    egc.object[[4]],
                                    egc.object[[5]],
                                    egc.object[[6]],
                                    egc.object[[7]],
                                    egc.object[[8]],
                                    egc.object[[9]],
                                    egc.object[[10]],
                                    egc.object[[11]],
                                    egc.object[[12]],
                                    egc.object[[13]]))

egc[["sample"]] <- egc[["orig.ident"]]
egc[["orig.ident"]] <- "egc.data"
egc[["disease.status"]] <- gsub("[1234]", "", egc@meta.data$sample)
egc[["percent.mt"]] <- PercentageFeatureSet(egc, pattern="^MT-")

pdf("violin plot.pdf",width=20,height=10)
VlnPlot(egc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, group.by="sample")
#, group.by = "orig.ident"
dev.off()
pdf("violin plot mt disease status.pdf")
VlnPlot(egc, features="percent.mt", group.by="disease.status") #, group.by = "orig.ident"
dev.off()

#qc following the procedure in the paper
#not helpful, don't use
# egc.2 <- subset(egc, subset= nFeature_RNA<=7000 & nFeature_RNA>=400 & percent.mt<=20)
# cell.names.egc.2 <- gsub("\\..*", "", colnames(egc.2))
# duplicated.barcodes <- cell.names.egc.2[duplicated(cell.names.egc.2) |
#                                         duplicated(cell.names.egc.2, fromLast=T)]
# egc.2 <- egc.2[, !cell.names.egc.2 %in% duplicated.barcodes]

#annotation
annt <- read.csv("annotation.csv")
annt$X <- gsub("-", ".", annt$X)
rownames(annt) <- annt$X
annt <- annt[,-1]

patient.annt <- data.frame(sample=project.names, patient=c(1,2,9,3,4,4,5,6,7,7,8,8,8))

egc.annt <- egc[,colnames(egc) %in% rownames(annt)]
egc.annt[["cell.type"]] <- annt$cell.type
egc.annt[["subject"]] <- NA
for(i in 1:dim(egc.annt)[2]){
  egc.annt@meta.data$subject[i] = patient.annt[which(patient.annt$sample==egc.annt@meta.data$sample[i]),2]
}
egc.annt[["sample"]] <- paste(egc.annt@meta.data$sample, egc.annt@meta.data$subject, sep=".")
pdf("violin plot2.pdf",width=20,height=10)
VlnPlot(egc.annt, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, group.by="sample") #, group.by = "orig.ident"
dev.off()

#heatmap celltype by sample
sample.info <- egc.annt@meta.data[,c("sample","cell.type")]
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

pdf("violin plot2 mt disease status.pdf")
VlnPlot(egc.annt, features="percent.mt", group.by="disease.status") #, group.by = "orig.ident"
dev.off()

pdf("violin plot2 mt cell type.pdf")
VlnPlot(egc.annt, features="percent.mt", group.by="cell.type") #, group.by = "orig.ident"
dev.off()

egc.annt <- NormalizeData(egc.annt, normalization.method="LogNormalize", scale.factor=1e6)
egc.annt <- FindVariableFeatures(egc.annt, selection.method="vst", nfeatures=2000)
egc.annt <- ScaleData(egc.annt, features=rownames(egc.annt))

egc.annt <- RunPCA(egc.annt, features=VariableFeatures(object=egc.annt))
pdf("elbow plot.pdf")
ElbowPlot(egc.annt)
dev.off()

# egc.annt <- FindNeighbors(egc.annt, dims=1:11)
# egc.annt <- FindClusters(egc.annt, resolution=0.5)
egc.annt <- RunUMAP(egc.annt, dims=1:15)
# egc.annt <- RunTSNE(egc.annt, dims=1:11) #, check_duplicates=F # didn't run this

pdf("pca.pdf")
DimPlot(egc.annt, reduction="pca")
dev.off()
pdf("umap.pdf")
DimPlot(egc.annt, reduction="umap")
dev.off()
# pdf("tsne.pdf")
# DimPlot(egc.annt, reduction="tsne")
# dev.off()

pdf("pca by patient.pdf")
DimPlot(egc.annt, reduction="pca", group.by="subject", label=T, repel=T)
dev.off()
pdf("umap by patient.pdf")
DimPlot(egc.annt, reduction="umap", group.by="subject", label=T, repel=T)
dev.off()
# pdf("tsne by patient.pdf")
# DimPlot(egc.annt, reduction="tsne", group.by="subject", label=T, repel=T)
# dev.off()

pdf("pca by disease status.pdf")
DimPlot(egc.annt, reduction="pca", group.by="disease.status", label=T, repel=T)
dev.off()
pdf("umap by disease status.pdf")
DimPlot(egc.annt, reduction="umap", group.by="disease.status", label=T, repel=T)
dev.off()
# pdf("tsne by disease.status.pdf")
# DimPlot(egc.annt, reduction="tsne", group.by="disease.status", label=T, repel=T)
# dev.off()

pdf("pca by cell type.pdf")
DimPlot(egc.annt, reduction="pca", group.by="cell.type", label=T, repel=T)
dev.off()
pdf("umap by cell type.pdf")
DimPlot(egc.annt, reduction="umap", group.by="cell.type", label=T, repel=T)
dev.off()
# pdf("tsne by cell type.pdf")
# DimPlot(egc.annt, reduction="tsne", group.by="cell.type", label=T, repel=T)
# dev.off()

saveRDS(gc, "../gc.rds")
saveRDS(egc.annt, "egc.annt.rds")

#regressing out the disease status, non-linear
#don't use
# egc.annt.list <- SplitObject(egc.annt, split.by="disease.status")
# egc.annt.list <- lapply(X=egc.annt.list, FUN=function(x){
#   x = NormalizeData(x)
#   x = FindVariableFeatures(x, selection.method="vst", nfeatures=2000) #
# })
# features <- SelectIntegrationFeatures(object.list=egc.annt.list) #
#
# egc.annt.anchors <- FindIntegrationAnchors(object.list=egc.annt.list, anchor.features=features)
# egc.annt.combined <- IntegrateData(anchorset=egc.annt.anchors)
#
# DefaultAssay(egc.annt.combined) <- "integrated"
# egc.annt.combined <- ScaleData(egc.annt.combined, verbose=F)
# egc.annt.combined <- RunPCA(egc.annt.combined, npcs=30, verbose=F)
# egc.annt.combined <- RunUMAP(egc.annt.combined, reduction="pca", dims=1:30)
# # egc.annt.combined <- RunTSNE(egc.annt.combined, reduction="pca", dims=1:30)
# egc.annt.combined <- FindNeighbors(egc.annt.combined, reduction="pca", dims=1:30)
# egc.annt.combined <- FindClusters(egc.annt.combined, resolution=0.5)
#
# pdf("regression disease status umap.pdf", width=15, height=8)
# p1 <- DimPlot(egc.annt.combined, reduction="umap", group.by="disease.status")
# p2 <- DimPlot(egc.annt.combined, reduction="umap", group.by="cell.type", label=T, repel=T)
# p1+p2
# dev.off()
# pdf("regression disease status umap split.pdf", width=30, height=10)
# DimPlot(egc.annt.combined, reduction="umap", split.by="disease.status", group.by="cell.type",
#         label=T, repel=T)
# dev.off()









# #regress out disease status, linear, same as limma
# #too slow
# egc.annt <- ScaleData(egc.annt, vars.to.regress="disease.status", features=rownames(egc.annt))
# egc.annt <- RunPCA(egc.annt)
# egc.annt <- RunUMAP(egc.annt, dims=1:30)
# pdf("regression disease status umap3.pdf", width=15, height=8)
# p1 <- DimPlot(egc.annt, reduction="umap", group.by="disease.status")
# p2 <- DimPlot(egc.annt, reduction="umap", group.by="cell.type", label=T, repel=T)
# p1+p2
# dev.off()
# pdf("regression disease status umap3 split.pdf", width=30, height=10)
# DimPlot(egc.annt, reduction="umap", split.by="disease.status", group.by="cell.type",
#         label=T, repel=T)
# dev.off()
#
# shared.gene <- intersect(rownames(gc), rownames(egc.annt))
# gc.new <- gc[shared.gene,]
# egc.annt.new <- egc.annt[shared.gene,]
# gc.annt.combined <- merge(gc.new, egc.annt.new)
# gc.annt.combined@assays$RNA@scale.data <- cbind(gc.new@assays$RNA@scale.data, egc.annt.new@assays$RNA@scale.data)
#
# gc.annt.combined <- ScaleData(egc.annt, vars.to.regress="disease.status", features=rownames(egc.annt))




#gene selection difference
#downstream analysis only on our data, not on itegrated data
#simple methods integration
