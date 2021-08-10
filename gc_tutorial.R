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

getwd()
data <- Read10X(data.dir="/share/quonlab/workspaces/fangyiwang/filtered_feature_bc_matrix/")
#data
gc <- CreateSeuratObject(counts=data, project="gcdata", min.cells=3, min.features=200)
#gc
#dim(gc)

#QC
gc[["percent.mt"]] <- PercentageFeatureSet(gc, pattern="^MT-")
gc[["sample"]] <- gsub('[ACGT]+-','',rownames(gc@meta.data))
# head(gc@meta.data)

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

#Normalization
gc <- NormalizeData(gc, normalization.method="LogNormalize", scale.factor=10000)

#Feature Selection
gc <- FindVariableFeatures(gc, selection.method="vst", nfeatures=2000)
top10 <- head(VariableFeatures(gc), 10)
# plot1 <- VariableFeaturePlot(gc)
# plot2 <- LabelPoints(plot=plot1, points=top10, repel=T, xnudge=0, ynudge=0)
# plot1 + plot2

#Scaling
all.genes <- rownames(gc)
gc <- ScaleData(gc, features=all.genes) #, features=all.genes

#Linear Dimension Reduction
gc <- RunPCA(gc, features=VariableFeatures(object=gc))
# print(gc[["pca"]], dims=1:5, nfeatures=5)
# VizDimLoadings(gc, dims=1:2, reduction="pca")
# DimPlot(gc, reduction="pca") # pc2 looks werid

# DimHeatmap(gc, dims=1, cells=500, balanced=T)

# gc <- JackStraw(gc, num.replicate=100)
# gc <- ScoreJackStraw(gc, dims=1:20)
# JackStrawPlot(gc, dims=1:15)
ElbowPlot(gc)

#Clustering
gc <- FindNeighbors(gc, dims=1:10)
gc <- FindClusters(gc, resolution=0.5)
# head(Idents(gc),5)

#Nonlinear Dimension Reduction
gc <- RunUMAP(gc, dims=1:10)
# DimPlot(gc, reduction="umap", label=T)

gc <- RunTSNE(gc, dims=1:10)
# DimPlot(gc, reduction="tsne", label=T)

#Differentially Expressed Features
gc.markers <- FindAllMarkers(gc, only.pos=T, min.pct=0.25, logfc.threshold=0.25)
gc.markers %>% group_by(cluster) %>% top_n(n=2, wt=avg_log2FC)

VlnPlot(gc, features=c("MS4A1"))
VlnPlot(gc, features=c("NKG7", "PF4"), slot="counts", log=T)
FeaturePlot(gc, features=c("MS4A1", "GNLY"))
gc.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC) -> top10
DoHeatmap(gc, features=top10$gene) + NoLegend()

cluster0.markers <- FindMarkers(gc, ident.1=0, min.pct=0.25)
head(cluster0.markers, n=10)

#Assigning Cell Type Identity to Clusters
levels(gc)
pdf("gc_plots/pca by sample.pdf")
DimPlot(gc, reduction="pca", group.by="sample")
dev.off()
pdf("gc_plots/umap by sample.pdf")
DimPlot(gc, reduction="umap", group.by="sample")
dev.off()
pdf("gc_plots/tsne by sample.pdf")
DimPlot(gc, reduction="tsne", group.by="sample")
dev.off()






## cellranger




#inflammatory bowel disease data

MajorMarkers <- read.csv("ib_plots/Major Markers.csv", skip=1)
iden.markers.auroc <- MajorMarkers[,6:7]
iden.markers.auroc <- iden.markers.auroc[complete.cases(iden.markers.auroc),]
iden.markers.nb <- MajorMarkers[,15:16]
iden.markers.nb <- iden.markers.nb[complete.cases(iden.markers.nb),]

data <- read.table("/share/quonlab/workspaces/fangyiwang/GSE116222_Expression_matrix.txt.gz")

ib <- CreateSeuratObject(counts=data, project="ibdata")
# head(ib@meta.data)

ib <- FindVariableFeatures(ib, selection.method="vst")
top10.ib <- head(VariableFeatures(ib), 10)
# pdf("ib_plots/variable feature.pdf", width=20,height=10)
# plot1 <- VariableFeaturePlot(ib)
# plot2 <- LabelPoints(plot=plot1, points=top10.ib, repel=T, xnudge=0, ynudge=0)
# plot1 + plot2
# dev.off()

all.genes.ib <- rownames(ib)
ib <- ScaleData(ib, features=all.genes.ib)

ib <- RunPCA(ib, features=VariableFeatures(object=ib))
ib <- RunPCA(ib, features=iden.markers.nb$'gene.1')
# ib <- JackStraw(ib, num.replicate=100)
# ib <- ScoreJackStraw(ib, dims=1:20)
pdf("ib_plots/elbow.pdf")
ElbowPlot(ib)
dev.off()
# #18 PC's

ib <- FindNeighbors(ib, dims=1:15)
ib <- FindClusters(ib, resolution=0.5)
ib <- RunUMAP(ib, dims=1:15)
ib <- RunTSNE(ib, dims=1:15)

pdf("ib_plots/pca.pdf")
DimPlot(ib, reduction="pca")
dev.off()
pdf("ib_plots/umap.pdf")
DimPlot(ib, reduction="umap")
dev.off()
pdf("ib_plots/tsne.pdf")
DimPlot(ib, reduction="tsne")
dev.off()

ib.markers <- FindAllMarkers(ib, only.pos=T, min.pct=0.25, logfc.threshold=0.25)
markers <- ib.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC
markers <- as.data.frame(markers)
# write.csv(markers, "ib_plots/markers.csv")

# markers$cell.type <- NA
# table(for (i in 1:300){
#   marker = markers$gene[i]
#   type = iden.markers.auroc[which(iden.markers.auroc$gene==marker),1]
#   print(length(type))
#
# })

cell.type <- unique(iden.markers.auroc$cluster)
top9.marker <- rep(NA, 9)
for (i in 1:10){
  for (j in 1:9){
    top9.marker[j]=as.character(iden.markers.auroc[which(iden.markers.auroc$cluster==cell.type[i])[j],2])
  }
  pdf(paste("ib_plots/featureplot", gsub("/", "-", as.character(cell.type[i])), ".pdf", sep=""))
  print(FeaturePlot(ib, features=top9.marker))
  dev.off()
}

# pdf("ib_plots/featureplot.pdf")
# FeaturePlot(ib, features = c("KRT19", "PHGR1", "SELENBP1"))
# dev.off()
#
# pdf("ib_plots/featureplot2.pdf")
# FeaturePlot(ib, features = c("PCSK1N", "CRYBA2", "SCGN"))
# dev.off()
#
# pdf("ib_plots/featureplot3.pdf")
# FeaturePlot(ib, features = c("TPSB2", "TPSAB1", "VIM"))
# dev.off()



#mucosal profiling paper data
mp.des <- read.table("/share/quonlab/workspaces/fangyiwang/mucosal profiling/GSE121380_series_matrix.txt.gz",
                      nrow=length(readLines("/share/quonlab/workspaces/fangyiwang/mucosal profiling/GSE121380_series_matrix.txt.gz"))-3)



#Crohn's disease paper data
data2009 <- Read10X(data.dir="/share/quonlab/workspaces/fangyiwang/Crohn/2009/")
Cd <- CreateSeuratObject(counts=data2009, project="Cddata", min.cells=3, min.features=200)

#colorectal tumor paper data
crt.id.match <- read.csv("/share/quonlab/workspaces/fangyiwang/colorectal tumor/GSE81861_GEO_EGA_ID_match.csv.gz")
crt.nm.c <- read.csv("/share/quonlab/workspaces/fangyiwang/colorectal tumor/GSE81861_CRC_NM_all_cells_COUNT.csv.gz")
# crt.nm.f <- read.csv("/share/quonlab/workspaces/fangyiwang/colorectal tumor/GSE81861_CRC_NM_all_cells_FPKM.csv.gz")
crt.tumor.c <- read.csv("/share/quonlab/workspaces/fangyiwang/colorectal tumor/GSE81861_CRC_tumor_all_cells_COUNT.csv.gz")
crt.nm.ep.c <- read.csv("/share/quonlab/workspaces/fangyiwang/colorectal tumor/GSE81861_CRC_NM_epithelial_cells_COUNT.csv.gz")
crt.tumor.ep.c <- read.csv("/share/quonlab/workspaces/fangyiwang/colorectal tumor/GSE81861_CRC_tumor_epithelial_cells_COUNT.csv.gz")

colnames(crt.nm.c) <- paste(colnames(crt.nm.c), "nm.all", sep="__")
colnames(crt.nm.ep.c) <- paste(colnames(crt.nm.ep.c), "nm.ep", sep="__")
colnames(crt.tumor.c) <- paste(colnames(crt.tumor.c), "tumor.all", sep="__")
colnames(crt.tumor.ep.c) <- paste(colnames(crt.tumor.ep.c), "tumor.ep", sep="__")

crt.data <- cbind(crt.nm.c,crt.nm.ep.c[,-1])
crt.data <- cbind(crt.data,crt.tumor.c[,-1])
crt.data <- cbind(crt.data,crt.tumor.ep.c[,-1])
rownames(crt.data) <- crt.data[,1]
crt.data <- crt.data[,-1]

crt <- CreateSeuratObject(counts=crt.data, project="crt.data", min.cells=3, min.features=200)

crt@meta.data <- cbind(crt@meta.data, colsplit(rownames(crt@meta.data), "__", c("sample", "cell_type", "other", "source")))
#267 cells from nm.c, 161 nm.pe.c, 376 tumor.c, 273 tumor.ep.c

pdf("crt_plots/violin plot.pdf")
VlnPlot(crt, features=c("nFeature_RNA", "nCount_RNA"), ncol=2, group.by = "orig.ident") #, group.by = "orig.ident"
dev.off()

crt <- NormalizeData(crt, normalization.method="LogNormalize", scale.factor=10000)

crt <- FindVariableFeatures(crt, selection.method="vst", nfeatures=2000)
top10 <- head(VariableFeatures(crt), 10)
plot1 <- VariableFeaturePlot(crt)
plot2 <- LabelPoints(plot=plot1, points=top10, repel=T, xnudge=0, ynudge=0)
pdf("crt_plots/features plot.pdf", width=20, height=10)
plot1 + plot2
dev.off()

all.genes <- rownames(crt)
crt <- ScaleData(crt, features=all.genes) #, features=all.genes

crt <- RunPCA(crt, features=VariableFeatures(object=crt))
# crt <- JackStraw(crt, num.replicate=100)
pdf("crt_plots/elbow plot.pdf")
ElbowPlot(crt)
dev.off()
#10 pcs

crt <- FindNeighbors(crt, dims=1:10)
crt <- FindClusters(crt, resolution=0.5)
crt <- RunUMAP(crt, dims=1:10)
crt <- RunTSNE(crt, dims=1:10, check_duplicates=F)

pdf("crt_plots/pca.pdf")
DimPlot(crt, reduction="pca")
dev.off()
pdf("crt_plots/umap.pdf")
DimPlot(crt, reduction="umap")
dev.off()
pdf("crt_plots/tsne.pdf")
DimPlot(crt, reduction="tsne")
dev.off()

pdf("crt_plots/pca by cell type.pdf")
DimPlot(crt, reduction="pca", group.by="cell_type")
dev.off()
pdf("crt_plots/umap by cell type.pdf")
DimPlot(crt, reduction="umap", group.by="cell_type")
dev.off()
pdf("crt_plots/tsne by cell type.pdf")
DimPlot(crt, reduction="tsne", group.by="cell_type")
dev.off()

pdf("crt_plots/pca by source.pdf")
DimPlot(crt, reduction="pca", group.by="source")
dev.off()
pdf("crt_plots/umap by source.pdf")
DimPlot(crt, reduction="umap", group.by="source")
dev.off()
pdf("crt_plots/tsne by source.pdf")
DimPlot(crt, reduction="tsne", group.by="source")
dev.off()

crt.markers <- FindAllMarkers(crt, only.pos=T, min.pct=0.25, logfc.threshold=0.25)
crt.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC) %>% write.csv("crt_plots/markers.csv")

pdf("crt_plots/pca by other.pdf")
DimPlot(crt, reduction="pca", group.by="other")
dev.off()
pdf("crt_plots/umap by other.pdf")
DimPlot(crt, reduction="umap", group.by="other")
dev.off()
pdf("crt_plots/tsne by other.pdf")
DimPlot(crt, reduction="tsne", group.by="other")
dev.off()

#regress out tumor/normal



#early gastric cancer paper data
# data.nag1 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/ GSM3954946_processed_NAG1.txt.gz")
# data.nag2 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954947_processed_NAG2.txt.gz")
# data.nag3 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954948_processed_NAG3.txt.gz")
# data.cag1 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954949_processed_CAG1.txt.gz")
# data.cag2 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954950_processed_CAG2.txt.gz")
# data.cag3 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954951_processed_CAG3.txt.gz")
# data.imw1 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954952_processed_IMW1.txt.gz")
# data.imw2 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954953_processed_IMW2.txt.gz")
# data.ims1 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954954_processed_IMS1.txt.gz")
# data.ims2 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954955_processed_IMS2.txt.gz")
# data.ims3 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954956_processed_IMS3.txt.gz")
# data.ims4 <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954957_processed_IMS4.txt.gz")
# data.egc <- read.table("/share/quonlab/workspaces/fangyiwang/early gastric cancer/GSM3954958_processed_EGC.txt.gz")
#
# data <- bind_cols(data.nag1, data.nag2, data.nag3, data.cag1, data.cag2, data.cag3, data.imw1, data.imw2,  data.ims1, data.ims2, data.ims3, data.ims4, data.egc)
#
# egc <- CreateSeuratObject(counts=data, project="egc.data", min.cells=3, min.features=200)
# egc[["sample"]] <- NA
# head(egc@meta.data)

setwd("./early gastric cancer")
temp = list.files(pattern="*.gz")
data <- list()
for (i in 1:13){
  data[[i]] = read.table(temp[i])
}

project.names <- c()
for (i in 1:13){
  project.names[i] = gsub(".*_", "", temp[i])
  project.names[i] = gsub("\\..*", "", project.names[i])
}

egc.object <- list()
for (i in 1:13){
  egc.object[[i]] = CreateSeuratObject(counts=data[[i]], project=project.names[i])
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

egc[["sample"]] <- gsub(".*\\.", "", rownames(egc@meta.data))
for(i in 1:dim(egc)[2]){
  egc@meta.data$sample[i] = patient.annt[which(patient.annt$batch==egc.3@meta.data$batch[i]),2]
}
egc[["percent.mt"]] <- PercentageFeatureSet(egc, pattern="^MT-")
egc[["disease.status"]] <- gsub("[1234]", "", egc@meta.data$orig.ident)

pdf("violin plot.pdf",width=20,height=10)
VlnPlot(egc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, group.by="sample") #, group.by = "orig.ident"
dev.off()

pdf("violin plot mt disease status.pdf")
VlnPlot(egc, features="percent.mt", group.by="disease.status") #, group.by = "orig.ident"
dev.off()

#qc following the procedure in the paper
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

patient.annt <- data.frame(batch=project.names, patient=c(1,2,9,3,4,4,5,6,7,7,8,8,8))

egc.3 <- egc[,colnames(egc) %in% rownames(annt)]
pdf("violin plot2.pdf",width=20,height=10)
VlnPlot(egc.3, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3, group.by="sample") #, group.by = "orig.ident"
dev.off()
egc.3[["batch"]] <- annt$batch
egc.3[["cell.type"]] <- annt$cell.type
egc.3[["subject"]] <- NA
for(i in 1:dim(egc.3)[2]){
  egc.3@meta.data$subject[i] = patient.annt[which(patient.annt$batch==egc.3@meta.data$batch[i]),2]
}
egc.3[["disease.status"]] <- gsub("[1234]", "", egc.3@meta.data$batch)
egc.3[["sample"]] <- paste(egc.3@meta.data$batch, egc.3@meta.data$subject, sep=".")

#heatmap celltype by sample/subject, color by proportion, text of percentage
sample.info <- egc.3@meta.data[,c("sample","cell.type")]
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

# subject.info <- egc.3@meta.data[,c("subject","cell.type")]
# subject.info <- as.matrix(data.frame(unclass(table(subject.info))))
# subject.info <- apply(subject.info, 1, FUN=function(x){
#   x/sum(x)
# })
# subject.info <- as.matrix(subject.info)
# subject.info <- t(subject.info)
# subject.info <- round(subject.info, 3)
#
# pdf("heatmap subject.pdf", width=20, height=15)
# image(1:ncol(subject.info), 1:nrow(subject.info), t(subject.info), col=terrain.colors(60), axes=F);
# axis(1, 1:ncol(subject.info), colnames(subject.info));
# axis(2, 1:nrow(subject.info), rownames(subject.info));
# for(x in 1:ncol(subject.info)){
#   for(y in 1:nrow(subject.info)){
#     text(x, y, subject.info[y,x])
#   }
# }
# dev.off()

pdf("violin plot2 mt disease status.pdf")
VlnPlot(egc.3, features="percent.mt", group.by="disease.status") #, group.by = "orig.ident"
dev.off()

pdf("violin plot2 mt cell type.pdf")
VlnPlot(egc.3, features="percent.mt", group.by="cell.type") #, group.by = "orig.ident"
dev.off()

egc.3 <- NormalizeData(egc.3, normalization.method="LogNormalize", scale.factor=10000)
egc.3 <- FindVariableFeatures(egc.3, selection.method="vst", nfeatures=2000)
egc.3 <- ScaleData(egc.3, features=rownames(egc.3))

egc.3 <- RunPCA(egc.3, features=VariableFeatures(object=egc.3))
pdf("elbow plot.pdf")
ElbowPlot(egc.3)
dev.off()

# egc.3 <- FindNeighbors(egc.3, dims=1:11)
# egc.3 <- FindClusters(egc.3, resolution=0.5)
egc.3 <- RunUMAP(egc.3, dims=1:11)
# egc.3 <- RunTSNE(egc.3, dims=1:11) #, check_duplicates=F

pdf("pca.pdf")
DimPlot(egc.3, reduction="pca")
dev.off()
pdf("umap.pdf")
DimPlot(egc.3, reduction="umap")
dev.off()
# pdf("tsne.pdf")
# DimPlot(egc.3, reduction="tsne")
# dev.off()

pdf("pca by patient.pdf")
DimPlot(egc.3, reduction="pca", group.by="subject", label=T, repel=T)
dev.off()
pdf("umap by patient.pdf")
DimPlot(egc.3, reduction="umap", group.by="subject", label=T, repel=T)
dev.off()
# pdf("tsne by patient.pdf")
# DimPlot(egc.3, reduction="tsne", group.by="subject", label=T, repel=T)
# dev.off()

pdf("pca by disease.status.pdf")
DimPlot(egc.3, reduction="pca", group.by="disease.status", label=T, repel=T)
dev.off()
pdf("umap by disease.status.pdf")
DimPlot(egc.3, reduction="umap", group.by="disease.status", label=T, repel=T)
dev.off()
# pdf("tsne by disease.status.pdf")
# DimPlot(egc.3, reduction="tsne", group.by="disease.status", label=T, repel=T)
# dev.off()

pdf("pca by cell type.pdf")
DimPlot(egc.3, reduction="pca", group.by="cell.type", label=T, repel=T)
dev.off()
pdf("umap by cell type.pdf")
DimPlot(egc.3, reduction="umap", group.by="cell.type", label=T, repel=T)
dev.off()
# pdf("tsne by cell type.pdf")
# DimPlot(egc.3, reduction="tsne", group.by="cell.type", label=T, repel=T)
# dev.off()

#regressing out the disease status, non-linear
egc.3.list <- SplitObject(egc.3, split.by="disease.status")
egc.3.list <- lapply(X=egc.3.list, FUN=function(x){
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method="vst", nfeatures=2000) #
})
features <- SelectIntegrationFeatures(object.list=egc.3.list) #

egc.3.anchors <- FindIntegrationAnchors(object.list=egc.3.list, anchor.features=features)
egc.3.combined <- IntegrateData(anchorset=egc.3.anchors)

DefaultAssay(egc.3.combined) <- "integrated"
egc.3.combined <- ScaleData(egc.3.combined, verbose=F)
egc.3.combined <- RunPCA(egc.3.combined, npcs=30, verbose=F)
egc.3.combined <- RunUMAP(egc.3.combined, reduction="pca", dims=1:30)
# egc.3.combined <- RunTSNE(egc.3.combined, reduction="pca", dims=1:30)
# egc.3.combined <- FindNeighbors(egc.3.combined, reduction="pca", dims=1:30)
# egc.3.combined <- FindClusters(egc.3.combined, resolution=0.5)

pdf("regression disease status umap.pdf", width=15, height=8)
p1 <- DimPlot(egc.3.combined, reduction="umap", group.by="disease.status")
p2 <- DimPlot(egc.3.combined, reduction="umap", group.by="cell.type", label=T, repel=T)
p1+p2
dev.off()
pdf("regression disease status umap split.pdf", width=30, height=10)
DimPlot(egc.3.combined, reduction="umap", split.by="disease.status", group.by="cell.type",
        label=T, repel=T)
dev.off()

#another way to regress out disease status, linear
egc.4 <- ScaleData(egc.3, vars.to.regress="disease.status", features=rownames(egc.4))
egc.4 <- RunPCA(egc.4)
egc.4 <- RunUMAP(egc.4, dims=1:30)
pdf("regression disease status umap4.pdf", width=15, height=8)
p1 <- DimPlot(egc.4, reduction="umap", group.by="disease.status")
p2 <- DimPlot(egc.4, reduction="umap", group.by="cell.type", label=T, repel=T)
p1+p2
dev.off()
pdf("regression disease status umap4 split.pdf", width=30, height=10)
DimPlot(egc.4, reduction="umap", split.by="disease.status", group.by="cell.type",
        label=T, repel=T)
dev.off()

shared.gene <- intersect(rownames(gc), rownames(egc.4))
gc.new <- gc[shared.gene,]
egc.4.new <- egc.4[shared.gene,]
gc.4.combined <- merge(gc.new, egc.4.new)
gc.4.combined@assays$RNA@scale.data <- cbind(gc.new@assays$RNA@scale.data, egc.4.new@assays$RNA@scale.data)




#
# egc.5 <- ScaleData(egc.3, vars.to.regress="sample")
# egc.5 <- RunPCA(egc.5)
# egc.5 <- RunUMAP(egc.5, dims=1:30)
pdf("regression sample umap5.pdf", width=15, height=8)
p1 <- DimPlot(egc.5, reduction="umap", group.by="sample")
p2 <- DimPlot(egc.5, reduction="umap", group.by="cell.type", label=T, repel=T)
p1+p2
dev.off()
pdf("regression sample umap split5.pdf", width=30, height=10)
DimPlot(egc.5, reduction="umap", split.by="sample", group.by="cell.type",
        label=T, repel=T, )
dev.off()

#integrating into our data
#repeat the code for our data
data <- Read10X(data.dir="/share/quonlab/workspaces/fangyiwang/filtered_feature_bc_matrix/")
gc <- CreateSeuratObject(counts=data, project="gcdata", min.cells=3, min.features=200)
gc[["percent.mt"]] <- PercentageFeatureSet(gc, pattern="^MT-")
gc[["sample"]] <- gsub('[ACGT]+-','',rownames(gc@meta.data))
gc <- subset(gc, subset= nFeature_RNA<6000& nCount_RNA<10000) # why didn't subset on percent mt
gc <- NormalizeData(gc, normalization.method="LogNormalize", scale.factor=10000)
gc <- FindVariableFeatures(gc, selection.method="vst", nfeatures=2000)
all.genes <- rownames(gc)
gc <- ScaleData(gc, features=all.genes)
gc <- RunPCA(gc, features=rownames(gc))
gc <- RunUMAP(gc, dims=1:30)

#integration
egc.3.anchors <- FindTransferAnchors(reference=egc.3.combined, query=gc, dims=1:30,
                                      reference.reduction="pca")
predictions <- TransferData(anchorset=egc.3.anchors, refdata=egc.3.combined$cell.type, dims=1:30)
gc <- AddMetaData(gc, metadata=predictions)

egc.3.combined <- RunUMAP(egc.3.combined, reduction="pca", dims=1:30, return.model=T)
gc <- MapQuery(anchorset=egc.3.anchors, reference=egc.3.combined, query=gc,
                refdata=list(celltype="cell.type"), reference.reduction="pca", reduction.model="umap")
pdf("../integration/umap integrated.pdf", width=25,height=15)
p1 <- DimPlot(egc.3.combined, reduction="umap", group.by="cell.type", label=T, repel=T) +
  ggtitle("reference annotation")
p2 <- DimPlot(gc, reduction="ref.umap", group.by="predicted.celltype", label=T, repel=T) +
  ggtitle("transferred labels")
p1+p2
dev.off()
#plot them together
#linear regress out disease status () integrate to our data plot
#why use linear instead non-linear

pdf("../integration/test.pdf", width=15,height=15)
DimPlot(egc.3.combined, reduction="umap", group.by="disease.status")
dev.off()
pdf("../integration/test2.pdf", width=15,height=15)
DimPlot(gc, reduction="ref.umap", group.by="sample")
dev.off()
pdf("../integration/test3.pdf", width=15,height=15)
DimPlot(gc, reduction="umap", group.by="sample")
dev.off()
pdf("../integration/test4.pdf", width=15,height=15)
DimPlot(gc, reduction="umap", group.by="predicted.celltype", label=T, repel=T)
dev.off()

#limma
egc.3@assays$RNA@scale.data <- removeBatchEffect(egc.3@assays$RNA@scale.data,
                                                batch=egc.3@meta.data$disease.status)

# gc@assays$RNA@scale.data <- removeBatchEffect(gc@assays$RNA@scale.data,
#                                                 batch=egc.3@meta.data$sample)

shared.gene <- intersect(rownames(gc), rownames(egc.3))
gc.new <- gc[shared.gene,]
egc.3.new <- egc.3[shared.gene,]

gc.combined <- merge(gc.new, egc.3.new)
gc.combined@assays$RNA@scale.data <- cbind(gc.new@assays$RNA@scale.data, egc.3.new@assays$RNA@scale.data)

gc.combined@assays$RNA@scale.data <- removeBatchEffect(gc.combined@assays$RNA@scale.data,
                                                batch=gc.combined@meta.data$orig.ident)
gc.combined <- FindVariableFeatures(gc.combined, nfeatures=2000)
gc.combined <- RunPCA(gc.combined, npcs=30, verbose=F, feature=rownames(gc.combined))
gc.combined <- RunUMAP(gc.combined, reduction="pca", dims=1:30)
# head(gc.combined@meta.data)

pdf("../integration/umap integrated limma.pdf", width=25,height=15)
p1 <- DimPlot(subset(gc.combined, subset=orig.ident=="egc.data"), reduction="umap",
              group.by="cell.type", label=T, repel=T) +
  ggtitle("reference annotation")
p2 <- DimPlot(subset(gc.combined, subset=orig.ident=="gcdata"), reduction="umap",
              group.by="predicted.celltype", label=T, repel=T) +
  ggtitle("transferred labels")
p1+p2
dev.off()

pdf("../integration/umap integrated limma test.pdf", width=15,height=15)
DimPlot(gc.combined, reduction="umap", group.by="orig.ident")
dev.off()

#fastmnn
gc.combined.fastmnn <- merge(gc.new, egc.3.new)
gc.combined.fastmnn@assays$RNA@scale.data <- cbind(gc.new@assays$RNA@scale.data,
                                                  egc.3.new@assays$RNA@scale.data)
# gc.combined.fastmnn <- NormalizeData(gc.combined.fastmnn, normalization.method="LogNormalize", scale.factor=10000)
gc.combined.fastmnn <- FindVariableFeatures(gc.combined.fastmnn, selection.method="vst", nfeatures=2000)
# gc.combined.fastmnn <- ScaleData(gc.combined.fastmnn, features=rownames(gc.combined.fastmnn))
gc.combined.fastmnn <- RunPCA(gc.combined.fastmnn, npcs=30, verbose=F)
gc.combined.fastmnn <- RunUMAP(gc.combined.fastmnn, reduction="pca", dims=1:30)

gc.combined.fastmnn <- RunFastMNN(object.list=SplitObject(gc.combined.fastmnn, split.by="orig.ident"))
gc.combined.fastmnn <- RunUmap(gc.combined.fastmnn, reduction="mnn", dims=1:30)
# gc.combined.fastmnn <- FindNeighbors(gc.combined.fastmnn, reduction="mnn", dims=1:30)
# gc.combined.fastmnn <- FindClusters(gc.combined.fastmnn, reduction="mnn", dims=1:30)

pdf("../integration/umap integrated fast.pdf", width=25,height=15)
p1 <- DimPlot(subset(gc.combined.fastmnn, subset=orig.ident=="egc.data"),
              group.by="cell.type", label=T, repel=T) +
  ggtitle("reference annotation")
p2 <- DimPlot(subset(gc.combined.fastmnn, subset=orig.ident=="gcdata"),
              group.by="cell.type", label=T, repel=T) +
  ggtitle("transferred labels")
p1+p2
dev.off()

pdf("../integration/umap integrated fast test.pdf", width=15,height=15)
DimPlot(gc.combined.fastmnn, reduction="mnn", group.by="orig.ident")
dev.off()


saveRDS(gc, "gc.rds")
saveRDS(egc.3, "egc3.rds")
