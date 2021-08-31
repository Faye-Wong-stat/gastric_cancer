library(tidyr)
library(dplyr)
library(plyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(grid)
library(limma)
library(batchelor)
library(scran)
library(scater)
# library(Biocmanager)

setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer")

#load our data
gc.data <- Read10X(data.dir="/share/quonlab/workspaces/fangyiwang/gastric cancer/filtered_feature_bc_matrix/")
gc <- CreateSeuratObject(counts=gc.data, project="gcdata", min.cells=3, min.features=200)

#qc
gc[["percent.mt"]] <- PercentageFeatureSet(gc, pattern="^MT-")
gc[["sample"]] <- gsub('[ACGT]+-','',rownames(gc@meta.data))
dim(gc)
#27144 35775
# plot1 <- FeatureScatter(gc, feature1="nCount_RNA", feature2="percent.mt")
# plot2 <- FeatureScatter(gc, feature1="nCount_RNA", feature2="nFeature_RNA")
# plot1 + plot2

# gc.25 <- subset(gc, subset=nFeature_RNA<6000 & percent.mt<25 & nCount_RNA<10000)
# dim(gc.25)
# dim(subset(gc, subset=nFeature_RNA<6000 & nCount_RNA<5000))

#add more meta data
gc[["disease.status"]] <- NA
gc[["ethnicity"]] <- NA
gc[["biopsy"]] <- NA
gc[["subject"]] <- NA
gc.patient <- data.frame(sample=1:7,
                        disease.status=c("NG", "NG", "NG", "NG", "MCAG", "MCAG", "MCAG"),
                        ethnicity=c(rep("Hispanic",3), rep("White",3), "Hispanic"),
                        biopsy=c("antrum", "antrum", "cardia", "antrum", "antrum", "cardia", "antrum"),
                        subject=c(2,3,3,4,6,6,7))
gc@meta.data[,c("disease.status", "ethnicity", "biopsy","subject")] = gc.patient[gc@meta.data$sample,2:5]

pdf("gc_plots/violin plot before QC.pdf", width=14)
VlnPlot(gc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#, group.by = "orig.ident"
dev.off()
pdf("gc_plots/violin plot before QC by disease status.pdf", width=14)
VlnPlot(gc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3,
        group.by = "disease.status")
dev.off()

# gc.adj <- gc[,((gc$disease.status=="MCAG"&gc$nFeature_RNA<9000) |
#                 (gc$disease.status=="NG"&gc$nFeature_RNA<8000)) &
#               ((gc$disease.status=="MCAG"&gc$nCount_RNA<8500) |
#                 (gc$disease.status=="NG"&gc$nCount_RNA<7500))]
pdf("gc_plots/violin plot adjusting.pdf", width=14)
VlnPlot(gc.adj,
        features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3,
        group.by = "disease.status")
dev.off()

gc <- gc[,((gc$disease.status=="MCAG"&gc$nFeature_RNA<9000) |
                (gc$disease.status=="NG"&gc$nFeature_RNA<8000)) &
              ((gc$disease.status=="MCAG"&gc$nCount_RNA<8500) |
                (gc$disease.status=="NG"&gc$nCount_RNA<7500))]
gc <- gc[!rowSums(gc@assays$RNA@counts)==0,]
dim(gc)
#26963 29210

# gc <- subset(gc, subset= nFeature_RNA<6000 & nCount_RNA<10000)
# #why didn't subset on percent mt

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
gc <- RunUMAP(gc, dims=1:50)
pdf("gc_plots/umap.pdf")
DimPlot(gc, reduction="umap", label=T)
dev.off()
# gc <- RunTSNE(gc, dims=1:10)
# pdf("gc_plots/tsne.pdf")
# DimPlot(gc, reduction="tsne", label=T)
# dev.off()



#load reference data
setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/early gastric cancer")
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

egc[["percent.mt"]] <- PercentageFeatureSet(egc, pattern="^MT-")
egc[["sample"]] <- egc[["orig.ident"]]
egc[["orig.ident"]] <- "egc.data"
egc[["disease.status"]] <- gsub("[1234]", "", egc@meta.data$sample)
dim(egc)
# 22910 56440

pdf("egc_plots/violin plot before QC.pdf",width=14)
VlnPlot(egc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3,
        group.by="orig.ident")
dev.off()
pdf("egc_plots/violin plot before QC by disease status.pdf",width=14)
VlnPlot(egc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by="disease.status")
dev.off()

# egc.adj <- egc[,((egc$disease.status=="EGC"&egc$nFeature_RNA<4000) |
#                  (egc$disease.status=="IMS"&egc$nFeature_RNA<3700) |
#                  (egc$disease.status=="CAG"&egc$nFeature_RNA<3000) |
#                  (egc$disease.status=="NAG"&egc$nFeature_RNA<2500) |
#                  (egc$disease.status=="IMW"&egc$nFeature_RNA<2000)) &
#                 ((egc$disease.status%in%c("EGC","IMS")&egc$nCount_RNA<20000) |
#                  (egc$disease.status=="CAG"&egc$nCount_RNA<12000) |
#                  (egc$disease.status=="NAG"&egc$nCount_RNA<10000) |
#                  (egc$disease.status=="IMW"&egc$nCount_RNA<7500)) &
#                 ((egc$disease.status=="EGC"&egc$percent.mt<75) |
#                  (egc$disease.status=="IMS"&egc$percent.mt<80) |
#                  (egc$disease.status=="CAG"&egc$percent.mt<45) |
#                  (egc$disease.status=="IMW"&egc$percent.mt<50) |
#                  (egc$disease.status=="NAG"&egc$percent.mt<30))]
# pdf("egc_plots/violin plot adjusting.pdf",width=14)
# VlnPlot(egc.adj, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#         group.by="disease.status")
# dev.off()

#qc following the procedure in the paper
#not helpful, don't use
# egc.2 <- subset(egc, subset= nFeature_RNA<=7000 & nFeature_RNA>=400 & percent.mt<=20)
# cell.names.egc.2 <- gsub("\\..*", "", colnames(egc.2))
# duplicated.barcodes <- cell.names.egc.2[duplicated(cell.names.egc.2) |
#                                         duplicated(cell.names.egc.2, fromLast=T)]
# egc.2 <- egc.2[, !cell.names.egc.2 %in% duplicated.barcodes]

egc <- egc[,((egc$disease.status=="EGC"&egc$nFeature_RNA<4000) |
                 (egc$disease.status=="IMS"&egc$nFeature_RNA<3700) |
                 (egc$disease.status=="CAG"&egc$nFeature_RNA<3000) |
                 (egc$disease.status=="NAG"&egc$nFeature_RNA<2500) |
                 (egc$disease.status=="IMW"&egc$nFeature_RNA<2000)) &
                ((egc$disease.status%in%c("EGC","IMS")&egc$nCount_RNA<20000) |
                 (egc$disease.status=="CAG"&egc$nCount_RNA<12000) |
                 (egc$disease.status=="NAG"&egc$nCount_RNA<10000) |
                 (egc$disease.status=="IMW"&egc$nCount_RNA<7500)) &
                ((egc$disease.status=="EGC"&egc$percent.mt<75) |
                 (egc$disease.status=="IMS"&egc$percent.mt<80) |
                 (egc$disease.status=="CAG"&egc$percent.mt<45) |
                 (egc$disease.status=="IMW"&egc$percent.mt<50) |
                 (egc$disease.status=="NAG"&egc$percent.mt<30))]
egc <- egc[!rowSums(egc@assays$RNA@counts)==0,]
dim(egc)
#22904 53727

#annotation
annt <- read.csv("annotation.csv")
annt$X <- gsub("-", ".", annt$X)
rownames(annt) <- annt$X
annt <- annt[,-1]

egc.annt <- egc

patient.annt <- data.frame(sample=project.names,
                            age=c(58,56,62,51,62,62,63,48,rep(68,2),rep(67,3)),
                            sex=c("m","f","m","m","f","f","m","f",rep("m",5)),
                            Hpylori=c(rep("N",6),"P","N",rep("P",2),rep("N",3)),
                            patient=c(1,2,9,3,4,4,5,6,7,7,8,8,8))

egc.annt[["age"]] <- NA
egc.annt[["sex"]] <- NA
egc.annt[["Hpylori"]] <- NA
egc.annt[["subject"]] <- NA
for(i in 1:dim(egc.annt)[2]){
  egc.annt@meta.data[i,c("age","sex","Hpylori","subject")] = patient.annt[which(patient.annt$sample==egc.annt@meta.data$sample[i]),2:5]
}
# egc.annt[["sample"]] <- paste(egc.annt@meta.data$sample, egc.annt@meta.data$subject, sep=".")

egc.annt[["cell.type"]] <- NA
egc.annt$cell.type <- annt[colnames(egc.annt),2]

dim(egc.annt)
#22904 53727



# pdf("egc_plots/violin plot2 before QC.pdf",width=14)
# VlnPlot(egc.annt, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3,
#         group.by="orig.ident")
# dev.off()
# pdf("egc_plots/violin plot2 before QC by disease status.pdf",width=14)
# VlnPlot(egc.annt, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3,
#         group.by="disease.status")
# dev.off()

#heatmap celltype by sample
sample.info <- egc.annt@meta.data[,c("sample","cell.type")] #!is.na(egc.annt@meta.data$cell.type)
sample.info$cell.type[is.na(sample.info$cell.type)] <- "NA"
sample.info <- as.matrix(data.frame(unclass(table(sample.info))))
sample.info <- apply(sample.info, 1, FUN=function(x){
  x/sum(x)
})
sample.info <- as.matrix(sample.info)
sample.info <- t(sample.info)
sample.info <- round(sample.info, 3)

pdf("egc_plots/heatmap sample.pdf",width=14,height=14)
image(1:ncol(sample.info), 1:nrow(sample.info), t(sample.info), axes=F,
      xlab="cell type", ylab="sample");
axis(1, 1:ncol(sample.info), labels=F);
axis(2, 1:nrow(sample.info), rownames(sample.info), las=2);
for(x in 1:ncol(sample.info)){
  for(y in 1:nrow(sample.info)){
    text(x, y, sample.info[y,x])
  }
}
text(x=1:ncol(sample.info), y=par("usr")[3]-0.45, xpd=NA,
  labels=colnames(sample.info),srt=25)
dev.off()

# pdf("violin plot2 mt disease status.pdf")
# VlnPlot(egc.annt, features="percent.mt", group.by="disease.status") #, group.by = "orig.ident"
# dev.off()
#
pdf("egc_plots/violin plot2 mt cell type.pdf")
VlnPlot(egc.annt[,!is.na(egc.annt$cell.type)],
        features="percent.mt", group.by="cell.type") #, group.by = "orig.ident"
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
egc.annt <- RunUMAP(egc.annt, dims=1:50)
# egc.annt <- RunTSNE(egc.annt, dims=1:11) #, check_duplicates=F # didn't run this

pdf("egc_plots/umap.pdf")
DimPlot(egc.annt, reduction="umap")
dev.off()
pdf("egc_plots/umap by patient.pdf")
DimPlot(egc.annt, reduction="umap", group.by="subject", label=T, repel=T)
dev.off()
pdf("egc_plots/umap by disease status.pdf")
DimPlot(egc.annt, reduction="umap", group.by="disease.status", label=T, repel=T)
dev.off()
pdf("egc_plots/umap by cell type.pdf")
DimPlot(egc.annt, reduction="umap", group.by="cell.type", label=T, repel=T)
dev.off()



#exclude patient 7 and 8
egc.annt.exc78 <- subset(egc.annt, subset=subject%in%c(1:6,9))
egc.annt.exc78 <- egc.annt.exc78[!rowSums(egc.annt.exc78@assays$RNA@counts)==0,]
dim(egc.annt.exc78)
#22597 33210

pdf("egc_plots/violin plot exclude 7&8.pdf",width=14)
VlnPlot(egc.annt.exc78, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by="disease.status")
dev.off()

egc.annt.exc78 <- NormalizeData(egc.annt.exc78, normalization.method="LogNormalize", scale.factor=1e6)
egc.annt.exc78 <- FindVariableFeatures(egc.annt.exc78, selection.method="vst", nfeatures=2000)
egc.annt.exc78 <- ScaleData(egc.annt.exc78, features=rownames(egc.annt.exc78))
egc.annt.exc78 <- RunPCA(egc.annt.exc78, features=VariableFeatures(object=egc.annt.exc78))
egc.annt.exc78 <- RunUMAP(egc.annt.exc78, dims=1:50)
pdf("egc_plots/umap exclude 7&8.pdf")
DimPlot(egc.annt.exc78, reduction="umap")
dev.off()
pdf("egc_plots/umap by patient exclude 7&8.pdf")
DimPlot(egc.annt.exc78, reduction="umap", group.by="subject", label=T, repel=T)
dev.off()
pdf("egc_plots/umap by disease status exclude 7&8.pdf")
DimPlot(egc.annt.exc78, reduction="umap", group.by="disease.status", label=T, repel=T)
dev.off()
pdf("egc_plots/umap by cell type exclude 7&8.pdf")
DimPlot(egc.annt.exc78, reduction="umap", group.by="cell.type", label=T, repel=T)
dev.off()




# pdf("pca.pdf")
# DimPlot(egc.annt, reduction="pca")
# dev.off()
# pdf("umap.pdf")
# DimPlot(egc.annt, reduction="umap")
# dev.off()
# # pdf("tsne.pdf")
# # DimPlot(egc.annt, reduction="tsne")
# # dev.off()
#
# pdf("pca by patient.pdf")
# DimPlot(egc.annt, reduction="pca", group.by="subject", label=T, repel=T)
# dev.off()
# pdf("umap by patient.pdf")
# DimPlot(egc.annt, reduction="umap", group.by="subject", label=T, repel=T)
# dev.off()
# # pdf("tsne by patient.pdf")
# # DimPlot(egc.annt, reduction="tsne", group.by="subject", label=T, repel=T)
# # dev.off()
#
# pdf("pca by disease status.pdf")
# DimPlot(egc.annt, reduction="pca", group.by="disease.status", label=T, repel=T)
# dev.off()
# pdf("umap by disease status.pdf")
# DimPlot(egc.annt, reduction="umap", group.by="disease.status", label=T, repel=T)
# dev.off()
# # pdf("tsne by disease.status.pdf")
# # DimPlot(egc.annt, reduction="tsne", group.by="disease.status", label=T, repel=T)
# # dev.off()
#
# pdf("pca by cell type.pdf")
# DimPlot(egc.annt, reduction="pca", group.by="cell.type", label=T, repel=T)
# dev.off()
# pdf("umap by cell type.pdf")
# DimPlot(egc.annt, reduction="umap", group.by="cell.type", label=T, repel=T)
# dev.off()
# # pdf("tsne by cell type.pdf")
# # DimPlot(egc.annt, reduction="tsne", group.by="cell.type", label=T, repel=T)
# # dev.off()

saveRDS(gc, "../gc.rds")
saveRDS(egc.annt, "egc.annt.rds")
saveRDS(egc.annt.exc78, "egc.annt.exc78.rds")

gc <- readRDS("../gc.rds")
egc.annt <- readRDS("egc.annt.rds")
egc.annt.exc78 <- readRDS("egc.annt.exc78.rds")



#plotting for presentation
# pdf("violin plot.pdf",width=20,height=10)
# VlnPlot(egc.annt, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3,
#         group.by="orig.ident") #, group.by = "orig.ident"
# dev.off()
# pdf("violin plot mt disease status.pdf")
# VlnPlot(egc.annt, features="percent.mt", group.by="disease.status")
# dev.off()

setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/gc_plots")
gc.heatmap <- gc@meta.data[,c("sample","disease.status","ethnicity","biopsy","subject")]
gc.heatmap$sample <- as.character(gc.heatmap$sample)
gc.heatmap$subject <- as.character(gc.heatmap$subject)
gc.heatmap <- plyr::count(gc.heatmap)
gc.heatmap <- gc.heatmap[order(gc.heatmap$freq,decreasing=T),]
rownames(gc.heatmap) <- 1:7

gc.heatmap.trim <- gc.heatmap
gc.heatmap.trim <- gc.heatmap.trim[,2:4]
vector.0 <- gc.heatmap.trim[1,]
# vector.1 <- c()
# for (i in 1:length(vector.0)){
#   vector.1[i] <- gc.heatmap.trim[gc.heatmap.trim[,i]!=vector.0[i],i][1]
# }
gc.heatmap.trim <- as.data.frame(apply(gc.heatmap.trim,2,FUN=function(x){ifelse(x==x[1],1,0)}))
rownames(gc.heatmap.trim) <- gc.heatmap$sample
gc.heatmap.trim$sample <- rownames(gc.heatmap.trim)
gc.heatmap.trim <- gather(gc.heatmap.trim,key="y",value="bi",disease.status:biopsy)
ordered <- gc.heatmap$sample

plot1 <- ggplot(gc.heatmap, aes(x=reorder(sample,-freq), y=freq)) +
  geom_bar(stat="identity") +
  xlab("sample") + ylab("cell_counts") + coord_flip() +
  theme(axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),)

plot2 <- ggplot(gc.heatmap.trim,
                aes(x=y,y=factor(sample,levels=ordered),fill=as.logical(bi))) +
         scale_x_discrete(labels=vector.0) +
         xlab("") + ylab("sample") + guides(fill=guide_legend(title="")) +
         geom_tile() + theme(legend.position="left")
pdf("descriptive analysis.pdf",width=14)
grid.arrange(plot2,plot1,ncol=2)
# grid.draw(rbind(ggplotGrob(plot1),ggplotGrob(plot2), size="first"))
dev.off()

setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/early gastric cancer/egc_plots")
egc.annt.heatmap <- egc.annt@meta.data[,c("sample","sex","Hpylori","subject")]
egc.annt.heatmap$subject <- as.character(egc.annt.heatmap$subject)
egc.annt.heatmap <- plyr::count(egc.annt.heatmap)
egc.annt.heatmap <- egc.annt.heatmap[order(egc.annt.heatmap$freq,decreasing=T),]
rownames(egc.annt.heatmap) <- 1:13

egc.annt.heatmap.trim <- egc.annt.heatmap
egc.annt.heatmap.trim <- egc.annt.heatmap.trim[,2:3]
vector.0 <- egc.annt.heatmap.trim[1,]

egc.annt.heatmap.trim <- as.data.frame(apply(egc.annt.heatmap.trim,2,FUN=function(x){ifelse(x==x[1],0,1)}))
rownames(egc.annt.heatmap.trim) <- egc.annt.heatmap$sample
egc.annt.heatmap.trim$sample <- rownames(egc.annt.heatmap.trim)
egc.annt.heatmap.trim <- gather(egc.annt.heatmap.trim,key="y",value="bi",sex:Hpylori)
ordered <- egc.annt.heatmap$sample

plot1 <- ggplot(egc.annt.heatmap, aes(x=reorder(sample,-freq), y=freq)) +
  geom_bar(stat="identity") +
  xlab("sample") + ylab("cell_counts") + coord_flip() +
  theme(axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),)

plot2 <- ggplot(egc.annt.heatmap.trim,
                aes(x=y,y=factor(sample,levels=ordered),fill=as.logical(bi))) +
         scale_x_discrete(labels=c("H.pylori Infection","Male")) +
         xlab("") + ylab("sample") + guides(fill=guide_legend(title="")) +
         geom_tile() + theme(legend.position="left")
pdf("descriptive analysis.pdf",width=14)
grid.arrange(plot2,plot1,ncol=2)
# grid.draw(rbind(ggplotGrob(plot1),ggplotGrob(plot2), size="first"))
dev.off()

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
