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
library(scPred)
library(magrittr)
library(doParallel)
library(lattice)

setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/integrated/scPred/")

gc.combined.seurat <- readRDS("../limma.regressed.mnn.integrated.seurat.rds")

reference <- gc.combined.seurat[,gc.combined.seurat$orig.ident=="egc.data"]
query <- gc.combined.seurat[,gc.combined.seurat$orig.ident=="gcdata"]

# reference <- reference %>%
#   NormalizeData() %>%
#   FindVariableFeatures() %>%
#   ScaleData() %>%
#   RunPCA() %>%
#   RunUMAP(dims = 1:30)

pdf("umap reference by cell type.pdf")
DimPlot(reference, group.by = "cell.type", label = TRUE, repel = TRUE)
dev.off()



#training
reference <- getFeatureSpace(reference, "cell.type")
reference <- trainModel(reference)
get_probabilities(reference) %>% head()

get_scpred(reference)
pdf("predicted probability.pdf")
plot_probabilities(reference)
dev.off()

reference <- trainModel(reference, model = "mda", reclassify = c("cMono", "ncMono"))
get_scpred(reference)
pdf("predicted probability new.pdf")
plot_probabilities(reference)
dev.off()



#classification
query <- NormalizeData(query)
query <- scPredict(query, reference, recompute_alignment = FALSE)

pdf("scPred/scpred query.pdf")
DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
dev.off()

query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
pdf("scPred/umap query predicted.pdf")
DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
dev.off()
pdf("scPred/umap query original.pdf")
DimPlot(query, group.by = "cell_type", label = TRUE, repel = TRUE)
dev.off()
pdf("scPred/feature plot query.pdf")
FeaturePlot(query, c("scpred_B.cell", "scpred_CD4.T.cell", "scpred_CD8.T.cell",
                     "scpred_cMono", "scpred_ncMono", "scpred_Plasma.cell",
                     "scpred_cDC", "scpred_pDC"))
dev.off()

crossTab(query, "cell_type", "scpred_prediction")
crossTab(query, "cell_type", "scpred_prediction", output = "prop")

# #accessing classifiers
# get_classifiers(reference)
# caret::plot.train(get_classifiers(reference)[["NK cell"]])
# #using a different prediction model approach
# reference <- trainModel(reference, model = "glm")
# get_scpred(reference)

#avoid aligning the data multiple times
query <- scPredict(query, reference, recompute_alignment = FALSE) # *

# #using a different probablity threshold
# query <- scPredict(query, reference, recompute_alignment = FALSE, threshold = 0.9)
# #parallel training
# cl <- makePSOCKcluster(2)
# registerDoParallel(cl)
# reference <- trainModel(reference, model = "mda", allowParallel = TRUE)
# stopCluster(cl)

#applying scPred classifiers without Seurat object
scpred <- get_scpred(reference)
query <- scPredict(query, scpred)

#simpler ways
#
