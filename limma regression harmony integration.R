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
gc.combined.seurat <- merge(gc.new, egc.annt.new)
gc.combined.seurat@assays$RNA@scale.data <- cbind(gc.new@assays$RNA@scale.data,
                                                  egc.annt.new@assays$RNA@scale.data)
