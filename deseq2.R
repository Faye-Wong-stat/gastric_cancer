library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(ashr)
library(png)
library(DESeq2)
library(RColorBrewer)
# BiocManager::install(version = '3.12')


setwd("/share/quonlab/workspaces/fangyiwang/gastric cancer/gc_plots")
gc.combined.exc78.seurat <- readRDS("../integrated/neighbor/limma.regressed.mnn.integrated.neighbor.labeled.exc78.seurat.rds")
gc.labeled <- gc.combined.exc78.seurat[,
              gc.combined.exc78.seurat$orig.ident=="gcdata" &
              !is.na(gc.combined.exc78.seurat$neighbor_cell.type_hvg10)]



#try deseq2
counts <- gc.labeled@assays$RNA@counts
metadata <- gc.labeled@meta.data
gc.labeled.sce <- SingleCellExperiment(assays=list(counts=counts),colData=metadata)
# gc.labeled.sce.meta <- colData(gc.labeled.sce)[,
# c("sample","disease.status","ethnicity","biopsy","subject")]
# gc.labeled.sce.meta.sample <- plyr::count(gc.labeled.sce.meta)



kids <- purrr::set_names(levels(as.factor(gc.labeled.sce$neighbor_cell.type_hvg10)))
nk <- length(kids)
sids <- purrr::set_names(levels(as.factor(gc.labeled.sce$sample)))
ns <- length(sids)

table(gc.labeled.sce$sample)
# 1    2    3    4    5    6    7
# 3675 1472 4823 1688 4793 5518 3635
table(gc.labeled.sce$neighbor_cell.type_hvg10,gc.labeled.sce$sample)
#                    1    2    3    4    5    6    7
# B cell             8   16  426   44    1  104   63
# EC                20   50  111   44   25   92   60
# Enterocyte         1    0    0   43    0    0    0
# Enteroendocrine  121   40  236   25  157  157   81
# Fibroblast        32    9   42   22    4   40   15
# GMC              254  196  602  142  280  187  258
# Macrophage        17   28   39   24   15   16   12
# Mast cell         54   77  424   59   61  488  108
# PC                67    8   14   28  177  134  132
# PMC             2815  394 1446  965 3585 2192 1734
# SM cell           11    9   69    8   18   61   17
# T cell           275  645 1414  284  470 2047 1155

n_cells <- as.numeric(table(gc.labeled.sce$sample))
m <- match(sids, gc.labeled.sce$sample)
ei <- data.frame(colData(gc.labeled.sce)[m, ],
                  n_cells, row.names = NULL) %>%
                select(-"neighbor_cell.type_hvg10")

groups <- colData(gc.labeled.sce)[, c("neighbor_cell.type_hvg10", "sample")]
pb <- aggregate.Matrix(t(counts(gc.labeled.sce)),
                       groupings = groups, fun = "sum")
splitf <- sapply(stringr::str_split(rownames(pb),pattern = "_",n = 2),
                 `[`, 1)
pb <- split.data.frame(pb,factor(splitf)) %>%
        lapply(function(u)set_colnames(t(u),
               stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))



get_sample_ids <- function(x){pb[[x]] %>%colnames()}
de_samples <- map(1:length(kids), get_sample_ids) %>% unlist()
samples_list <- map(1:length(kids), get_sample_ids)
get_cluster_ids <- function(x){rep(names(pb)[x],each = length(samples_list[[x]]))}
de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>% unlist()



setwd("deseq2_plots")
interest <- data.frame(ethnicity=levels(as.factor(gc.labeled.sce$ethnicity)),
                       biopsy=levels(as.factor(gc.labeled.sce$biopsy)),
                       disease.status=levels(as.factor(gc.labeled.sce$disease.status)))
clusters <- levels(as.factor(metadata$neighbor_cell.type_hvg10))
padj_cutoff <- 0.05










#ethnicity
gg_df <- data.frame(neighbor_cell.type_hvg10 = de_cluster_ids,
                    sample = de_samples)
gg_df <- left_join(gg_df, ei[, c("sample", "ethnicity")])
metadata <- gg_df %>%
        dplyr::select(neighbor_cell.type_hvg10, sample, ethnicity)

#2,6,8,9,10 work
for(j in 1:12) {
  tryCatch({
  print(j)

  cluster_metadata <- metadata[which(metadata$neighbor_cell.type_hvg10==clusters[j]), ]
  rownames(cluster_metadata) <- cluster_metadata$sample
  counts <- pb[[clusters[j]]]
  cluster_counts <- as.data.frame(counts[,
                                    which(colnames(counts) %in% rownames(cluster_metadata))])

  dds <- DESeqDataSetFromMatrix(cluster_counts,
                                colData = cluster_metadata,
                                design = ~ ethnicity)

  #QC
  rld <- rlog(dds, blind=TRUE)
  pdf(paste0("deseq2 ","ethnicity ",clusters[j]," pca plot.pdf"))
  print(DESeq2::plotPCA(rld, intgroup = "ethnicity"))
  dev.off()

  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  pdf(paste0("deseq2 ","ethnicity ",clusters[j]," heatmap.pdf"))
  print(pheatmap(rld_cor, annotation = cluster_metadata[, c("ethnicity"), drop=F]))
  dev.off()

  dds <- DESeq(dds)
  dds <- estimateSizeFactors(dds)
  colnames(dds) <- as.character(1:7)
  pdf(paste0("deseq2 ","ethnicity ",clusters[j]," dispersion plot.pdf"))
  print(plotDispEsts(dds))
  dev.off()

  #results
  contrast <- c("ethnicity", interest[,"ethnicity"])
  res <- results(dds, contrast = contrast, alpha = 0.05)
  res <- lfcShrink(dds, contrast =  contrast, res=res, type="ashr")

  res_tbl <- res %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
  write.csv(res_tbl,
            paste0("results_", clusters[j], "_", interest[1,"ethnicity"], "_vs_", interest[2,"ethnicity"], "_all_genes.csv"),
            quote = FALSE,
            row.names = FALSE)

  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>% dplyr::arrange(padj)
  write.csv(sig_res,
            paste0("results_", clusters[j], "_", interest[1,"ethnicity"], "_vs_", interest[2,"ethnicity"], "_sig_genes.csv"),
            quote = FALSE,
            row.names = FALSE)

  #plotting
  normalized_counts <- counts(dds, normalized = TRUE)
  log_nor_counts <- log(normalized_counts+1)
  top20_sig_genes <- sig_res %>% dplyr::arrange(padj) %>% dplyr::pull(gene) %>% head(n=20)
  top20_sig_norm <- as.data.frame(log_nor_counts) %>% rownames_to_column(var = "gene") %>%
          dplyr::filter(gene %in% top20_sig_genes)
  gathered_top20_sig <- top20_sig_norm %>%
          gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "log_nor_counts")
  gathered_top20_sig <- inner_join(ei[, c("sample", "ethnicity" )], gathered_top20_sig, by = c("sample" = "samplename"))

  pdf(paste0("deseq2 ","ethnicity ",clusters[j]," top 20 genes.pdf"))
  print(ggplot(gathered_top20_sig) +
          geom_point(aes(x = gene,
                         y = log_nor_counts,
                         color = ethnicity),
                     position=position_jitter(w=0.1,h=0)) +
          xlab("Genes") +
          ylab("log Normalized Counts") +
          ggtitle(paste0("ethnicity ",clusters[j]," Top 20 Significant DE Genes")) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(plot.title = element_text(hjust = 0.5)))
  dev.off()

  sig_norm <- as.data.frame(log_nor_counts) %>% rownames_to_column(var = "gene") %>%
          dplyr::filter(gene %in% sig_res$gene)

  pdf(paste0("deseq2 ","ethnicity ",clusters[j]," heatmap top 20 genes.pdf"))
  print(pheatmap(sig_norm[ , 2:length(colnames(sig_norm))],
      # color = heat_colors,
      cluster_rows = T,
      show_rownames = F,
      annotation = cluster_metadata[, c("ethnicity", "neighbor_cell.type_hvg10")],
      border_color = NA,
      # fontsize = 10,
      # scale = "row",
      # fontsize_row = 10,
      # height = 20
    ))
  dev.off()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



#biopsy
gg_df <- data.frame(neighbor_cell.type_hvg10 = de_cluster_ids,
                    sample = de_samples)
gg_df <- left_join(gg_df, ei[, c("sample", "biopsy")])
metadata <- gg_df %>%
        dplyr::select(neighbor_cell.type_hvg10, sample, biopsy)

#
for(j in 1:12) {
tryCatch({
  print(j)

  cluster_metadata <- metadata[which(metadata$neighbor_cell.type_hvg10==clusters[j]), ]
  rownames(cluster_metadata) <- cluster_metadata$sample
  counts <- pb[[clusters[j]]]
  cluster_counts <- as.data.frame(counts[,
                                    which(colnames(counts) %in% rownames(cluster_metadata))])

  dds <- DESeqDataSetFromMatrix(cluster_counts,
                                colData = cluster_metadata,
                                design = ~ biopsy)

  #QC
  rld <- rlog(dds, blind=TRUE)
  pdf(paste0("deseq2 ","biopsy ",clusters[j]," pca plot.pdf"))
  print(DESeq2::plotPCA(rld, intgroup = "biopsy"))
  dev.off()

  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  pdf(paste0("deseq2 ","biopsy ",clusters[j]," heatmap.pdf"))
  print(pheatmap(rld_cor, annotation = cluster_metadata[, c("biopsy"), drop=F]))
  dev.off()

  dds <- DESeq(dds)
  dds <- estimateSizeFactors(dds)
  colnames(dds) <- as.character(1:7)
  pdf(paste0("deseq2 ","biopsy ",clusters[j]," dispersion plot.pdf"))
  print(plotDispEsts(dds))
  dev.off()

  #results
  contrast <- c("biopsy", interest[,"biopsy"])
  res <- results(dds, contrast = contrast, alpha = 0.05)
  res <- lfcShrink(dds, contrast =  contrast, res=res, type="ashr")

  res_tbl <- res %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
  write.csv(res_tbl,
            paste0("results_", clusters[j], "_", interest[1,"biopsy"], "_vs_", interest[2,"biopsy"], "_all_genes.csv"),
            quote = FALSE,
            row.names = FALSE)

  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>% dplyr::arrange(padj)
  write.csv(sig_res,
            paste0("results_", clusters[j], "_", interest[1,"biopsy"], "_vs_", interest[2,"biopsy"], "_sig_genes.csv"),
            quote = FALSE,
            row.names = FALSE)

  #plotting
  normalized_counts <- counts(dds, normalized = TRUE)
  log_nor_counts <- log(normalized_counts+1)
  top20_sig_genes <- sig_res %>% dplyr::arrange(padj) %>% dplyr::pull(gene) %>% head(n=20)
  top20_sig_norm <- as.data.frame(log_nor_counts) %>% rownames_to_column(var = "gene") %>%
          dplyr::filter(gene %in% top20_sig_genes)
  gathered_top20_sig <- top20_sig_norm %>%
          gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "log_nor_counts")
  gathered_top20_sig <- inner_join(ei[, c("sample", "biopsy" )], gathered_top20_sig, by = c("sample" = "samplename"))

  pdf(paste0("deseq2 ","biopsy ",clusters[j]," top 20 genes.pdf"))
  print(ggplot(gathered_top20_sig) +
          geom_point(aes(x = gene,
                         y = log_nor_counts,
                         color = biopsy),
                     position=position_jitter(w=0.1,h=0)) +
          xlab("Genes") +
          ylab("log Normalized Counts") +
          ggtitle(paste0("biopsy ",clusters[j]," Top 20 Significant DE Genes")) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(plot.title = element_text(hjust = 0.5)))
  dev.off()

  sig_norm <- as.data.frame(log_nor_counts) %>% rownames_to_column(var = "gene") %>%
          dplyr::filter(gene %in% sig_res$gene)

  pdf(paste0("deseq2 ","biopsy ",clusters[j]," heatmap top 20 genes.pdf"))
  print(pheatmap(sig_norm[ , 2:length(colnames(sig_norm))],
      # color = heat_colors,
      cluster_rows = T,
      show_rownames = F,
      annotation = cluster_metadata[, c("biopsy", "neighbor_cell.type_hvg10")],
      border_color = NA,
      # fontsize = 10,
      # scale = "row",
      # fontsize_row = 10,
      # height = 20
    ))
  dev.off()
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


#disease.status
gg_df <- data.frame(neighbor_cell.type_hvg10 = de_cluster_ids,
                    sample = de_samples)
gg_df <- left_join(gg_df, ei[, c("sample", "disease.status")])
metadata <- gg_df %>%
        dplyr::select(neighbor_cell.type_hvg10, sample, disease.status)

#
for(j in 1:12) {
tryCatch({
  print(j)

  cluster_metadata <- metadata[which(metadata$neighbor_cell.type_hvg10==clusters[j]), ]
  rownames(cluster_metadata) <- cluster_metadata$sample
  counts <- pb[[clusters[j]]]
  cluster_counts <- as.data.frame(counts[,
                                    which(colnames(counts) %in% rownames(cluster_metadata))])

  dds <- DESeqDataSetFromMatrix(cluster_counts,
                                colData = cluster_metadata,
                                design = ~ disease.status)

  #QC
  rld <- rlog(dds, blind=TRUE)
  pdf(paste0("deseq2 ","disease.status ",clusters[j]," pca plot.pdf"))
  print(DESeq2::plotPCA(rld, intgroup = "disease.status"))
  dev.off()

  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  pdf(paste0("deseq2 ","disease.status ",clusters[j]," heatmap.pdf"))
  print(pheatmap(rld_cor, annotation = cluster_metadata[, c("disease.status"), drop=F]))
  dev.off()

  dds <- DESeq(dds)
  dds <- estimateSizeFactors(dds)
  colnames(dds) <- as.character(1:7)
  pdf(paste0("deseq2 ","disease.status ",clusters[j]," dispersion plot.pdf"))
  print(plotDispEsts(dds))
  dev.off()

  #results
  contrast <- c("disease.status", interest[,"disease.status"])
  res <- results(dds, contrast = contrast, alpha = 0.05)
  res <- lfcShrink(dds, contrast =  contrast, res=res, type="ashr")

  res_tbl <- res %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
  write.csv(res_tbl,
            paste0("results_", clusters[j], "_", interest[1,"disease.status"], "_vs_", interest[2,"disease.status"], "_all_genes.csv"),
            quote = FALSE,
            row.names = FALSE)

  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>% dplyr::arrange(padj)
  write.csv(sig_res,
            paste0("results_", clusters[j], "_", interest[1,"disease.status"], "_vs_", interest[2,"disease.status"], "_sig_genes.csv"),
            quote = FALSE,
            row.names = FALSE)

  #plotting
  normalized_counts <- counts(dds, normalized = TRUE)
  log_nor_counts <- log(normalized_counts+1)
  top20_sig_genes <- sig_res %>% dplyr::arrange(padj) %>% dplyr::pull(gene) %>% head(n=20)
  top20_sig_norm <- as.data.frame(log_nor_counts) %>% rownames_to_column(var = "gene") %>%
          dplyr::filter(gene %in% top20_sig_genes)
  gathered_top20_sig <- top20_sig_norm %>%
          gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "log_nor_counts")
  gathered_top20_sig <- inner_join(ei[, c("sample", "disease.status" )], gathered_top20_sig, by = c("sample" = "samplename"))

  pdf(paste0("deseq2 ","disease.status ",clusters[j]," top 20 genes.pdf"))
  print(ggplot(gathered_top20_sig) +
          geom_point(aes(x = gene,
                         y = log_nor_counts,
                         color = disease.status),
                     position=position_jitter(w=0.1,h=0)) +
          scale_y_log10() +
          xlab("Genes") +
          ylab("log10 Normalized Counts") +
          ggtitle(paste0("disease.status ",clusters[j]," Top 20 Significant DE Genes")) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(plot.title = element_text(hjust = 0.5)))
  dev.off()

  sig_norm <- as.data.frame(log_nor_counts) %>% rownames_to_column(var = "gene") %>%
          dplyr::filter(gene %in% sig_res$gene)

  pdf(paste0("deseq2 ","disease.status ",clusters[j]," heatmap top 20 genes.pdf"))
  print(pheatmap(sig_norm[ , 2:length(colnames(sig_norm))],
      # color = heat_colors,
      cluster_rows = T,
      show_rownames = F,
      annotation = cluster_metadata[, c("disease.status", "neighbor_cell.type_hvg10")],
      border_color = NA,
      # fontsize = 10,
      # scale = "row",
      # fontsize_row = 10,
      # height = 20
    ))
  dev.off()
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}











gg_df <- data.frame(neighbor_cell.type_hvg10 = de_cluster_ids,
                    sample = de_samples)
gg_df <- left_join(gg_df, ei[, c("sample", "disease.status")])
metadata <- gg_df %>%
        dplyr::select(neighbor_cell.type_hvg10, sample, disease.status)


cluster_metadata <- metadata[which(metadata$neighbor_cell.type_hvg10==clusters[1]), ]
rownames(cluster_metadata) <- cluster_metadata$sample
counts <- pb[[clusters[1]]]
cluster_counts <- as.data.frame(counts[,
                                  which(colnames(counts) %in% rownames(cluster_metadata))])
#####data.frame => as.data.frame

dds <- DESeqDataSetFromMatrix(cluster_counts,
                              colData = cluster_metadata,
                              design = ~ disease.status)


#QC
rld <- rlog(dds, blind=TRUE)
pdf("deseq2 pca plot.pdf")
DESeq2::plotPCA(rld, intgroup = "disease.status")
dev.off()

rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pdf("deseq2 heatmap.pdf")
pheatmap(rld_cor, annotation = cluster_metadata[, c("disease.status"), drop=F])
dev.off()

dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
# colnames(dds) <- as.character(1:7)
pdf("deseq2 dispersion plot.pdf")
plotDispEsts(dds)
dev.off()

#results
levels(as.factor(cluster_metadata$disease.status))[2]
levels(as.factor(cluster_metadata$disease.status))[1]

contrast <- c("disease.status", levels(as.factor(cluster_metadata$disease.status))[2], levels(as.factor(cluster_metadata$disease.status))[1])

res <- results(dds, contrast = contrast, alpha = 0.05)
res <- lfcShrink(dds, contrast =  contrast, res=res, type="ashr")

res_tbl <- res %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
head(res_tbl)
write.csv(res_tbl,
          paste0("results", clusters[1], "_", levels(as.factor(cluster_metadata$sample))[2], "_vs_", levels(as.factor(cluster_metadata$sample))[1], "_all_genes.csv"),
          quote = FALSE,
          row.names = FALSE)

# padj_cutoff <- 0.5
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>% dplyr::arrange(padj)
head(sig_res)
write.csv(sig_res,
          paste0("results", clusters[1], "_", levels(as.factor(cluster_metadata$sample))[2], "_vs_", levels(as.factor(cluster_metadata$sample))[1], "_sig_genes.csv"),
          quote = FALSE,
          row.names = FALSE)

normalized_counts <- counts(dds, normalized = TRUE)
log_nor_counts <- log(normalized_counts+1)
top20_sig_genes <- sig_res %>% dplyr::arrange(padj) %>% dplyr::pull(gene) %>% head(n=20)
top20_sig_norm <- as.data.frame(log_nor_counts) %>% rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% top20_sig_genes)
gathered_top20_sig <- top20_sig_norm %>%
        gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "log_nor_counts")
gathered_top20_sig <- inner_join(ei[, c("sample", "disease.status" )], gathered_top20_sig, by = c("sample" = "samplename"))

pdf("deseq2 top 20 genes.pdf")
ggplot(gathered_top20_sig) +
        geom_point(aes(x = gene,
                       y = log_nor_counts,
                       color = disease.status),
                   position=position_jitter(w=0.1,h=0)) +
        xlab("Genes") +
        ylab("log Normalized Counts") +
        ggtitle("Top 20 Significant DE Genes") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()

sig_norm <- as.data.frame(log_nor_counts) %>% rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% sig_res$gene)
# colnames(sig_norm)[2:length(colnames(sig_norm))] <- as.character(1:7)
# heat_colors <- brewer.pal(6, "YlOrRd")

pdf("deseq2 heatmap top 20 genes.pdf")
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))],
    # color = heat_colors,
    cluster_rows = T,
    show_rownames = F,
    annotation = cluster_metadata[, c("disease.status", "neighbor_cell.type_hvg10")],
    border_color = NA,
    # fontsize = 10,
    # scale = "row",
    # fontsize_row = 10,
    # height = 20
  )
dev.off()
