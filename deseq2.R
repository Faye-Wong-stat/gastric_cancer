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
gc.labeled.sce <- as.SingleCellExperiment(gc.labeled)
gc.labeled.sce.meta <- colData(gc.labeled.sce)[,
c("sample","disease.status","ethnicity","biopsy","subject")]
gc.labeled.sce.meta.sample <- plyr::count(gc.labeled.sce.meta)



groups <- colData(gc.labeled.sce)[, c("neighbor_cell.type_hvg10", "sample")]
pb <- aggregate.Matrix(t(counts(gc.labeled.sce)),
                       groupings = groups, fun = "sum")
splitf <- sapply(stringr::str_split(rownames(pb),
                                    pattern = "_",
                                    n = 2),
                 `[`, 1)
pb <- split.data.frame(pb,
                       factor(splitf)) %>%
        lapply(function(u)
                set_colnames(t(u),
                             stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))



kids <- purrr::set_names(levels(as.factor(gc.labeled.sce$neighbor_cell.type_hvg10)))
nk <- length(kids)
sids <- purrr::set_names(levels(as.factor(gc.labeled.sce$sample)))
ns <- length(sids)

n_cells <- as.numeric(table(gc.labeled.sce$sample))
m <- match(sids, gc.labeled.sce$sample)
ei <- data.frame(colData(gc.labeled.sce)[m, ],
                  n_cells, row.names = NULL) %>%
                select(-"neighbor_cell.type_hvg10")

get_sample_ids <- function(x){
        pb[[x]] %>%
                colnames()
}
de_samples <- map(1:length(kids), get_sample_ids) %>%
        unlist()
samples_list <- map(1:length(kids), get_sample_ids)
get_cluster_ids <- function(x){
        rep(names(pb)[x],
            each = length(samples_list[[x]]))
}
de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
        unlist()



gg_df <- data.frame(neighbor_cell.type_hvg10 = de_cluster_ids,
                    sample = de_samples)
gg_df <- left_join(gg_df, ei[, c("sample", "ethnicity")])
metadata <- gg_df %>%
        dplyr::select(neighbor_cell.type_hvg10, sample, ethnicity)



clusters <- levels(as.factor(metadata$neighbor_cell.type_hvg10))
cluster_metadata <- metadata[which(metadata$neighbor_cell.type_hvg10==clusters[2]), ]
rownames(cluster_metadata) <- cluster_metadata$sample
counts <- pb[[clusters[2]]]
cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

dds <- DESeqDataSetFromMatrix(cluster_counts,
                              colData = cluster_metadata,
                              design = ~ ethnicity)


#QC
rld <- rlog(dds, blind=TRUE)
pdf("deseq2 pca plot.pdf")
DESeq2::plotPCA(rld, intgroup = "ethnicity")
dev.off()

rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pdf("deseq2 heatmap.pdf")
pheatmap(rld_cor, annotation = cluster_metadata[, c("ethnicity"), drop=F])
dev.off()

dds <- DESeq(dds)
colnames(dds) <- as.character(1:7)
pdf("deseq2 dispersion plot.pdf")
plotDispEsts(dds)
dev.off()

#results
levels(as.factor(cluster_metadata$ethnicity))[2]
levels(as.factor(cluster_metadata$ethnicity))[1]

contrast <- c("ethnicity", levels(as.factor(cluster_metadata$ethnicity))[2], levels(as.factor(cluster_metadata$ethnicity))[1])

res <- results(dds,
               contrast = contrast,
               alpha = 0.05)
res <- lfcShrink(dds,
                 contrast =  contrast,
                 res=res, type="ashr")

res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
head(res_tbl)
write.csv(res_tbl,
          paste0("results", clusters[2], "_", levels(as.factor(cluster_metadata$sample))[2], "_vs_", levels(as.factor(cluster_metadata$sample))[1], "_all_genes.csv"),
          quote = FALSE,
          row.names = FALSE)

padj_cutoff <- 0.9
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
        dplyr::arrange(padj)
head(sig_res)
write.csv(sig_res,
          paste0("results", clusters[2], "_", levels(as.factor(cluster_metadata$sample))[2], "_vs_", levels(as.factor(cluster_metadata$sample))[1], "_sig_genes.csv"),
          quote = FALSE,
          row.names = FALSE)



normalized_counts <- counts(dds,
                            normalized = TRUE)

top20_sig_genes <- sig_res %>%
        dplyr::arrange(padj) %>%
        dplyr::pull(gene) %>%
        head(n=20)

top20_sig_norm <- as.data.frame(normalized_counts) %>%
        rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% top20_sig_genes)

gathered_top20_sig <- top20_sig_norm %>%
        gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

gathered_top20_sig <- inner_join(ei[, c("sample", "ethnicity" )], gathered_top20_sig, by = c("sample" = "samplename"))

pdf("deseq2 top 20 genes.pdf")
ggplot(gathered_top20_sig) +
        geom_point(aes(x = gene,
                       y = normalized_counts,
                       color = ethnicity),
                   position=position_jitter(w=0.1,h=0)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("log10 Normalized Counts") +
        ggtitle("Top 20 Significant DE Genes") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5))
dev.off()
