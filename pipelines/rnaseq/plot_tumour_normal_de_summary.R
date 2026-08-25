#!/usr/bin/env Rscript

# Manuscript-ready summary of a paired tumour-normal RNA-seq comparison.
# The plot makes method sensitivity visible rather than presenting DESeq2-only
# discoveries as cross-method consensus findings.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: plot_tumour_normal_de_summary.R <normalized_counts.tsv> <metadata.tsv> <deseq2.tsv> <edger.tsv> <limma.tsv> <outdir>")
}

normalized_file <- args[[1]]
metadata_file <- args[[2]]
deseq_file <- args[[3]]
edge_file <- args[[4]]
limma_file <- args[[5]]
outdir <- args[[6]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

normalized <- fread(normalized_file, check.names = FALSE)
metadata <- fread(metadata_file, colClasses = "character")
deseq <- fread(deseq_file)
edge <- fread(edge_file)
limma <- fread(limma_file)

required_metadata <- c("sample_id", "patient_id", "condition")
if (!all(required_metadata %in% names(metadata))) stop("Metadata must contain sample_id, patient_id and condition")

gene_column <- names(normalized)[1]
genes <- normalized[[gene_column]]
normalized[[gene_column]] <- NULL
sample_order <- intersect(metadata$sample_id, names(normalized))
if (length(sample_order) != nrow(metadata)) stop("Normalized counts and metadata sample sets do not match")
matrix <- as.matrix(normalized[, ..sample_order])
storage.mode(matrix) <- "numeric"
rownames(matrix) <- genes
log_matrix <- log2(matrix + 1)
variances <- apply(log_matrix, 1, var)
top_n <- min(500L, length(variances))
top_genes <- names(sort(variances, decreasing = TRUE))[seq_len(top_n)]
pca <- prcomp(t(log_matrix[top_genes, , drop = FALSE]), scale. = TRUE)
variance_explained <- 100 * summary(pca)$importance[2, 1:2]

pca_dt <- data.table(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)
pca_dt <- merge(pca_dt, metadata[, ..required_metadata], by = "sample_id", all.x = TRUE)
p_pca <- ggplot(pca_dt, aes(PC1, PC2)) +
  geom_line(aes(group = patient_id), colour = "#B8B8B8", linewidth = 0.45) +
  geom_point(aes(fill = condition, shape = condition), size = 2.8, colour = "#303030") +
  scale_fill_manual(values = c("Normal" = "#D8D1C3", "Tumour" = "#B44C43")) +
  scale_shape_manual(values = c("Normal" = 21, "Tumour" = 22)) +
  labs(
    title = "A  Paired RNA-seq PCA",
    subtitle = "Lines connect matched tumour-normal samples",
    x = sprintf("PC1 (%.1f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", variance_explained[2]),
    fill = NULL, shape = NULL
  ) +
  theme_minimal(base_family = "serif", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
    legend.position = c(0.84, 0.14),
    legend.background = element_rect(fill = scales::alpha("white", 0.82), colour = NA)
  )

merged <- merge(
  deseq[, .(gene_id, gene_symbol, baseMean, DESeq2_log2FC = DESeq2_log2FoldChange, DESeq2_P = DESeq2_pvalue, DESeq2_FDR = DESeq2_padj)],
  edge[, .(gene_id, edgeR_log2FC = logFC, edgeR_P = PValue, edgeR_FDR = FDR)],
  by = "gene_id"
)
merged <- merge(
  merged,
  limma[, .(gene_id, limma_log2FC = logFC, limma_P = P.Value, limma_FDR = adj.P.Val)],
  by = "gene_id"
)
merged[, DESeq2_significant := !is.na(DESeq2_FDR) & DESeq2_FDR < 0.05 & abs(DESeq2_log2FC) > 1]
merged[, sensitivity_nominal_n := rowSums(cbind(edgeR_P < 0.05, limma_P < 0.05), na.rm = TRUE)]
merged[, support_category := fifelse(
  DESeq2_significant & sensitivity_nominal_n == 2,
    "DESeq2 + both nominal",
  fifelse(
    DESeq2_significant,
    "DESeq2 + incomplete",
    "Below threshold"
  )
)]
merged[, minus_log10_fdr := -log10(pmax(DESeq2_FDR, .Machine$double.xmin))]

labels <- merged[
  DESeq2_significant == TRUE & sensitivity_nominal_n == 2 & !is.na(gene_symbol) & gene_symbol != ""
][order(DESeq2_P)][seq_len(min(10L, .N))]

category_colours <- c(
  "Below threshold" = "#BDBDBD",
  "DESeq2 + incomplete" = "#D79A43",
  "DESeq2 + both nominal" = "#B44C43"
)
p_volcano <- ggplot(merged, aes(DESeq2_log2FC, minus_log10_fdr, colour = support_category)) +
  geom_point(size = 0.7, alpha = 0.55) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "#777777") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#777777") +
  geom_text(data = labels, aes(label = gene_symbol), size = 2.5, check_overlap = TRUE, vjust = -0.6, show.legend = FALSE) +
  scale_colour_manual(values = category_colours) +
  labs(
    title = "B  DESeq2 tumour-normal contrast",
    subtitle = "No gene met the joint DESeq2-edgeR FDR rule",
    x = "DESeq2 log2 fold change", y = "-log10(DESeq2 FDR)", colour = "Statistical support"
  ) +
  theme_minimal(base_family = "serif", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
    legend.position = c(0.19, 0.84),
    legend.direction = "vertical",
    legend.background = element_rect(fill = scales::alpha("white", 0.86), colour = NA),
    legend.text = element_text(size = 7.5)
  )

correlation <- cor(merged$DESeq2_log2FC, merged$edgeR_log2FC, use = "complete.obs")
p_concordance <- ggplot(merged, aes(DESeq2_log2FC, edgeR_log2FC, colour = support_category)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, colour = "#303030", linewidth = 0.7) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, label = sprintf("Pearson r = %.3f", correlation), size = 3.1) +
  scale_colour_manual(values = category_colours) +
  coord_equal() +
  labs(
    title = "C  Count-model effect concordance",
    subtitle = "DESeq2 and edgeR quasi-likelihood",
    x = "DESeq2 log2 fold change", y = "edgeR log2 fold change", colour = "Statistical support"
  ) +
  theme_minimal(base_family = "serif", base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(), legend.position = "none")

combined <- p_pca + p_volcano + p_concordance +
  plot_layout(ncol = 3, widths = c(0.95, 1.15, 1.0), guides = "keep")

ggsave(file.path(outdir, "tumour_normal_de_summary.png"), combined, width = 15.5, height = 5.9, dpi = 320)
ggsave(file.path(outdir, "tumour_normal_de_summary.pdf"), combined, width = 15.5, height = 5.9)

summary <- data.table(
  metric = c(
    "samples", "patients", "genes_tested", "DESeq2_FDR_absLFC", "edgeR_FDR_absLFC",
    "limma_FDR_absLFC", "DESeq2_with_both_sensitivity_nominal", "DESeq2_edgeR_logFC_Pearson"
  ),
  value = c(
    nrow(metadata), uniqueN(metadata$patient_id), nrow(merged),
    sum(merged$DESeq2_significant),
    sum(merged$edgeR_FDR < 0.05 & abs(merged$edgeR_log2FC) > 1, na.rm = TRUE),
    sum(merged$limma_FDR < 0.05 & abs(merged$limma_log2FC) > 1, na.rm = TRUE),
    sum(merged$DESeq2_significant & merged$sensitivity_nominal_n == 2),
    correlation
  )
)
fwrite(summary, file.path(outdir, "tumour_normal_de_figure_summary.tsv"), sep = "\t")
writeLines(
  c(
    paste("normalized_counts_md5:", unname(tools::md5sum(normalized_file))),
    paste("metadata_md5:", unname(tools::md5sum(metadata_file))),
    paste("DESeq2_md5:", unname(tools::md5sum(deseq_file))),
    paste("edgeR_md5:", unname(tools::md5sum(edge_file))),
    paste("limma_md5:", unname(tools::md5sum(limma_file))),
    "PCA_input: log2(DESeq2 normalized counts + 1), top 500 variable filtered genes",
    "DE_threshold: FDR < 0.05 and absolute log2 fold change > 1",
    paste("R_version:", R.version.string)
  ),
  file.path(outdir, "tumour_normal_de_figure_method_notes.txt")
)

message("Wrote tumour-normal DE summary figure to: ", normalizePath(outdir))
