#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript de_diagnostics_sensitivity.R ",
    "<counts.tsv> <metadata.tsv> <limma_results.tsv> <deseq2_results.tsv> <outdir>\n",
    "Expected metadata columns: sample_id, patient_id, condition"
  )
}

counts_file <- args[[1]]
metadata_file <- args[[2]]
limma_file <- args[[3]]
deseq2_file <- args[[4]]
outdir <- args[[5]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

pick_col <- function(df, candidates, required = TRUE) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  if (required) {
    stop("None of these columns were found: ", paste(candidates, collapse = ", "))
  }
  NA_character_
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}

counts <- read_table(counts_file)
metadata <- read_table(metadata_file)
limma <- read_table(limma_file)
deseq2 <- read_table(deseq2_file)

gene_col_counts <- colnames(counts)[[1]]
rownames(counts) <- counts[[gene_col_counts]]
count_matrix <- as.matrix(counts[, setdiff(colnames(counts), gene_col_counts), drop = FALSE])
mode(count_matrix) <- "numeric"

sample_id_col <- pick_col(metadata, c("sample_id", "sample", "Sample", "sample_name"))
patient_id_col <- pick_col(metadata, c("patient_id", "patient", "Patient"))
condition_col <- pick_col(metadata, c("condition", "status", "Condition"))

common_samples <- intersect(colnames(count_matrix), metadata[[sample_id_col]])
if (length(common_samples) < 2) {
  stop("Fewer than two samples overlap between count matrix and metadata.")
}
count_matrix <- count_matrix[, common_samples, drop = FALSE]
metadata <- metadata[match(common_samples, metadata[[sample_id_col]]), , drop = FALSE]

condition <- metadata[[condition_col]]
patient_id <- metadata[[patient_id_col]]

library_size <- colSums(count_matrix, na.rm = TRUE)
detected_genes <- colSums(count_matrix > 0, na.rm = TRUE)
sample_qc <- data.frame(
  sample_id = common_samples,
  patient_id = patient_id,
  condition = condition,
  library_size = library_size,
  detected_genes = detected_genes,
  stringsAsFactors = FALSE
)
write_tsv(sample_qc, file.path(outdir, "de_diagnostics_sample_qc.tsv"))

pairing <- aggregate(
  sample_qc$sample_id,
  by = list(patient_id = sample_qc$patient_id, condition = sample_qc$condition),
  FUN = length
)
colnames(pairing)[[3]] <- "n_samples"
write_tsv(pairing, file.path(outdir, "de_diagnostics_pairing_summary.tsv"))

lib <- pmax(colSums(count_matrix, na.rm = TRUE), 1)
cpm <- sweep(count_matrix, 2, lib, "/") * 1e6
log_cpm <- log2(cpm + 1)

top_var_n <- min(500, nrow(log_cpm))
vars <- apply(log_cpm, 1, var, na.rm = TRUE)
top_genes <- names(sort(vars, decreasing = TRUE))[seq_len(top_var_n)]
pca <- prcomp(t(log_cpm[top_genes, , drop = FALSE]), scale. = TRUE)
pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  patient_id = patient_id,
  condition = condition,
  stringsAsFactors = FALSE
)
write_tsv(pca_df, file.path(outdir, "de_diagnostics_pca_coordinates.tsv"))

pdf(file.path(outdir, "de_diagnostics_pca.pdf"), width = 7, height = 6)
plot(
  pca_df$PC1, pca_df$PC2,
  pch = ifelse(tolower(pca_df$condition) %in% c("tumour", "tumor", "t"), 19, 1),
  col = as.factor(pca_df$condition),
  xlab = sprintf("PC1 (%.1f%%)", 100 * summary(pca)$importance[2, 1]),
  ylab = sprintf("PC2 (%.1f%%)", 100 * summary(pca)$importance[2, 2]),
  main = "RNA-seq logCPM PCA"
)
legend("topright", legend = unique(pca_df$condition), col = seq_along(unique(pca_df$condition)), pch = 19, bty = "n")
dev.off()

gene_col_limma <- pick_col(limma, c("gene_id", "gene", "Gene", "genes", "ID", colnames(limma)[[1]]))
gene_col_deseq2 <- pick_col(deseq2, c("gene_id", "gene", "Gene", "genes", "ID", colnames(deseq2)[[1]]))
lfc_limma_col <- pick_col(limma, c("logFC", "log2FoldChange", "log2FC", "coef"))
lfc_deseq_col <- pick_col(deseq2, c("DESeq2_log2FoldChange", "log2FoldChange", "logFC", "log2FC"))
padj_limma_col <- pick_col(limma, c("adj.P.Val", "FDR", "padj", "qvalue"), required = FALSE)
padj_deseq_col <- pick_col(deseq2, c("DESeq2_padj", "padj", "adj.P.Val", "FDR", "qvalue"), required = FALSE)

limma_small <- data.frame(
  gene_id = limma[[gene_col_limma]],
  limma_logFC = as.numeric(limma[[lfc_limma_col]]),
  limma_padj = if (!is.na(padj_limma_col)) as.numeric(limma[[padj_limma_col]]) else NA_real_,
  stringsAsFactors = FALSE
)
deseq_small <- data.frame(
  gene_id = deseq2[[gene_col_deseq2]],
  deseq2_log2FoldChange = as.numeric(deseq2[[lfc_deseq_col]]),
  deseq2_padj = if (!is.na(padj_deseq_col)) as.numeric(deseq2[[padj_deseq_col]]) else NA_real_,
  stringsAsFactors = FALSE
)
merged <- merge(limma_small, deseq_small, by = "gene_id")
merged$abs_lfc_difference <- abs(merged$limma_logFC - merged$deseq2_log2FoldChange)
merged$limma_significant <- !is.na(merged$limma_padj) & merged$limma_padj < 0.05 & abs(merged$limma_logFC) > 1
merged$deseq2_significant <- !is.na(merged$deseq2_padj) & merged$deseq2_padj < 0.05 & abs(merged$deseq2_log2FoldChange) > 1
write_tsv(merged, file.path(outdir, "de_diagnostics_limma_deseq2_merged.tsv"))

diagnostics <- data.frame(
  metric = c(
    "n_samples",
    "n_patients",
    "n_genes_in_counts",
    "n_genes_compared_between_limma_and_deseq2",
    "limma_significant_FDR_0.05_absLFC_1",
    "deseq2_significant_FDR_0.05_absLFC_1",
    "overlap_significant_both",
    "logFC_pearson_correlation",
    "logFC_spearman_correlation",
    "median_abs_logFC_difference"
  ),
  value = c(
    length(common_samples),
    length(unique(patient_id)),
    nrow(count_matrix),
    nrow(merged),
    sum(merged$limma_significant, na.rm = TRUE),
    sum(merged$deseq2_significant, na.rm = TRUE),
    sum(merged$limma_significant & merged$deseq2_significant, na.rm = TRUE),
    suppressWarnings(cor(merged$limma_logFC, merged$deseq2_log2FoldChange, use = "complete.obs", method = "pearson")),
    suppressWarnings(cor(merged$limma_logFC, merged$deseq2_log2FoldChange, use = "complete.obs", method = "spearman")),
    median(merged$abs_lfc_difference, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write_tsv(diagnostics, file.path(outdir, "de_diagnostics_summary.tsv"))

discordant <- merged[order(merged$abs_lfc_difference, decreasing = TRUE), , drop = FALSE]
write_tsv(head(discordant, 100), file.path(outdir, "de_diagnostics_top_discordant_genes.tsv"))

pdf(file.path(outdir, "de_diagnostics_limma_vs_deseq2_logfc.pdf"), width = 7, height = 7)
plot(
  merged$limma_logFC, merged$deseq2_log2FoldChange,
  pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.25),
  xlab = "limma-voom logFC",
  ylab = "DESeq2 log2FoldChange",
  main = "Tumour-normal effect-size concordance"
)
abline(h = 0, v = 0, col = "grey70")
abline(0, 1, col = "red", lwd = 2)
dev.off()

pdf(file.path(outdir, "de_diagnostics_library_sizes.pdf"), width = 8, height = 5)
barplot(
  sample_qc$library_size / 1e6,
  names.arg = sample_qc$sample_id,
  las = 2,
  col = as.factor(sample_qc$condition),
  ylab = "Library size (millions of assigned reads)",
  main = "RNA-seq count library sizes"
)
dev.off()
