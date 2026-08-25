#!/usr/bin/env Rscript

# Compare paired tumour-normal DESeq2, edgeR quasi-likelihood and limma-voom
# results generated from the same filtered count matrix and design.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: de_three_method_diagnostics.R <deseq2.tsv> <edger.tsv> <limma.tsv> <outdir>")
}

deseq_file <- args[[1]]
edge_file <- args[[2]]
limma_file <- args[[3]]
outdir <- args[[4]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

deseq <- fread(deseq_file)
edge <- fread(edge_file)
limma <- fread(limma_file)

required <- list(
  DESeq2 = c("gene_id", "baseMean", "DESeq2_log2FoldChange", "DESeq2_pvalue", "DESeq2_padj"),
  edgeR = c("gene_id", "logFC", "PValue", "FDR"),
  limma = c("gene_id", "logFC", "P.Value", "adj.P.Val")
)
tables <- list(DESeq2 = deseq, edgeR = edge, limma = limma)
for (name in names(required)) {
  table <- tables[[name]]
  missing <- setdiff(required[[name]], names(table))
  if (length(missing) > 0) stop(name, " columns missing: ", paste(missing, collapse = ", "))
}

merged <- merge(
  deseq[, .(
    gene_id, gene_symbol, baseMean,
    DESeq2_log2FC = DESeq2_log2FoldChange,
    DESeq2_P = DESeq2_pvalue,
    DESeq2_FDR = DESeq2_padj
  )],
  edge[, .(
    gene_id,
    edgeR_log2FC = logFC,
    edgeR_P = PValue,
    edgeR_FDR = FDR
  )],
  by = "gene_id"
)
merged <- merge(
  merged,
  limma[, .(
    gene_id,
    limma_log2FC = logFC,
    limma_P = P.Value,
    limma_FDR = adj.P.Val
  )],
  by = "gene_id"
)

merged[, DESeq2_significant := !is.na(DESeq2_FDR) & DESeq2_FDR < 0.05 & abs(DESeq2_log2FC) > 1]
merged[, edgeR_significant := !is.na(edgeR_FDR) & edgeR_FDR < 0.05 & abs(edgeR_log2FC) > 1]
merged[, limma_significant := !is.na(limma_FDR) & limma_FDR < 0.05 & abs(limma_log2FC) > 1]
merged[, direction_concordant :=
  sign(DESeq2_log2FC) == sign(edgeR_log2FC) & sign(DESeq2_log2FC) == sign(limma_log2FC)]
merged[, method_call := fifelse(
  DESeq2_significant & edgeR_significant & limma_significant,
  "all_three",
  fifelse(
    DESeq2_significant & edgeR_significant,
    "DESeq2_edgeR",
    fifelse(DESeq2_significant, "DESeq2_only", "not_significant")
  )
)]

safe_cor <- function(x, y, method = "pearson") {
  suppressWarnings(cor(x, y, use = "complete.obs", method = method))
}

top_overlap <- function(n) {
  top_deseq <- head(merged[order(DESeq2_P)]$gene_id, n)
  top_edge <- head(merged[order(edgeR_P)]$gene_id, n)
  top_limma <- head(merged[order(limma_P)]$gene_id, n)
  length(Reduce(intersect, list(top_deseq, top_edge, top_limma)))
}

deseq_sig <- merged[DESeq2_significant == TRUE]
summary <- data.table(
  metric = c(
    "genes_compared",
    "DESeq2_significant_FDR_0.05_absLFC_1",
    "edgeR_significant_FDR_0.05_absLFC_1",
    "limma_significant_FDR_0.05_absLFC_1",
    "all_three_significant",
    "DESeq2_edgeR_significant",
    "DESeq2_significant_with_edgeR_nominal_P_0.05",
    "DESeq2_significant_with_limma_nominal_P_0.05",
    "DESeq2_significant_direction_concordant_all_three",
    "logFC_pearson_DESeq2_edgeR",
    "logFC_pearson_DESeq2_limma",
    "logFC_pearson_edgeR_limma",
    "logFC_spearman_DESeq2_edgeR",
    "logFC_spearman_DESeq2_limma",
    "top_100_gene_overlap_all_three_by_P",
    "top_500_gene_overlap_all_three_by_P",
    "minimum_FDR_DESeq2",
    "minimum_FDR_edgeR",
    "minimum_FDR_limma",
    "median_baseMean_DESeq2_significant"
  ),
  value = c(
    nrow(merged),
    sum(merged$DESeq2_significant),
    sum(merged$edgeR_significant),
    sum(merged$limma_significant),
    sum(merged$DESeq2_significant & merged$edgeR_significant & merged$limma_significant),
    sum(merged$DESeq2_significant & merged$edgeR_significant),
    sum(deseq_sig$edgeR_P < 0.05, na.rm = TRUE),
    sum(deseq_sig$limma_P < 0.05, na.rm = TRUE),
    sum(deseq_sig$direction_concordant, na.rm = TRUE),
    safe_cor(merged$DESeq2_log2FC, merged$edgeR_log2FC),
    safe_cor(merged$DESeq2_log2FC, merged$limma_log2FC),
    safe_cor(merged$edgeR_log2FC, merged$limma_log2FC),
    safe_cor(merged$DESeq2_log2FC, merged$edgeR_log2FC, "spearman"),
    safe_cor(merged$DESeq2_log2FC, merged$limma_log2FC, "spearman"),
    top_overlap(100),
    top_overlap(500),
    min(merged$DESeq2_FDR, na.rm = TRUE),
    min(merged$edgeR_FDR, na.rm = TRUE),
    min(merged$limma_FDR, na.rm = TRUE),
    median(deseq_sig$baseMean, na.rm = TRUE)
  )
)

fwrite(summary, file.path(outdir, "de_three_method_summary.tsv"), sep = "\t")
fwrite(merged, file.path(outdir, "de_three_method_merged.tsv"), sep = "\t")

png(file.path(outdir, "de_three_method_effect_size_concordance.png"), width = 2400, height = 1100, res = 220)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
plot(
  merged$DESeq2_log2FC, merged$edgeR_log2FC,
  pch = 16, cex = 0.35, col = rgb(0.1, 0.2, 0.25, 0.2),
  xlab = "DESeq2 log2 fold change", ylab = "edgeR QL log2 fold change",
  main = "DESeq2 versus edgeR"
)
abline(0, 1, col = "#B44C43", lwd = 2)
plot(
  merged$DESeq2_log2FC, merged$limma_log2FC,
  pch = 16, cex = 0.35, col = rgb(0.1, 0.2, 0.25, 0.2),
  xlab = "DESeq2 log2 fold change", ylab = "limma-voom log2 fold change",
  main = "DESeq2 versus limma-voom"
)
abline(0, 1, col = "#B44C43", lwd = 2)
dev.off()

writeLines(
  c(
    paste("DESeq2_file:", normalizePath(deseq_file)),
    paste("DESeq2_md5:", unname(tools::md5sum(deseq_file))),
    paste("edgeR_file:", normalizePath(edge_file)),
    paste("edgeR_md5:", unname(tools::md5sum(edge_file))),
    paste("limma_file:", normalizePath(limma_file)),
    paste("limma_md5:", unname(tools::md5sum(limma_file))),
    "threshold: FDR < 0.05 and absolute log2 fold change > 1",
    "interpretation_rule: cross-method significance requires agreement of count-model sensitivity analyses",
    paste("R_version:", R.version.string)
  ),
  file.path(outdir, "de_three_method_method_notes.txt")
)

message("Wrote three-method DE diagnostics to: ", normalizePath(outdir))
