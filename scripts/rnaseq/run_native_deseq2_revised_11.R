#!/usr/bin/env Rscript
# PDAC2026 analysis: Salmon gene counts/TPM, paired DESeq2, response PCA and enrichment.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

source(file.path(Sys.getenv("PDAC_CODE_ROOT"), "scripts/shared/runtime.R"))
root <- pdac_path("project_root")
quant_dir <- pdac_path("salmon_quant_dir")
out_dir <- file.path(root, "results", "rnaseq_revised_11_native_deseq2_20260917")
mapping_path <- file.path(root, "reproducible_pipeline_share_20260713_with_outputs", "raw_outputs", "results", "figures", "categories", "04_cross_cohort_concordance", "rnaseq_ensembl_to_symbol.csv")
patients <- pdac_cohort("primary_pdac")
samples <- as.vector(t(outer(patients, c("N", "T"), paste0)))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
files <- vapply(samples, function(sample) {
  matches <- Sys.glob(file.path(quant_dir, paste0(sample, "RNA_*.quant.genes.sf")))
  if (length(matches) != 1) stop("Expected exactly one gene-level Salmon quantification for ", sample)
  matches
}, character(1))

read_counts <- function(path) {
  # DESeq2 receives estimated counts, never TPM or tumour-normal response scores.
  quant <- fread(path, select = c("Name", "NumReads"))
  setNames(quant$NumReads, quant$Name)
}

counts_list <- lapply(files, read_counts)
genes <- Reduce(union, lapply(counts_list, names))
counts <- do.call(cbind, lapply(counts_list, function(values) values[match(genes, names(values))]))
counts[is.na(counts)] <- 0
colnames(counts) <- samples
rownames(counts) <- sub("\\..*$", "", genes)
counts <- rowsum(counts, group = rownames(counts), reorder = FALSE)
keep <- rowSums(counts >= 10) >= length(patients)
counts <- round(counts[keep, , drop = FALSE])

coldata <- data.frame(
  row.names = samples,
  patient = factor(sub("[NT]$", "", samples)),
  condition = factor(ifelse(grepl("T$", samples), "tumour", "normal"), levels = c("normal", "tumour"))
)
# Patient blocks pair the tissues; normal is the reference for this contrast.
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ patient + condition)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("condition", "tumour", "normal"), alpha = 0.05))
res$gene <- rownames(res)
mapping <- fread(mapping_path)
setnames(mapping, names(mapping)[1:2], c("gene", "gene_symbol"))
res <- merge(res, mapping, by = "gene", all.x = TRUE, sort = FALSE)
res$direction <- ifelse(res$log2FoldChange > 0, "up_in_tumour", "down_in_tumour")
res <- res[, c("gene", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "direction")]
fwrite(res, file.path(out_dir, "native_deseq2_results_11_pdac.csv"))

degs <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 1, ]
fwrite(degs, file.path(out_dir, "native_deseq2_degs_FDR05_abslog2FC1_11_pdac.csv"))
summary <- data.frame(
  metric = c("input_genes", "tested_genes", "significant_FDR05", "significant_FDR05_abslog2FC1", "up_in_tumour_FDR05_abslog2FC1", "down_in_tumour_FDR05_abslog2FC1"),
  value = c(nrow(counts), nrow(res), sum(res$padj < 0.05, na.rm = TRUE), nrow(degs), sum(degs$direction == "up_in_tumour"), sum(degs$direction == "down_in_tumour"))
)
fwrite(summary, file.path(out_dir, "native_deseq2_summary_11_pdac.csv"))

png(file.path(out_dir, "native_deseq2_volcano_11_pdac.png"), width = 2100, height = 1800, res = 300)
plot_data <- res[!is.na(res$padj), ]
plot(plot_data$log2FoldChange, -log10(pmax(plot_data$padj, 1e-300)), pch = 16, cex = 0.35, col = "grey75", xlab = "DESeq2 log2 fold change (tumour / normal)", ylab = "-log10 adjusted p-value", main = "Native paired DESeq2: 11 PDAC patients")
points(degs$log2FoldChange, -log10(pmax(degs$padj, 1e-300)), pch = 16, cex = 0.4, col = ifelse(degs$log2FoldChange > 0, "#d1495b", "#0077b6"))
abline(v = c(-1, 1), h = -log10(0.05), lty = 2, col = "grey40")
dev.off()
