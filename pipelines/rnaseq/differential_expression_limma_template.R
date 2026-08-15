#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: differential_expression_limma_template.R <counts.tsv> <metadata.tsv> <outdir>\n",
    "metadata.tsv must contain sample_id, patient_id and condition columns.\n",
    "condition must contain Tumour and Normal.",
    call. = FALSE
  )
}

counts_file <- args[[1]]
metadata_file <- args[[2]]
outdir <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts_dt <- fread(counts_file)
gene_col <- names(counts_dt)[1]
gene_ids <- counts_dt[[gene_col]]
counts <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(counts) <- gene_ids
storage.mode(counts) <- "numeric"

metadata <- fread(metadata_file)
required_cols <- c("sample_id", "patient_id", "condition")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

metadata <- metadata[condition %in% c("Tumour", "Normal")]
metadata[, condition := factor(condition, levels = c("Normal", "Tumour"))]
metadata[, patient_id := factor(patient_id)]

common_samples <- intersect(colnames(counts), metadata$sample_id)
if (length(common_samples) < 4) {
  stop("Too few overlapping samples between count matrix and metadata.", call. = FALSE)
}
counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[match(common_samples, sample_id)]

keep <- filterByExpr(counts, group = metadata$condition)
counts_filtered <- counts[keep, , drop = FALSE]

dge <- DGEList(counts = counts_filtered)
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ patient_id + condition, data = metadata)
voom_fit <- voom(dge, design, plot = FALSE)
fit <- lmFit(voom_fit, design)
fit <- eBayes(fit)

coef_name <- "conditionTumour"
if (!coef_name %in% colnames(design)) {
  stop("Could not find tumour-vs-normal coefficient: ", coef_name, call. = FALSE)
}

results <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
results <- data.table(gene_id = rownames(results), results)
results[, significant_FDR_0_05_logFC_1 := adj.P.Val < 0.05 & abs(logFC) > 1]

fwrite(
  data.table(gene_id = rownames(counts_filtered), counts_filtered),
  file.path(outdir, "filtered_counts.tsv"),
  sep = "\t"
)
fwrite(
  data.table(gene_id = rownames(voom_fit$E), voom_fit$E),
  file.path(outdir, "voom_logCPM_expression.tsv"),
  sep = "\t"
)
fwrite(
  results,
  file.path(outdir, "DE_tumour_vs_normal_paired_limma.tsv"),
  sep = "\t"
)

summary_dt <- data.table(
  metric = c("input_genes", "filtered_genes", "samples", "tumour_samples", "normal_samples", "significant_genes"),
  value = c(
    nrow(counts),
    nrow(counts_filtered),
    ncol(counts_filtered),
    sum(metadata$condition == "Tumour"),
    sum(metadata$condition == "Normal"),
    sum(results$significant_FDR_0_05_logFC_1, na.rm = TRUE)
  )
)
fwrite(summary_dt, file.path(outdir, "DE_tumour_vs_normal_summary.tsv"), sep = "\t")

message("Differential expression complete: ", normalizePath(outdir))
