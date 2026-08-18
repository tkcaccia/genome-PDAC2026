#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: phenotype_group_comparison_limma_template.R <counts.tsv> <metadata.tsv> <outdir>\n",
    "metadata.tsv must contain sample_id, patient_id, condition and phenotype_group.\n",
    "condition must contain Tumour and Normal.\n",
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
required_cols <- c("sample_id", "patient_id", "condition", "phenotype_group")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

normalize_phenotype_group <- function(x) {
  x <- as.character(x)
  x[x %in% c("Immune_high_Stromal_low", "immune_high_stromal_low")] <- "ImmuneHigh_StromalLow"
  x[x %in% c("Stromal_EMT_high_Immune_low", "stromal_emt_high_immune_low")] <- "StromalHigh_EMTHigh_ImmuneLow"
  x
}

metadata <- metadata[condition %in% c("Tumour", "Normal")]
metadata[, condition := factor(condition, levels = c("Normal", "Tumour"))]
metadata[, patient_id := factor(patient_id)]
metadata[, phenotype_group := normalize_phenotype_group(phenotype_group)]
metadata[, phenotype_group := factor(phenotype_group)]

common_samples <- intersect(colnames(counts), metadata$sample_id)
if (length(common_samples) < 4) {
  stop("Too few overlapping samples between count matrix and metadata.", call. = FALSE)
}

counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[match(common_samples, sample_id)]

run_limma_voom <- function(count_matrix, meta, design, coef_name, output_file) {
  keep <- filterByExpr(count_matrix, design = design)
  dge <- DGEList(counts = count_matrix[keep, , drop = FALSE])
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  if (!coef_name %in% colnames(design)) {
    stop("Coefficient not found: ", coef_name, call. = FALSE)
  }
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  res <- data.table(gene_id = rownames(res), res)
  res[, significant_FDR_0_05_logFC_1 := adj.P.Val < 0.05 & abs(logFC) > 1]
  fwrite(res, output_file, sep = "\t")
  invisible(res)
}

# 1. Paired all-patient tumour versus normal comparison.
paired_design <- model.matrix(~ patient_id + condition, data = metadata)
run_limma_voom(
  counts,
  metadata,
  paired_design,
  "conditionTumour",
  file.path(outdir, "DE_all_patients_paired_tumour_vs_normal_limma.tsv")
)

# 2. Tumour-only phenotype contrast between the two extreme groups.
extreme_groups <- c("ImmuneHigh_StromalLow", "StromalHigh_EMTHigh_ImmuneLow")
tumour_meta <- metadata[condition == "Tumour" & phenotype_group %in% extreme_groups]
if (length(unique(tumour_meta$phenotype_group)) == 2 && nrow(tumour_meta) >= 4) {
  tumour_meta[, phenotype_group := factor(phenotype_group, levels = extreme_groups)]
  tumour_counts <- counts[, tumour_meta$sample_id, drop = FALSE]
  tumour_design <- model.matrix(~ phenotype_group, data = tumour_meta)
  run_limma_voom(
    tumour_counts,
    tumour_meta,
    tumour_design,
    "phenotype_groupStromalHigh_EMTHigh_ImmuneLow",
    file.path(outdir, "DE_tumour_stromal_emt_high_vs_immune_high_limma.tsv")
  )
} else {
  fwrite(
    data.table(reason = "Tumour-only extreme phenotype contrast skipped: fewer than two groups or too few samples."),
    file.path(outdir, "DE_tumour_stromal_emt_high_vs_immune_high_SKIPPED.tsv"),
    sep = "\t"
  )
}

# 3. Optional paired delta comparison: tumour-normal expression change by phenotype group.
v_all_design <- model.matrix(~ condition, data = metadata)
keep_all <- filterByExpr(counts, design = v_all_design)
dge_all <- DGEList(counts = counts[keep_all, , drop = FALSE])
dge_all <- calcNormFactors(dge_all, method = "TMM")
logcpm <- cpm(dge_all, log = TRUE, prior.count = 1)

paired_ids <- metadata[, .N, by = .(patient_id, condition)][, .N, by = patient_id][N == 2, patient_id]
delta_rows <- list()
for (pid in paired_ids) {
  normal_sample <- metadata[patient_id == pid & condition == "Normal", sample_id][1]
  tumour_sample <- metadata[patient_id == pid & condition == "Tumour", sample_id][1]
  group_value <- metadata[patient_id == pid & condition == "Tumour", phenotype_group][1]
  if (normal_sample %in% colnames(logcpm) && tumour_sample %in% colnames(logcpm)) {
    delta_rows[[as.character(pid)]] <- data.table(
      gene_id = rownames(logcpm),
      patient_id = as.character(pid),
      phenotype_group = as.character(group_value),
      delta_tumour_minus_normal = logcpm[, tumour_sample] - logcpm[, normal_sample]
    )
  }
}

delta_dt <- rbindlist(delta_rows)
if (nrow(delta_dt) > 0) {
  fwrite(delta_dt, file.path(outdir, "paired_logCPM_delta_by_patient.tsv"), sep = "\t")
}

message("Phenotype-group limma template complete: ", normalizePath(outdir))
