#!/usr/bin/env Rscript

# Compare a primary paired differential-expression analysis with a prespecified
# anatomical-subset rerun. The subset DE table must be produced independently
# from the restricted samples, not by filtering the primary result table.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    "Usage: compare_anatomical_sensitivity.R <full_deseq.tsv> <subset_deseq.tsv> ",
    "<full_feature_deltas.tsv|-> <subset_feature_deltas.tsv|-> <outdir>",
    call. = FALSE
  )
}

full_file <- normalizePath(args[[1]], mustWork = TRUE)
subset_file <- normalizePath(args[[2]], mustWork = TRUE)
full_feature_file <- args[[3]]
subset_feature_file <- args[[4]]
outdir <- args[[5]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required <- function(table, columns, label) {
  missing <- setdiff(columns, names(table))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "))
}

full <- fread(full_file)
subset <- fread(subset_file)
columns <- c("gene_id", "DESeq2_log2FoldChange", "DESeq2_padj")
required(full, columns, "Full-cohort DE table")
required(subset, columns, "Anatomical-subset DE table")

joined <- merge(
  full[, .(
    gene_id,
    full_log2FC = DESeq2_log2FoldChange,
    full_FDR = DESeq2_padj
  )],
  subset[, .(
    gene_id,
    subset_log2FC = DESeq2_log2FoldChange,
    subset_FDR = DESeq2_padj
  )],
  by = "gene_id"
)

full_sig <- joined$full_FDR < 0.05 & abs(joined$full_log2FC) > 1
subset_sig <- joined$subset_FDR < 0.05 & abs(joined$subset_log2FC) > 1
same <- sign(joined$full_log2FC) == sign(joined$subset_log2FC)
summary <- data.table(
  common_genes = nrow(joined),
  pearson_r_log2FC = cor(joined$full_log2FC, joined$subset_log2FC, use = "complete.obs"),
  spearman_rho_log2FC = cor(joined$full_log2FC, joined$subset_log2FC, method = "spearman", use = "complete.obs"),
  same_direction_percent = 100 * mean(same, na.rm = TRUE),
  full_significant_in_overlap = sum(full_sig, na.rm = TRUE),
  full_significant_same_direction = sum(full_sig & same, na.rm = TRUE),
  significant_in_both = sum(full_sig & subset_sig, na.rm = TRUE),
  significant_in_both_same_direction = sum(full_sig & subset_sig & same, na.rm = TRUE)
)
fwrite(summary, file.path(outdir, "anatomical_DE_sensitivity_summary.tsv"), sep = "\t")

if (full_feature_file != "-" || subset_feature_file != "-") {
  if (full_feature_file == "-" || subset_feature_file == "-") {
    stop("Supply both feature-delta tables or use '-' for both", call. = FALSE)
  }
  full_features <- fread(normalizePath(full_feature_file, mustWork = TRUE))
  subset_features <- fread(normalizePath(subset_feature_file, mustWork = TRUE))
  feature_columns <- c("method", "feature", "median_delta", "p_value", "FDR")
  required(full_features, feature_columns, "Full-cohort feature table")
  required(subset_features, feature_columns, "Subset feature table")
  feature_join <- merge(
    full_features[, .(
      method,
      feature,
      full_median_delta = median_delta,
      full_p_value = p_value,
      full_FDR = FDR
    )],
    subset_features[, .(
      method,
      feature,
      subset_median_delta = median_delta,
      subset_p_value = p_value,
      subset_FDR = FDR
    )],
    by = c("method", "feature")
  )
  feature_join[, same_direction :=
    sign(full_median_delta) == sign(subset_median_delta)]
  feature_summary <- data.table(
    features_compared = nrow(feature_join),
    features_same_direction = sum(feature_join$same_direction, na.rm = TRUE),
    features_same_direction_percent = 100 * mean(feature_join$same_direction, na.rm = TRUE),
    full_FDR_significant = sum(feature_join$full_FDR < 0.05, na.rm = TRUE),
    subset_FDR_significant = sum(feature_join$subset_FDR < 0.05, na.rm = TRUE)
  )
  fwrite(feature_join, file.path(outdir, "anatomical_feature_sensitivity.tsv"), sep = "\t")
  fwrite(feature_summary, file.path(outdir, "anatomical_feature_sensitivity_summary.tsv"), sep = "\t")
}

writeLines(
  c(
    paste("full DE MD5:", unname(tools::md5sum(full_file))),
    paste("subset DE MD5:", unname(tools::md5sum(subset_file))),
    paste("R version:", R.version.string),
    "Interpretation: retained direction is descriptive and does not replace FDR-controlled inference."
  ),
  file.path(outdir, "anatomical_sensitivity_provenance.txt")
)

message("Sensitivity outputs written to: ", normalizePath(outdir))
