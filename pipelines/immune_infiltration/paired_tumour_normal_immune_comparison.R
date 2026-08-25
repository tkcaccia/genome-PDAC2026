#!/usr/bin/env Rscript

# Paired tumour-normal comparison for immune/stromal deconvolution scores.
#
# This script is intentionally generic and data-safe: it expects a score matrix
# and a metadata file supplied by the user, then writes only derived statistics.
# It does not contain patient data, sample identifiers from the PDAC cohort, or
# absolute private paths. The same approach was used for ESTIMATE, MCP-counter,
# CIBERSORT LM22, EPIC, xCell and quanTIseq score tables in the manuscript.
#
# Required input files:
#   1. scores.tsv
#      - Rows are immune/stromal features, e.g. "T cells" or "StromalScore".
#      - Columns are sample IDs.
#      - The first column must contain the feature name.
#   2. metadata.tsv
#      - Must contain sample_id, patient_id and condition.
#      - condition must use "Tumour" and "Normal".
#
# Output files:
#   - paired_tumour_normal_paired_tests.tsv:
#       one row per feature with tumour-minus-normal deltas, Wilcoxon P values,
#       paired t-test P values and Benjamini-Hochberg FDR values.
#   - paired_tumour_normal_deltas.tsv:
#       one row per patient-feature pair with the tumour-minus-normal delta.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: paired_tumour_normal_immune_comparison.R <scores.tsv> <metadata.tsv> <outdir> [feature_column_name]\n",
    "metadata.tsv must contain sample_id, patient_id and condition columns.",
    call. = FALSE
  )
}

scores_file <- args[[1]]
metadata_file <- args[[2]]
outdir <- args[[3]]
feature_col <- if (length(args) >= 4) args[[4]] else "feature"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load score matrix. The first column stores feature names; all other columns
# should be numeric sample-level deconvolution scores.
scores_dt <- fread(scores_file)
input_feature_col <- names(scores_dt)[1]
setnames(scores_dt, input_feature_col, feature_col)

# Convert wide score matrix to long format so each row is one feature/sample.
score_long <- melt(
  scores_dt,
  id.vars = feature_col,
  variable.name = "sample_id",
  value.name = "score"
)
score_long[, score := as.numeric(score)]

# Load and validate metadata. We keep only samples labelled Tumour or Normal,
# because the paired test requires a tumour-normal contrast within each patient.
metadata <- fread(metadata_file)
required_cols <- c("sample_id", "patient_id", "condition")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}
metadata <- metadata[condition %in% c("Tumour", "Normal")]

# Join scores to metadata. Samples missing from either file naturally drop out.
score_long <- merge(score_long, metadata[, ..required_cols], by = "sample_id")

# Keep complete tumour-normal pairs for each patient and feature. This prevents
# unmatched samples from entering the paired statistics.
# Build the formula from the requested feature-column name. For example, if the
# feature column is "cell_type", this becomes:
#   patient_id + cell_type ~ condition
paired_formula <- as.formula(paste("patient_id +", feature_col, "~ condition"))
paired <- dcast(score_long, paired_formula, value.var = "score")
paired <- paired[!is.na(Tumour) & !is.na(Normal)]
paired[, delta_tumour_minus_normal := Tumour - Normal]

# Run paired tests feature by feature. Wilcoxon signed-rank tests are used as
# the primary non-parametric comparison; paired t-tests are reported as a
# sensitivity statistic because some deconvolution scores are approximately
# continuous.
test_one_feature <- function(dt) {
  n_pairs <- nrow(dt)
  wilcox_p <- if (n_pairs >= 2) {
    suppressWarnings(wilcox.test(dt$Tumour, dt$Normal, paired = TRUE, exact = FALSE)$p.value)
  } else {
    NA_real_
  }
  paired_t_p <- if (n_pairs >= 2 && stats::sd(dt$delta_tumour_minus_normal) > 0) {
    t.test(dt$Tumour, dt$Normal, paired = TRUE)$p.value
  } else {
    NA_real_
  }
  data.table(
    n_pairs = n_pairs,
    mean_delta = mean(dt$delta_tumour_minus_normal, na.rm = TRUE),
    median_delta = median(dt$delta_tumour_minus_normal, na.rm = TRUE),
    wilcoxon_p = wilcox_p,
    paired_t_p = paired_t_p
  )
}

tests <- paired[, test_one_feature(.SD), by = feature_col]

# Correct the Wilcoxon P values across features within the supplied method.
# The manuscript interpreted nominal trends separately from FDR-significant
# findings because these methods test multiple immune/stromal features.
tests[, fdr := p.adjust(wilcoxon_p, method = "BH")]
setorderv(tests, c("wilcoxon_p", feature_col), na.last = TRUE)

fwrite(
  tests,
  file.path(outdir, "paired_tumour_normal_paired_tests.tsv"),
  sep = "\t"
)
fwrite(
  paired,
  file.path(outdir, "paired_tumour_normal_deltas.tsv"),
  sep = "\t"
)

message("Paired immune/stromal comparison complete: ", normalizePath(outdir))
