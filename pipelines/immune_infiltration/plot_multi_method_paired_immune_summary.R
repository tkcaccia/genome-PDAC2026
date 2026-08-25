#!/usr/bin/env Rscript

# Plot paired tumour-minus-normal changes from several immune-deconvolution
# methods without treating their incompatible native scales as percentages.
# Each row is standardized by that feature's across-sample standard deviation.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: plot_multi_method_paired_immune_summary.R ",
    "<score_dir> <paired_test_dir> <metadata.tsv> <outdir>"
  )
}

score_dir <- args[[1]]
paired_dir <- args[[2]]
metadata_file <- args[[3]]
outdir <- args[[4]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

metadata <- fread(metadata_file, colClasses = "character")
required_metadata <- c("sample_id", "patient_id", "condition")
if (!all(required_metadata %in% names(metadata))) {
  stop("Metadata must contain sample_id, patient_id and condition")
}

method_labels <- c(
  estimate = "ESTIMATE",
  mcp_counter = "MCP-counter",
  epic = "EPIC",
  xcell = "xCell",
  quantiseq = "quanTIseq"
)

# These populations provide a compact cross-method view. The script retains
# exact method labels because similarly named estimates are not interchangeable.
selected_features <- list(
  estimate = c("immune score", "stroma score", "tumor purity"),
  mcp_counter = c(
    "T cell", "T cell CD8+", "cytotoxicity score", "B cell",
    "Macrophage/Monocyte", "Endothelial cell", "Cancer associated fibroblast"
  ),
  epic = c(
    "B cell", "T cell CD4+", "T cell CD8+", "Macrophage",
    "Endothelial cell", "Cancer associated fibroblast"
  ),
  xcell = c(
    "immune score", "stroma score", "T cell CD8+", "Macrophage M2",
    "Cancer associated fibroblast"
  ),
  quantiseq = c(
    "B cell", "T cell CD8+", "T cell regulatory (Tregs)",
    "Macrophage M1", "Macrophage M2", "Myeloid dendritic cell",
    "uncharacterized cell"
  )
)

patient_values <- unique(metadata$patient_id)
patient_numeric <- suppressWarnings(as.numeric(patient_values))
patient_order <- order(is.na(patient_numeric), patient_numeric, patient_values)
patient_values <- patient_values[patient_order]
patient_map <- data.table(
  patient_id = patient_values,
  anonymous_patient = sprintf("P%02d", seq_along(patient_values))
)

delta_parts <- list()
test_parts <- list()
missing_features <- list()

for (method in names(method_labels)) {
  score_file <- file.path(score_dir, paste0(method, "_scores.tsv"))
  delta_file <- file.path(paired_dir, method, "paired_tumour_normal_deltas.tsv")
  test_file <- file.path(paired_dir, method, "paired_tumour_normal_paired_tests.tsv")
  if (!all(file.exists(c(score_file, delta_file, test_file)))) {
    stop("Missing score, delta or test file for method: ", method)
  }

  scores <- fread(score_file, check.names = FALSE)
  deltas <- fread(delta_file, colClasses = list(character = c("patient_id", "feature")))
  tests <- fread(test_file)
  wanted <- selected_features[[method]]
  available <- intersect(wanted, scores$feature)
  missing_features[[method]] <- setdiff(wanted, scores$feature)
  if (!length(available)) next

  score_long <- melt(
    scores[feature %in% available],
    id.vars = "feature",
    variable.name = "sample_id",
    value.name = "score"
  )
  score_scale <- score_long[, .(feature_sd = sd(score, na.rm = TRUE)), by = feature]
  score_scale[!is.finite(feature_sd) | feature_sd == 0, feature_sd := NA_real_]

  current <- merge(
    deltas[feature %in% available],
    score_scale,
    by = "feature",
    all.x = TRUE
  )
  current[, standardized_delta := delta_tumour_minus_normal / feature_sd]
  current <- merge(
    current,
    tests[, .(feature, n_pairs, wilcoxon_p, paired_t_p, fdr)],
    by = "feature",
    all.x = TRUE
  )
  current <- merge(current, patient_map, by = "patient_id", all.x = TRUE)
  current[, `:=`(method = method_labels[[method]], method_key = method)]
  delta_parts[[method]] <- current
}

delta_data <- rbindlist(delta_parts, use.names = TRUE, fill = TRUE)
if (!nrow(delta_data)) stop("No selected immune features were available")

method_order <- unname(method_labels)
delta_data[, method := factor(method, levels = method_order)]
delta_data[, anonymous_patient := factor(
  anonymous_patient,
  levels = patient_map$anonymous_patient
)]

feature_order <- unlist(lapply(names(method_labels), function(method_value) {
  present <- unique(delta_data[method_key == method_value, feature])
  paste(method_value, intersect(selected_features[[method_value]], present), sep = "::")
}), use.names = FALSE)
delta_data[, feature_plot := factor(
  paste(method_key, feature, sep = "::"),
  levels = rev(unique(feature_order))
)]

summary_data <- delta_data[, .(
  n_pairs = unique(n_pairs)[1],
  median_standardized_delta = median(standardized_delta, na.rm = TRUE),
  q25 = quantile(standardized_delta, 0.25, na.rm = TRUE),
  q75 = quantile(standardized_delta, 0.75, na.rm = TRUE),
  wilcoxon_p = unique(wilcoxon_p)[1],
  paired_t_p = unique(paired_t_p)[1],
  fdr = unique(fdr)[1]
), by = .(method, method_key, feature, feature_plot)]
summary_data[, support := fifelse(
  !is.na(fdr) & fdr < 0.05,
  "FDR < 0.05",
  fifelse(!is.na(wilcoxon_p) & wilcoxon_p < 0.05, "Nominal P < 0.05 only", "No nominal signal")
)]

p_heatmap <- ggplot(
  delta_data,
  aes(anonymous_patient, feature_plot, fill = standardized_delta)
) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_grid(method ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = function(value) sub("^.*::", "", value)) +
  scale_fill_gradient2(
    low = "#3273A8", mid = "#F7F3EA", high = "#B84337", midpoint = 0,
    limits = c(-3, 3), oob = scales::squish,
    name = "Tumour - normal\n(feature SD units)"
  ) +
  labs(
    title = "A  Paired immune and stromal changes",
    subtitle = "Rows are standardized within each method-feature; columns are matched patients",
    x = "Anonymized patient", y = NULL
  ) +
  theme_minimal(base_family = "serif", base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom"
  )

p_summary <- ggplot(
  summary_data,
  aes(median_standardized_delta, feature_plot)
) +
  geom_vline(xintercept = 0, colour = "#555555", linetype = "dashed") +
  geom_errorbar(
    aes(xmin = q25, xmax = q75), width = 0,
    orientation = "y", colour = "#666666"
  ) +
  geom_point(aes(fill = support), shape = 21, size = 3, colour = "#303030") +
  facet_grid(method ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = function(value) sub("^.*::", "", value)) +
  scale_fill_manual(values = c(
    "FDR < 0.05" = "#B84337",
    "Nominal P < 0.05 only" = "#E4A64A",
    "No nominal signal" = "#D5D2CA"
  )) +
  labs(
    title = "B  Cohort-level paired effects",
    subtitle = "Wilcoxon tests; BH FDR calculated within each method",
    x = "Median tumour - normal change (feature SD units)", y = NULL,
    fill = "Statistical support"
  ) +
  theme_minimal(base_family = "serif", base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom"
  )

combined <- p_heatmap + p_summary + plot_layout(widths = c(1.25, 1.0))
ggsave(file.path(outdir, "multi_method_paired_immune_summary.png"), combined, width = 15.5, height = 10.5, dpi = 320)
ggsave(file.path(outdir, "multi_method_paired_immune_summary.pdf"), combined, width = 15.5, height = 10.5)

anonymized_delta <- copy(delta_data)
anonymized_delta[, patient_id := NULL]
setcolorder(anonymized_delta, c(
  "anonymous_patient", "method", "feature", "Normal", "Tumour",
  "delta_tumour_minus_normal", "standardized_delta", "n_pairs",
  "wilcoxon_p", "paired_t_p", "fdr"
))
fwrite(anonymized_delta, file.path(outdir, "multi_method_paired_immune_anonymized_values.tsv"), sep = "\t")
fwrite(summary_data, file.path(outdir, "multi_method_paired_immune_summary.tsv"), sep = "\t")

missing_text <- unlist(lapply(names(missing_features), function(method) {
  missing <- missing_features[[method]]
  if (!length(missing)) return(NULL)
  paste(method, paste(missing, collapse = ";"), sep = "\t")
}))
writeLines(
  c(
    paste("score_directory:", normalizePath(score_dir)),
    paste("paired_test_directory:", normalizePath(paired_dir)),
    paste("metadata_md5:", unname(tools::md5sum(metadata_file))),
    "primary_test: paired two-sided Wilcoxon signed-rank test",
    "multiple_testing: Benjamini-Hochberg correction within each method",
    "visual_scale: tumour-minus-normal delta divided by the feature SD across all samples",
    "interpretation: native scales differ among methods and are not percentages unless a method explicitly returns fractions",
    paste("missing_selected_features:", ifelse(length(missing_text), paste(missing_text, collapse = " | "), "none")),
    paste("R_version:", R.version.string)
  ),
  file.path(outdir, "multi_method_paired_immune_method_notes.txt")
)

message("Wrote paired multi-method immune figure to: ", normalizePath(outdir))
