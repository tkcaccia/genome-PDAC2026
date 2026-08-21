#!/usr/bin/env Rscript

# Patient-data-safe table/figure generator for GSVA/ssGSEA-style programme scores.
#
# Input: an integrated tumour-level score table containing programme score
# columns, and optionally phenotype-group labels. Output: a programme summary
# table and a z-scaled programme heatmap.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[[idx + 1]]
}

score_table <- get_arg("--score-table")
out_dir <- get_arg("--out-dir")
sample_column <- get_arg("--sample-column", "sample_id")
phenotype_column <- get_arg("--phenotype-column", "phenotype_group")
msi_column <- get_arg("--msi-column", NULL)
if (is.null(score_table) || is.null(out_dir)) {
  stop("Usage: Rscript create_programme_score_table_figure.R --score-table scores.tsv --out-dir results/")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dt <- fread(score_table, check.names = TRUE)
if (!sample_column %in% names(dt)) stop("Missing sample column: ", sample_column)
if (!phenotype_column %in% names(dt)) dt[, (phenotype_column) := "not_provided"]
if (!is.null(msi_column) && !msi_column %in% names(dt)) stop("Missing MSI/MMR column: ", msi_column)

programme_map <- data.table(
  score_col = c(
    "PDAC_classical_progenitor_score", "PDAC_basal_squamous_score", "stromal_rich_score",
    "mesenchymal_score", "CAF_iCAF_score", "CAF_myCAF_score", "CAF_apCAF_score",
    "CAF_ECM_score", "EMT_invasion_score", "Hypoxia_score", "Angiogenesis_score",
    "Metabolism_Glycolysis_score", "Metabolism_OxPhos_score", "Metabolism_Pentose_phosphate_score",
    "Immune_Cytolytic_score", "Immune_Antigen_presentation_score", "Immune_Checkpoint_score",
    "Immune_IFNG_response_score"
  ),
  class_col = c(
    rep(NA_character_, 4), "CAF_iCAF_class", "CAF_myCAF_class", "CAF_apCAF_class",
    "CAF_ECM_class", "EMT_invasion_class", "Hypoxia_class", "Angiogenesis_class",
    "Metabolism_Glycolysis_class", "Metabolism_OxPhos_class", "Metabolism_Pentose_phosphate_class",
    "Immune_Cytolytic_class", "Immune_Antigen_presentation_class", "Immune_Checkpoint_class",
    "Immune_IFNG_response_class"
  ),
  programme = c(
    "PDAC classical/progenitor", "PDAC basal/squamous", "Stromal-rich programme",
    "Mesenchymal programme", "iCAF programme", "myCAF programme", "apCAF programme",
    "CAF extracellular matrix", "EMT/invasion", "Hypoxia", "Angiogenesis",
    "Glycolysis", "Oxidative phosphorylation", "Pentose phosphate",
    "Cytolytic activity", "Antigen presentation", "Checkpoint expression", "IFN-gamma response"
  ),
  category = c(
    rep("PDAC subtype", 4), rep("CAF/stroma", 4), rep("EMT/hypoxia/angiogenesis", 3),
    rep("Metabolism", 3), rep("Immune programme", 4)
  )
)
programme_map <- programme_map[score_col %in% names(dt)]
if (nrow(programme_map) == 0) stop("No recognized programme score columns were found.")

label_map <- data.table(
  sample_id_tmp = as.character(dt[[sample_column]]),
  tumour_label = sprintf("Tumour %02d", seq_len(nrow(dt)))
)
setnames(label_map, "sample_id_tmp", sample_column)
dt <- merge(dt, label_map, by = sample_column, all.x = TRUE, sort = FALSE)

long <- melt(
  dt,
  id.vars = c(sample_column, "tumour_label", phenotype_column),
  measure.vars = programme_map$score_col,
  variable.name = "score_col",
  value.name = "score"
)
setnames(long, phenotype_column, "phenotype_group")
long <- merge(long, programme_map, by = "score_col", all.x = TRUE)
long[, score := as.numeric(score)]
long[, z_score := {
  if (all(is.na(score)) || sd(score, na.rm = TRUE) == 0) rep(0, .N) else as.numeric(scale(score))
}, by = programme]

class_long <- rbindlist(lapply(seq_len(nrow(programme_map)), function(i) {
  class_col <- programme_map$class_col[i]
  if (is.na(class_col) || !class_col %in% names(dt)) {
    return(data.table(score_col = programme_map$score_col[i], sample_value = dt[[sample_column]], programme_class = NA_character_))
  }
  data.table(score_col = programme_map$score_col[i], sample_value = dt[[sample_column]], programme_class = as.character(dt[[class_col]]))
}))
setnames(class_long, "sample_value", sample_column)
long <- merge(long, class_long, by = c("score_col", sample_column), all.x = TRUE)

summary <- long[, .(
  category = category[1],
  median_score = median(score, na.rm = TRUE),
  min_score = min(score, na.rm = TRUE),
  max_score = max(score, na.rm = TRUE),
  high_n = sum(tolower(programme_class) == "high", na.rm = TRUE),
  intermediate_n = sum(tolower(programme_class) == "intermediate", na.rm = TRUE),
  low_n = sum(tolower(programme_class) == "low", na.rm = TRUE)
), by = .(programme, score_col)]
setorder(summary, category, programme)
fwrite(summary, file.path(out_dir, "programme_score_summary.tsv"), sep = "\t")
fwrite(long, file.path(out_dir, "programme_scores_long.tsv"), sep = "\t")

annotation <- unique(long[, .(tumour_label, phenotype_group)])
if (!is.null(msi_column)) {
  annotation <- merge(
    annotation,
    unique(dt[, .(tumour_label, MSI_MMR_status = as.character(get(msi_column)))]),
    by = "tumour_label",
    all.x = TRUE
  )
} else {
  annotation[, MSI_MMR_status := NA_character_]
}
phenotype_levels <- c("StromalHigh_EMTHigh_ImmuneLow", "Intermediate_or_mixed", "ImmuneHigh_StromalLow", "not_provided")
annotation[, phenotype_group := factor(phenotype_group, levels = phenotype_levels)]
annotation <- annotation[order(phenotype_group, tumour_label)]
tumour_levels <- annotation$tumour_label

pretty_phenotype <- c(
  "StromalHigh_EMTHigh_ImmuneLow" = "Stromal/EMT-high\nimmune-low",
  "Intermediate_or_mixed" = "Intermediate/\nmixed",
  "ImmuneHigh_StromalLow" = "Immune-high\nstromal-low",
  "not_provided" = "Not provided"
)
annotation[, TME_phenotype := pretty_phenotype[as.character(phenotype_group)]]
annotation[is.na(TME_phenotype), TME_phenotype := as.character(phenotype_group)]
annotation[, MSI_MMR := fifelse(
  is.na(MSI_MMR_status),
  NA_character_,
  fifelse(
    grepl("MSI-high|MMR-deficient", MSI_MMR_status, ignore.case = TRUE),
    "MSI-high/\nMMRd",
    fifelse(
      grepl("Borderline", MSI_MMR_status, ignore.case = TRUE),
      "Borderline\nMSI",
      fifelse(
        grepl("variant review", MSI_MMR_status, ignore.case = TRUE),
        "MSS + MMR\nvariant",
        "MSS/low\nMSI"
      )
    )
  )
)]

annotation_cols <- c("TME_phenotype")
if (!all(is.na(annotation$MSI_MMR))) annotation_cols <- c(annotation_cols, "MSI_MMR")
ann_long <- melt(
  annotation[, c("tumour_label", annotation_cols), with = FALSE],
  id.vars = "tumour_label",
  variable.name = "track",
  value.name = "annotation"
)
ann_long[, track := factor(track, levels = rev(annotation_cols), labels = rev(c("TME phenotype", "MSI/MMR")[seq_along(annotation_cols)]))]
ann_long[, tumour_label := factor(tumour_label, levels = tumour_levels)]

annotation_colours <- c(
  "Stromal/EMT-high\nimmune-low" = "#C44E52",
  "Intermediate/\nmixed" = "#B0B0B0",
  "Immune-high\nstromal-low" = "#4C72B0",
  "Not provided" = "#DDDDDD",
  "MSI-high/\nMMRd" = "#7B3294",
  "Borderline\nMSI" = "#C2A5CF",
  "MSS + MMR\nvariant" = "#F6E8C3",
  "MSS/low\nMSI" = "#E6E6E6"
)
annotation_colours <- annotation_colours[names(annotation_colours) %in% ann_long$annotation]

p_ann <- ggplot(ann_long, aes(tumour_label, track, fill = annotation)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = annotation_colours, name = "Annotation") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 9),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(5, 25, 0, 15)
  )

long[, tumour_label := factor(tumour_label, levels = tumour_levels)]
long[, programme := factor(programme, levels = rev(programme_map$programme))]

p <- ggplot(long, aes(tumour_label, programme, fill = pmax(pmin(z_score, 2.5), -2.5))) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-2.5, 2.5), name = "Row z-score") +
  labs(title = "RNA gene-set programme scores across tumours", subtitle = "GSVA/ssGSEA-style curated programme scores, z-scaled within each programme", x = "Anonymised tumour", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    strip.text.y = element_text(angle = 0, face = "bold", size = 8),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(15, 25, 15, 15)
  )

combined <- p_ann / p + plot_layout(heights = c(0.13, 0.87), guides = "collect")

ggsave(file.path(out_dir, "programme_score_heatmap.png"), combined, width = 12.5, height = 9.2, dpi = 300)
ggsave(file.path(out_dir, "programme_score_heatmap.pdf"), combined, width = 12.5, height = 9.2)

cat("Wrote programme summary and heatmap to:", out_dir, "\n")
