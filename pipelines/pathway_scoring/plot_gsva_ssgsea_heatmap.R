#!/usr/bin/env Rscript

# Plot a true GSVA or ssGSEA feature-by-sample score matrix. Patient IDs are
# joined from restricted metadata and replaced with ordered manuscript labels.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index)) return(default)
  if (index == length(args)) stop("Missing value for ", flag)
  args[[index + 1]]
}

scores_file <- get_arg("--scores")
metadata_file <- get_arg("--metadata")
out_dir <- get_arg("--out-dir")
score_type <- tolower(get_arg("--score-type", "ssgsea"))
feature_column <- get_arg("--feature-column", "programme")
sample_column <- get_arg("--sample-column", "sample_id")
patient_column <- get_arg("--patient-column", "patient_id")
phenotype_column <- get_arg("--phenotype-column", "phenotype_group")
msi_file <- get_arg("--msi-table", NULL)
msi_patient_column <- get_arg("--msi-patient-column", "patient_id")
msi_column <- get_arg("--msi-column", "validated_MSI_MMR_status")

if (any(vapply(list(scores_file, metadata_file, out_dir), is.null, logical(1)))) {
  stop("Required: --scores FILE --metadata FILE --out-dir DIR")
}
if (!score_type %in% c("ssgsea", "gsva")) stop("--score-type must be ssgsea or gsva")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scores <- fread(scores_file, check.names = FALSE)
if (!feature_column %in% names(scores)) stop("Missing feature column: ", feature_column)
metadata <- fread(metadata_file, colClasses = "character")
required_metadata <- c(sample_column, patient_column, phenotype_column)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0) stop("Metadata columns missing: ", paste(missing_metadata, collapse = ", "))

sample_columns <- setdiff(names(scores), feature_column)
metadata <- unique(metadata[get(sample_column) %in% sample_columns, ..required_metadata])
if (nrow(metadata) != length(sample_columns)) {
  stop("Every score-matrix sample must have exactly one metadata row")
}

numeric_patient <- suppressWarnings(as.integer(metadata[[patient_column]]))
metadata[, patient_sort := ifelse(is.na(numeric_patient), Inf, numeric_patient)]
setorderv(metadata, c("patient_sort", patient_column, sample_column), na.last = TRUE)
metadata[, tumour_label := sprintf("Tumour %02d", seq_len(.N))]

if (!is.null(msi_file)) {
  msi <- fread(msi_file, colClasses = "character")
  if (!all(c(msi_patient_column, msi_column) %in% names(msi))) {
    stop("MSI table does not contain the declared patient/status columns")
  }
  metadata <- merge(
    metadata,
    unique(msi[, .(patient_join = get(msi_patient_column), MSI_MMR_status = get(msi_column))]),
    by.x = patient_column,
    by.y = "patient_join",
    all.x = TRUE,
    sort = FALSE
  )
} else {
  metadata[, MSI_MMR_status := NA_character_]
}
setorderv(metadata, c("patient_sort", patient_column, sample_column), na.last = TRUE)

fwrite(
  metadata[, c(sample_column, patient_column, "tumour_label", phenotype_column, "MSI_MMR_status"), with = FALSE],
  file.path(out_dir, "anonymised_sample_annotation.tsv"),
  sep = "\t"
)

long <- melt(scores, id.vars = feature_column, variable.name = sample_column, value.name = "score")
long[, score := as.numeric(score)]
long <- merge(
  long,
  metadata[, c(sample_column, "tumour_label", phenotype_column, "MSI_MMR_status"), with = FALSE],
  by = sample_column,
  all.x = TRUE,
  sort = FALSE
)
setnames(long, c(feature_column, phenotype_column), c("programme", "phenotype_group"))
long[, z_score := if (sd(score, na.rm = TRUE) == 0) 0 else as.numeric(scale(score)), by = programme]
long[, programme_class := {
  cutoffs <- quantile(score, c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
  fifelse(score <= cutoffs[1], "low", fifelse(score >= cutoffs[2], "high", "intermediate"))
}, by = programme]

category_map <- list(
  "PDAC subtype" = c("PDAC_CLASSICAL_PROGENITOR", "PDAC_BASAL_SQUAMOUS", "PDAC_MESENCHYMAL", "PDAC_STROMAL_RICH"),
  "CAF/stroma" = c("CAF_ICAF", "CAF_MYCAF", "CAF_APCAF", "CAF_ECM", "TGF_BETA_FIBROBLAST"),
  "EMT/hypoxia/angiogenesis" = c("EMT_INVASION", "HYPOXIA", "ANGIOGENESIS"),
  "Metabolism" = c("GLYCOLYSIS", "OXIDATIVE_PHOSPHORYLATION", "PENTOSE_PHOSPHATE"),
  "Immune programme" = c("CYTOLYTIC_ACTIVITY", "ANTIGEN_PRESENTATION", "CHECKPOINT_EXPRESSION", "IFNG_RESPONSE")
)
category_lookup <- rbindlist(lapply(names(category_map), function(category) {
  data.table(programme = category_map[[category]], category = category)
}))
long <- merge(long, category_lookup, by = "programme", all.x = TRUE, sort = FALSE)
long[is.na(category), category := "Other"]

summary <- long[, .(
  category = category[1],
  median_score = median(score, na.rm = TRUE),
  min_score = min(score, na.rm = TRUE),
  max_score = max(score, na.rm = TRUE),
  low_n = sum(programme_class == "low"),
  intermediate_n = sum(programme_class == "intermediate"),
  high_n = sum(programme_class == "high")
), by = programme]
setorder(summary, category, programme)
fwrite(summary, file.path(out_dir, paste0(score_type, "_programme_summary.tsv")), sep = "\t")
fwrite(long, file.path(out_dir, paste0(score_type, "_programme_scores_long.tsv")), sep = "\t")

pretty_programme <- function(value) {
  labels <- c(
    "PDAC_CLASSICAL_PROGENITOR" = "PDAC classical/progenitor",
    "PDAC_BASAL_SQUAMOUS" = "PDAC basal/squamous",
    "PDAC_MESENCHYMAL" = "PDAC mesenchymal",
    "PDAC_STROMAL_RICH" = "PDAC stromal-rich",
    "CAF_ICAF" = "CAF iCAF",
    "CAF_MYCAF" = "CAF myCAF",
    "CAF_APCAF" = "CAF apCAF",
    "CAF_ECM" = "CAF ECM",
    "TGF_BETA_FIBROBLAST" = "TGF-beta fibroblast",
    "EMT_INVASION" = "EMT/invasion",
    "HYPOXIA" = "Hypoxia",
    "ANGIOGENESIS" = "Angiogenesis",
    "GLYCOLYSIS" = "Glycolysis",
    "OXIDATIVE_PHOSPHORYLATION" = "Oxidative phosphorylation",
    "PENTOSE_PHOSPHATE" = "Pentose phosphate pathway",
    "CYTOLYTIC_ACTIVITY" = "Cytolytic activity",
    "ANTIGEN_PRESENTATION" = "Antigen presentation",
    "CHECKPOINT_EXPRESSION" = "Checkpoint expression",
    "IFNG_RESPONSE" = "IFN-gamma response"
  )
  unname(labels[value])
}
long[, programme_label := pretty_programme(programme)]

phenotype_labels <- c(
  "StromalHigh_EMTHigh_ImmuneLow" = "Stromal/EMT-high\nimmune-low",
  "Intermediate_or_mixed" = "Intermediate/mixed",
  "ImmuneHigh_StromalLow" = "Immune-high\nstromal-low"
)
annotation <- unique(long[, .(tumour_label, phenotype_group, MSI_MMR_status)])
annotation[, TME_phenotype := phenotype_labels[phenotype_group]]
annotation[is.na(TME_phenotype), TME_phenotype := "Not assigned"]
annotation[, MSI_MMR := fifelse(
  is.na(MSI_MMR_status),
  "Not provided",
  fifelse(
    grepl("MSI-high|MMR-deficient", MSI_MMR_status, ignore.case = TRUE),
    "MSI-high/MMRd",
    fifelse(grepl("Borderline", MSI_MMR_status, ignore.case = TRUE), "Borderline MSI", "MSS/low MSI")
  )
)]

phenotype_order <- c("StromalHigh_EMTHigh_ImmuneLow", "Intermediate_or_mixed", "ImmuneHigh_StromalLow")
annotation[, phenotype_rank := match(phenotype_group, phenotype_order)]
setorder(annotation, phenotype_rank, tumour_label)
tumour_levels <- annotation$tumour_label

ann_long <- melt(
  annotation[, .(tumour_label, `TME phenotype` = TME_phenotype, `MSI/MMR` = MSI_MMR)],
  id.vars = "tumour_label",
  variable.name = "track",
  value.name = "annotation"
)
ann_long[, tumour_label := factor(tumour_label, levels = tumour_levels)]

annotation_colours <- c(
  "Stromal/EMT-high\nimmune-low" = "#B44C43",
  "Intermediate/mixed" = "#B7B7B7",
  "Immune-high\nstromal-low" = "#326A8C",
  "Not assigned" = "#E2E2E2",
  "MSI-high/MMRd" = "#9B3B73",
  "Borderline MSI" = "#D28CB3",
  "MSS/low MSI" = "#E8DFC8",
  "Not provided" = "#F0F0F0"
)

p_annotation <- ggplot(ann_long, aes(tumour_label, track, fill = annotation)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = annotation_colours, name = "Annotation") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "Times New Roman", base_size = 10) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold"), panel.grid = element_blank(),
    plot.margin = margin(3, 10, 0, 10)
  )

category_levels <- names(category_map)
programme_levels <- unique(category_lookup$programme)
programme_label_levels <- pretty_programme(programme_levels)
long[, category := factor(category, levels = category_levels)]
long[, programme_label := factor(programme_label, levels = rev(programme_label_levels))]
long[, tumour_label := factor(tumour_label, levels = tumour_levels)]

p_scores <- ggplot(long, aes(tumour_label, programme_label, fill = pmax(pmin(z_score, 2.5), -2.5))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#2B6F92", mid = "#F7F3EA", high = "#A83E36", midpoint = 0,
    limits = c(-2.5, 2.5), name = "Within-programme\nz-score"
  ) +
  labs(
    title = paste0(toupper(score_type), " RNA programme scores"),
    subtitle = "Raw method scores converted to within-programme z-scores for visualization",
    x = "Anonymised tumour", y = NULL
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 9.5, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8.5),
    strip.text.y = element_text(angle = 0, face = "bold", size = 8),
    panel.grid = element_blank(), plot.margin = margin(5, 10, 10, 10)
  )

combined <- p_annotation / p_scores + plot_layout(heights = c(0.12, 0.88), guides = "collect")
stem <- file.path(out_dir, paste0(score_type, "_programme_score_heatmap"))
ggsave(paste0(stem, ".png"), combined, width = 12.5, height = 9.5, dpi = 320)
ggsave(paste0(stem, ".pdf"), combined, width = 12.5, height = 9.5)

writeLines(
  c(
    paste("scores_file:", normalizePath(scores_file)),
    paste("scores_md5:", unname(tools::md5sum(scores_file))),
    paste("metadata_file:", normalizePath(metadata_file)),
    paste("metadata_md5:", unname(tools::md5sum(metadata_file))),
    paste("msi_file:", if (is.null(msi_file)) "not_used" else normalizePath(msi_file)),
    paste("msi_md5:", if (is.null(msi_file)) "not_used" else unname(tools::md5sum(msi_file))),
    paste("score_type:", score_type),
    paste("programmes:", uniqueN(long$programme)),
    paste("tumours:", uniqueN(long$tumour_label)),
    paste("R_version:", R.version.string)
  ),
  file.path(out_dir, paste0(score_type, "_heatmap_method_notes.txt"))
)

message("Wrote true ", score_type, " heatmap and tables to: ", normalizePath(out_dir))
