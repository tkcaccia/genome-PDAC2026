#!/usr/bin/env Rscript

# Create a focused, anonymized oncoprint from one patient-level integrated table.
# Germline entries are deliberately labelled as exploratory candidates because
# clinical classification and orthogonal confirmation are outside this script.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: plot_integrated_driver_oncoprint.R ",
    "<genotype_phenotype.tsv> <phenotype_assignments.tsv> <msisensor_summary.tsv> <outdir>"
  )
}

genotype_file <- args[[1]]
phenotype_file <- args[[2]]
msi_file <- args[[3]]
outdir <- args[[4]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

genotype <- fread(genotype_file, check.names = FALSE, colClasses = list(character = "patient_id"))
phenotype <- fread(phenotype_file, colClasses = list(character = "patient_id"))
msi <- fread(msi_file, colClasses = list(character = "patient_id"))

genes <- c(
  "KRAS", "TP53", "CDKN2A", "SMAD4", "BRCA1", "BRCA2",
  "PALB2", "ATM", "ARID1A", "RNF43", "GNAS", "ERBB2"
)
alteration_columns <- paste0(genes, "_alteration")
missing_columns <- setdiff(c("patient_id", alteration_columns), names(genotype))
if (length(missing_columns)) stop("Missing genotype columns: ", paste(missing_columns, collapse = ", "))

patient_values <- unique(genotype$patient_id)
patient_numeric <- suppressWarnings(as.numeric(patient_values))
patient_values <- patient_values[order(is.na(patient_numeric), patient_numeric, patient_values)]
patient_map <- data.table(
  patient_id = patient_values,
  anonymous_patient = sprintf("P%02d", seq_along(patient_values))
)

long <- melt(
  genotype[, c("patient_id", alteration_columns), with = FALSE],
  id.vars = "patient_id",
  variable.name = "gene_column",
  value.name = "alteration"
)
long[, gene := sub("_alteration$", "", gene_column)]
long[, alteration := fifelse(is.na(alteration) | trimws(alteration) == "", "not_detected", alteration)]

classify_event <- function(value) {
  value <- as.character(value)
  if (is.na(value) || value == "" || value == "not_detected") return("Not detected")
  categories <- c(
    somatic = grepl("Somatic", value, fixed = TRUE),
    germline = grepl("Germline", value, fixed = TRUE),
    cn_loss = grepl("CN_loss", value, fixed = TRUE),
    cn_gain = grepl("CN_gain", value, fixed = TRUE),
    sv = grepl("SV", value, fixed = TRUE)
  )
  broad_classes <- c(
    categories[["somatic"]], categories[["germline"]],
    categories[["cn_loss"]] || categories[["cn_gain"]], categories[["sv"]]
  )
  if (sum(broad_classes) > 1) return("Multiple alteration classes")
  if (categories[["somatic"]]) return("Somatic coding variant")
  if (categories[["germline"]]) return("Exploratory germline candidate")
  if (categories[["cn_loss"]]) return("Copy-number loss")
  if (categories[["cn_gain"]]) return("Copy-number gain")
  if (categories[["sv"]]) return("Structural variant")
  "Other exploratory event"
}

long[, event_class := vapply(alteration, classify_event, character(1))]
long[, biallelic_candidate := grepl("Biallelic_candidate", alteration, fixed = TRUE)]
long <- merge(long, patient_map, by = "patient_id", all.x = TRUE)

prevalence <- long[event_class != "Not detected", .(altered_n = uniqueN(patient_id)), by = gene]
prevalence <- merge(data.table(gene = genes), prevalence, by = "gene", all.x = TRUE)
prevalence[is.na(altered_n), altered_n := 0L]
prevalence[, gene_label := sprintf("%s  %d/%d", gene, altered_n, nrow(patient_map))]
long <- merge(long, prevalence[, .(gene, gene_label)], by = "gene", all.x = TRUE)

long[, anonymous_patient := factor(anonymous_patient, levels = patient_map$anonymous_patient)]
gene_labels <- prevalence[match(genes, gene), gene_label]
long[, gene_label := factor(gene_label, levels = rev(gene_labels))]

phenotype_annotation <- phenotype[, .(patient_id, phenotype_group)]
msi_annotation <- msi[, .(patient_id, msisensor_interpretation)]
annotation <- merge(patient_map, phenotype_annotation, by = "patient_id", all.x = TRUE)
annotation <- merge(annotation, msi_annotation, by = "patient_id", all.x = TRUE)
annotation[is.na(phenotype_group) | phenotype_group == "", phenotype_group := "Not assigned"]
annotation[is.na(msisensor_interpretation) | msisensor_interpretation == "", msisensor_interpretation := "Not available"]
annotation[, phenotype_display := fcase(
  phenotype_group == "ImmuneHigh_StromalLow", "Immune-high / stromal-low",
  phenotype_group == "StromalHigh_EMTHigh_ImmuneLow", "Stromal/EMT-high / immune-low",
  phenotype_group == "Intermediate_or_mixed", "Intermediate or mixed",
  default = "Not assigned"
)]
annotation[, msi_display := fcase(
  msisensor_interpretation == "MSS_low_MSI_score", "Low MSIsensor-pro score",
  msisensor_interpretation == "borderline_elevated_MSI_score", "Borderline elevated score",
  msisensor_interpretation == "MSI_high_by_common_msisensor_threshold", "Above common MSI-high threshold",
  default = "Not available"
)]

annotation_long <- rbindlist(list(
  annotation[, .(
    anonymous_patient,
    track = "RNA-defined TME phenotype",
    value = phenotype_display
  )],
  annotation[, .(
    anonymous_patient,
    track = "MSIsensor-pro category",
    value = msi_display
  )]
))
annotation_long[, anonymous_patient := factor(anonymous_patient, levels = patient_map$anonymous_patient)]
annotation_long[, track := factor(track, levels = c("MSIsensor-pro category", "RNA-defined TME phenotype"))]

event_colours <- c(
  "Not detected" = "#F2F0EA",
  "Somatic coding variant" = "#B84337",
  "Exploratory germline candidate" = "#D7A33D",
  "Copy-number loss" = "#4F7CA5",
  "Copy-number gain" = "#D97940",
  "Structural variant" = "#6A8D73",
  "Multiple alteration classes" = "#563C62",
  "Other exploratory event" = "#9C9A96"
)

annotation_colours <- c(
  "Immune-high / stromal-low" = "#397C86",
  "Stromal/EMT-high / immune-low" = "#B8683E",
  "Intermediate or mixed" = "#B7B1A5",
  "Not assigned" = "#E4E0D8",
  "Low MSIsensor-pro score" = "#D8D4CB",
  "Borderline elevated score" = "#E5B957",
  "Above common MSI-high threshold" = "#A43D36",
  "Not available" = "#EEECE7"
)

p_annotation <- ggplot(annotation_long, aes(anonymous_patient, track, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_manual(values = annotation_colours, name = "Cohort annotation") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "serif", base_size = 8.5) +
  theme(
    panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8), legend.position = "right",
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  )

p_oncoprint <- ggplot(long, aes(anonymous_patient, gene_label, fill = event_class)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_point(
    data = long[biallelic_candidate == TRUE],
    shape = 8, size = 2.7, colour = "white", stroke = 0.8, show.legend = FALSE
  ) +
  scale_fill_manual(values = event_colours, name = "Integrated alteration") +
  labs(
    title = "Focused pancreatic-cancer alteration profile",
    subtitle = "Gene labels show the number of tumours with any exploratory event; stars mark biallelic candidates",
    x = "Anonymized tumour", y = NULL
  ) +
  theme_minimal(base_family = "serif", base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    legend.position = "right",
    plot.margin = margin(0, 5.5, 5.5, 5.5)
  )

combined <- p_annotation / p_oncoprint + plot_layout(heights = c(0.23, 1.0), guides = "collect")
ggsave(file.path(outdir, "integrated_driver_oncoprint.png"), combined, width = 12.5, height = 8.5, dpi = 320)
ggsave(file.path(outdir, "integrated_driver_oncoprint.pdf"), combined, width = 12.5, height = 8.5)

safe_long <- long[, .(
  anonymous_patient, gene, event_class, biallelic_candidate
)]
safe_annotation <- annotation[, .(
  anonymous_patient, phenotype_group, msisensor_interpretation
)]
fwrite(safe_long, file.path(outdir, "integrated_driver_oncoprint_anonymized_events.tsv"), sep = "\t")
fwrite(safe_annotation, file.path(outdir, "integrated_driver_oncoprint_anonymized_annotations.tsv"), sep = "\t")
writeLines(
  c(
    paste("genotype_table_md5:", unname(tools::md5sum(genotype_file))),
    paste("phenotype_assignments_md5:", unname(tools::md5sum(phenotype_file))),
    paste("msisensor_summary_md5:", unname(tools::md5sum(msi_file))),
    "germline_label: exploratory germline candidate; no clinical classification is inferred",
    "event_integration: a tile can combine somatic, germline-candidate, copy-number and structural-variant evidence",
    "biallelic_symbol: star marks a pre-existing computational biallelic-candidate flag",
    paste("R_version:", R.version.string)
  ),
  file.path(outdir, "integrated_driver_oncoprint_method_notes.txt")
)

message("Wrote integrated driver oncoprint to: ", normalizePath(outdir))
