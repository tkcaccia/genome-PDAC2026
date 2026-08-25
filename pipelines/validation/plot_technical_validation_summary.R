#!/usr/bin/env Rscript

# Summarize computational validation evidence without exposing source patient IDs.
# This is a sequencing-data review figure, not orthogonal clinical validation.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: plot_technical_validation_summary.R ",
    "<kras_validation.tsv> <hypermutation_validation.tsv> <msisensor_summary.tsv> <outdir>"
  )
}

kras_file <- args[[1]]
hypermutation_file <- args[[2]]
msi_file <- args[[3]]
outdir <- args[[4]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

kras <- fread(kras_file, colClasses = list(character = "patient_id"))
hypermutation <- fread(hypermutation_file, colClasses = list(character = "patient_id"))
msi <- fread(msi_file, colClasses = list(character = "patient_id"))

required_kras <- c("patient_id", "validation_verdict", "tumour_kras_region_median_depth")
required_hyper <- c(
  "patient_id", "previous_TMB_per_Mb_38Mb",
  "mutect2_strict_rare_TMB_per_Mb_38Mb", "relaxed_consensus_TMB_per_Mb_38Mb"
)
required_msi <- c("patient_id", "msi_score_percent", "msisensor_interpretation")
if (!all(required_kras %in% names(kras))) stop("KRAS validation input is missing required columns")
if (!all(required_hyper %in% names(hypermutation))) stop("Hypermutation input is missing required columns")
if (!all(required_msi %in% names(msi))) stop("MSIsensor-pro input is missing required columns")

patient_values <- Reduce(intersect, list(kras$patient_id, hypermutation$patient_id, msi$patient_id))
patient_numeric <- suppressWarnings(as.numeric(patient_values))
patient_values <- patient_values[order(is.na(patient_numeric), patient_numeric, patient_values)]
patient_map <- data.table(
  patient_id = patient_values,
  anonymous_patient = sprintf("P%02d", seq_along(patient_values))
)

kras <- merge(kras, patient_map, by = "patient_id")
hypermutation <- merge(hypermutation, patient_map, by = "patient_id")
msi <- merge(msi, patient_map, by = "patient_id")
for (object_name in c("kras", "hypermutation", "msi")) {
  object <- get(object_name)
  object[, anonymous_patient := factor(anonymous_patient, levels = patient_map$anonymous_patient)]
  assign(object_name, object)
}

kras[, verdict_label := fcase(
  validation_verdict == "KRAS_mutated_supported_by_VCF", "Coding variant supported by VCF",
  validation_verdict == "KRAS_no_coding_variant_detected_with_adequate_gene_coverage", "No coding variant; adequate coverage",
  validation_verdict == "KRAS_only_filtered_candidate_review_manually", "Filtered candidate only",
  default = gsub("_", " ", validation_verdict)
)]

kras_colours <- c(
  "Coding variant supported by VCF" = "#B84337",
  "No coding variant; adequate coverage" = "#397C86",
  "Filtered candidate only" = "#D9A441"
)
p_kras <- ggplot(kras, aes(anonymous_patient, tumour_kras_region_median_depth, fill = verdict_label)) +
  geom_col(width = 0.78, colour = "#303030", linewidth = 0.25) +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "#666666") +
  annotate("text", x = Inf, y = 20, label = "20x", hjust = 1.1, vjust = -0.4, size = 3) +
  scale_fill_manual(values = kras_colours) +
  labs(
    title = "A  KRAS coding-region review",
    subtitle = "Median tumour depth and VCF evidence category",
    x = NULL, y = "Median KRAS-region depth", fill = NULL
  ) +
  theme_minimal(base_family = "serif", base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom",
    legend.text = element_text(size = 7.5)
  )

msi[, msi_label := fcase(
  msisensor_interpretation == "MSS_low_MSI_score", "Low score",
  msisensor_interpretation == "borderline_elevated_MSI_score", "Borderline elevated",
  msisensor_interpretation == "MSI_high_by_common_msisensor_threshold", "Above common MSI-high threshold",
  default = gsub("_", " ", msisensor_interpretation)
)]
msi_colours <- c(
  "Low score" = "#BDB8AD",
  "Borderline elevated" = "#E0A943",
  "Above common MSI-high threshold" = "#A43D36"
)
p_msi <- ggplot(msi, aes(anonymous_patient, msi_score_percent, fill = msi_label)) +
  geom_col(width = 0.78, colour = "#303030", linewidth = 0.25) +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "#666666") +
  annotate("text", x = Inf, y = 20, label = "Common 20% threshold", hjust = 1.05, vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = msi_colours) +
  labs(
    title = "B  MSIsensor-pro screening",
    subtitle = "Computational score; tissue-based confirmation was unavailable",
    x = NULL, y = "Unstable microsatellite sites (%)", fill = NULL
  ) +
  theme_minimal(base_family = "serif", base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom",
    legend.text = element_text(size = 7.5)
  )

tmb_long <- melt(
  hypermutation,
  id.vars = c("patient_id", "anonymous_patient"),
  measure.vars = c(
    "previous_TMB_per_Mb_38Mb",
    "mutect2_strict_rare_TMB_per_Mb_38Mb",
    "relaxed_consensus_TMB_per_Mb_38Mb"
  ),
  variable.name = "filtering_strategy",
  value.name = "tmb_per_mb"
)
tmb_long[, filtering_strategy := factor(
  filtering_strategy,
  levels = c(
    "previous_TMB_per_Mb_38Mb",
    "relaxed_consensus_TMB_per_Mb_38Mb",
    "mutect2_strict_rare_TMB_per_Mb_38Mb"
  ),
  labels = c("Initial summary", "Mutect2-Strelka consensus", "Strict rare Mutect2")
)]

p_tmb <- ggplot(tmb_long, aes(filtering_strategy, tmb_per_mb, group = anonymous_patient)) +
  geom_line(colour = "#B9B5AC", linewidth = 0.45) +
  geom_point(aes(colour = filtering_strategy), size = 2.2) +
  scale_y_log10() +
  scale_colour_manual(values = c(
    "Initial summary" = "#8D8780",
    "Mutect2-Strelka consensus" = "#397C86",
    "Strict rare Mutect2" = "#B84337"
  )) +
  labs(
    title = "C  Mutation-burden sensitivity to filtering",
    subtitle = "Each line connects estimates for one anonymized tumour; y-axis is logarithmic",
    x = NULL, y = "Exploratory mutations per Mb (38-Mb denominator)", colour = NULL
  ) +
  theme_minimal(base_family = "serif", base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "none"
  )

combined <- (p_kras | p_msi) / p_tmb + plot_layout(heights = c(1.0, 0.9))
ggsave(file.path(outdir, "technical_validation_summary.png"), combined, width = 14.5, height = 10.2, dpi = 320)
ggsave(file.path(outdir, "technical_validation_summary.pdf"), combined, width = 14.5, height = 10.2)

aggregate_summary <- rbindlist(list(
  kras[, .(analysis = "KRAS review", category = verdict_label, n = .N), by = verdict_label][, verdict_label := NULL],
  msi[, .(analysis = "MSIsensor-pro", category = msi_label, n = .N), by = msi_label][, msi_label := NULL]
), use.names = TRUE)
fwrite(aggregate_summary, file.path(outdir, "technical_validation_aggregate_summary.tsv"), sep = "\t")
writeLines(
  c(
    paste("kras_validation_md5:", unname(tools::md5sum(kras_file))),
    paste("hypermutation_validation_md5:", unname(tools::md5sum(hypermutation_file))),
    paste("msisensor_summary_md5:", unname(tools::md5sum(msi_file))),
    "scope: computational sequencing-data review; not orthogonal pathology or clinical validation",
    "TMB_denominator: fixed exploratory 38-Mb exome denominator supplied by the source table",
    "MSI_threshold: dashed 20% line is the common MSIsensor threshold encoded by the source interpretation",
    "interpretation: mutation burden changed substantially with caller/filtering strategy and should not support broad hypermutation claims",
    paste("R_version:", R.version.string)
  ),
  file.path(outdir, "technical_validation_method_notes.txt")
)

message("Wrote technical validation figure to: ", normalizePath(outdir))
