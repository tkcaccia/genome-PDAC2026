#!/usr/bin/env Rscript

# Compare genomic and technical features between RNA-defined extreme TME groups.
# With three tumours per extreme, results are exploratory and are corrected
# across the displayed metrics rather than interpreted as confirmatory tests.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: plot_tme_phenotype_genomic_qc_comparison.R ",
    "<genotype_phenotype.tsv> <validated_summary.tsv> <phenotype_assignments.tsv> <outdir>"
  )
}

genotype_file <- args[[1]]
validated_file <- args[[2]]
phenotype_file <- args[[3]]
outdir <- args[[4]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

genotype <- fread(genotype_file, check.names = FALSE, colClasses = list(character = "patient_id"))
validated <- fread(validated_file, colClasses = list(character = "patient_id"))
phenotype <- fread(phenotype_file, colClasses = list(character = "patient_id"))

genotype_columns <- c(
  "patient_id", "cnv_burden_altered_segments", "SV_total_records",
  "ASCAT_aberrant_cell_fraction_purity_estimate",
  "tumour_wes_mosdepth_mean_depth_genomewide_recal",
  "tumour_normal_depth_ratio_genomewide_recal"
)
validated_columns <- c("patient_id", "strict_rare_coding_TMB_per_Mb")
if (!all(genotype_columns %in% names(genotype))) stop("Genotype table is missing required columns")
if (!all(validated_columns %in% names(validated))) stop("Validated table is missing required columns")
if (!all(c("patient_id", "phenotype_group") %in% names(phenotype))) stop("Phenotype table is missing required columns")

data <- merge(genotype[, ..genotype_columns], validated[, ..validated_columns], by = "patient_id")
data <- merge(data, phenotype[, .(patient_id, phenotype_group)], by = "patient_id")
data <- data[phenotype_group %in% c("ImmuneHigh_StromalLow", "StromalHigh_EMTHigh_ImmuneLow")]
if (nrow(data) != 6L) warning("Expected six extreme-group tumours; observed ", nrow(data))

data[, group_label := fcase(
  phenotype_group == "ImmuneHigh_StromalLow", "Immune-high / stromal-low",
  phenotype_group == "StromalHigh_EMTHigh_ImmuneLow", "Stromal/EMT-high / immune-low"
)]

metric_labels <- c(
  strict_rare_coding_TMB_per_Mb = "Strict rare-coding TMB\n(mutations/Mb)",
  cnv_burden_altered_segments = "Altered CNA\nsegments",
  SV_total_records = "SV records\n(summary level)",
  ASCAT_aberrant_cell_fraction_purity_estimate = "ASCAT purity\nestimate",
  tumour_wes_mosdepth_mean_depth_genomewide_recal = "Tumour WES\nmean depth",
  tumour_normal_depth_ratio_genomewide_recal = "Tumour/normal\ndepth ratio"
)

for (column in names(metric_labels)) {
  set(data, j = column, value = as.numeric(data[[column]]))
}

long <- melt(
  data,
  id.vars = c("patient_id", "phenotype_group", "group_label"),
  measure.vars = names(metric_labels),
  variable.name = "metric",
  value.name = "value"
)
long[, value := as.numeric(value)]

tests <- long[, {
  first <- value[group_label == "Immune-high / stromal-low"]
  second <- value[group_label == "Stromal/EMT-high / immune-low"]
  p <- if (sum(is.finite(first)) >= 2L && sum(is.finite(second)) >= 2L) {
    wilcox.test(first, second, exact = FALSE)$p.value
  } else {
    NA_real_
  }
  .(
    n_immune = sum(is.finite(first)),
    n_stromal = sum(is.finite(second)),
    median_immune = median(first, na.rm = TRUE),
    median_stromal = median(second, na.rm = TRUE),
    wilcoxon_p = p
  )
}, by = metric]
tests[, fdr := p.adjust(wilcoxon_p, method = "BH")]
tests[, label := sprintf("P=%.2g; FDR=%.2g", wilcoxon_p, fdr)]

long <- merge(long, tests[, .(metric, label)], by = "metric")
long[, metric_label := factor(metric_labels[metric], levels = unname(metric_labels))]
long[, group_label := factor(
  group_label,
  levels = c("Immune-high / stromal-low", "Stromal/EMT-high / immune-low")
)]

group_colours <- c(
  "Immune-high / stromal-low" = "#397C86",
  "Stromal/EMT-high / immune-low" = "#B8683E"
)

plot <- ggplot(long, aes(group_label, value, colour = group_label)) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 2.7, alpha = 0.9
  ) +
  stat_summary(fun = median, geom = "crossbar", width = 0.48, linewidth = 0.6) +
  facet_wrap(~metric_label, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = group_colours) +
  labs(
    title = "Genomic and technical features across expression-defined TME extremes",
    subtitle = "Three tumours per group; Wilcoxon tests with BH correction across displayed metrics",
    x = NULL, y = NULL, colour = NULL
  ) +
  theme_minimal(base_family = "serif", base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  ) +
  geom_text(
    data = unique(long[, .(metric_label, label)]),
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.4, size = 3.0, family = "serif"
  )

ggsave(file.path(outdir, "tme_phenotype_genomic_qc_comparison.png"), plot, width = 12.5, height = 8.2, dpi = 320)
ggsave(file.path(outdir, "tme_phenotype_genomic_qc_comparison.pdf"), plot, width = 12.5, height = 8.2)
fwrite(tests, file.path(outdir, "tme_phenotype_genomic_qc_comparison_tests.tsv"), sep = "\t")
writeLines(
  c(
    paste("genotype_table_md5:", unname(tools::md5sum(genotype_file))),
    paste("validated_summary_md5:", unname(tools::md5sum(validated_file))),
    paste("phenotype_assignments_md5:", unname(tools::md5sum(phenotype_file))),
    "groups: RNA-defined cohort-relative extremes; three tumours per group",
    "statistics: two-sided Wilcoxon rank-sum tests; BH correction across six displayed metrics",
    "interpretation: exploratory; the group labels were derived from RNA expression, not from these genomic/QC metrics",
    "resolution: CNA and SV measures are WES-compatible summary-level estimates, not WGS-grade structural profiles",
    paste("R_version:", R.version.string)
  ),
  file.path(outdir, "tme_phenotype_genomic_qc_comparison_method_notes.txt")
)

message("Wrote TME phenotype genomic/QC figure to: ", normalizePath(outdir))
