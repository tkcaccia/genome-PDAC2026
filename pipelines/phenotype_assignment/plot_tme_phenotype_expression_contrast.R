#!/usr/bin/env Rscript

# Characterize the expression-defined extreme TME groups. Because the groups
# were derived from RNA-seq scores, this contrast is descriptive rather than an
# independent biomarker-discovery test.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: plot_tme_phenotype_expression_contrast.R <limma_results.tsv> <outdir>")
}

results_file <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

results <- fread(results_file)
required <- c("gene_id", "gene_symbol", "logFC", "P.Value", "adj.P.Val", "t")
if (!all(required %in% names(results))) stop("Input is missing required limma result columns")

results[, gene_label := fifelse(
  !is.na(gene_symbol) & gene_symbol != "",
  gene_symbol,
  sub("\\..*$", "", gene_id)
)]
results[, significant := !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1]
results[, direction := fcase(
  significant & logFC > 0, "Higher in stromal/EMT-high group",
  significant & logFC < 0, "Higher in immune-high group",
  default = "Below threshold"
)]
results[, minus_log10_fdr := -log10(pmax(adj.P.Val, .Machine$double.xmin))]

top_labels <- results[significant == TRUE][order(P.Value)][seq_len(min(12L, .N))]
top_labels[, label_y := minus_log10_fdr + rep(c(0.12, 0.30, 0.48), length.out = .N)]
direction_colours <- c(
  "Higher in immune-high group" = "#397C86",
  "Below threshold" = "#C7C5BF",
  "Higher in stromal/EMT-high group" = "#B84337"
)

p_volcano <- ggplot(results, aes(logFC, minus_log10_fdr, colour = direction)) +
  geom_point(size = 0.65, alpha = 0.55) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "#777777") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#777777") +
  geom_text(
    data = top_labels,
    aes(y = label_y, label = gene_label),
    size = 2.6, vjust = 0, check_overlap = TRUE, show.legend = FALSE
  ) +
  scale_colour_manual(values = direction_colours) +
  labs(
    title = "A  Expression-defined TME group contrast",
    subtitle = "Stromal/EMT-high/immune-low versus immune-high/stromal-low tumours (3 versus 3)",
    x = "limma-voom log2 fold change",
    y = "-log10 FDR",
    colour = NULL
  ) +
  theme_minimal(base_family = "serif", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
    legend.position = "bottom", legend.text = element_text(size = 8)
  )

top_effects <- results[significant == TRUE][order(P.Value)][seq_len(min(20L, .N))]
top_effects[, standard_error := fifelse(is.finite(t) & t != 0, abs(logFC / t), NA_real_)]
top_effects[, lower_95 := logFC - 1.96 * standard_error]
top_effects[, upper_95 := logFC + 1.96 * standard_error]
top_effects[, gene_label := factor(gene_label, levels = rev(gene_label))]

p_effects <- ggplot(top_effects, aes(logFC, gene_label, colour = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#777777") +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95), width = 0,
    orientation = "y", linewidth = 0.45
  ) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = direction_colours) +
  labs(
    title = "B  Top descriptive expression differences",
    subtitle = "Approximate 95% intervals from moderated limma statistics",
    x = "Log2 fold change", y = NULL, colour = NULL
  ) +
  theme_minimal(base_family = "serif", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "italic", size = 8), legend.position = "none"
  )

footnote <- sprintf(
  "Descriptive result: %s of %s tested genes met FDR < 0.05 and |log2FC| > 1. Group labels were constructed from the same RNA-seq data and are not an independent validation cohort.",
  format(sum(results$significant, na.rm = TRUE), big.mark = ","),
  format(nrow(results), big.mark = ",")
)

combined <- p_volcano + p_effects +
  plot_annotation(caption = footnote) &
  theme(plot.caption = element_text(hjust = 0, size = 8.5, margin = margin(t = 8)))

ggsave(file.path(outdir, "tme_phenotype_expression_contrast.png"), combined, width = 14.5, height = 7.2, dpi = 320)
ggsave(file.path(outdir, "tme_phenotype_expression_contrast.pdf"), combined, width = 14.5, height = 7.2)

summary <- data.table(
  metric = c(
    "genes_tested", "FDR_lt_0.05", "FDR_lt_0.05_abs_log2FC_gt_1",
    "higher_stromal_EMT_group", "higher_immune_group", "minimum_FDR"
  ),
  value = c(
    nrow(results),
    sum(results$adj.P.Val < 0.05, na.rm = TRUE),
    sum(results$significant, na.rm = TRUE),
    sum(results$direction == "Higher in stromal/EMT-high group"),
    sum(results$direction == "Higher in immune-high group"),
    min(results$adj.P.Val, na.rm = TRUE)
  )
)
fwrite(summary, file.path(outdir, "tme_phenotype_expression_contrast_summary.tsv"), sep = "\t")
writeLines(
  c(
    paste("limma_result_md5:", unname(tools::md5sum(results_file))),
    "contrast: stromal/EMT-high/immune-low minus immune-high/stromal-low tumours",
    "sample_size: 3 tumours versus 3 tumours",
    "threshold: Benjamini-Hochberg FDR < 0.05 and absolute log2 fold change > 1",
    "critical_limitation: groups were derived from expression features in the same samples; the test is descriptive and circular for those defining programmes",
    paste("R_version:", R.version.string)
  ),
  file.path(outdir, "tme_phenotype_expression_contrast_method_notes.txt")
)

message("Wrote TME phenotype expression figure to: ", normalizePath(outdir))
