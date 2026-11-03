#!/usr/bin/env Rscript
# PDAC2026 analysis: Expression-to-deconvolution scores, paired comparisons and tumour-only Hallmark/TME groups.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
suppressPackageStartupMessages({
  library(quantiseqr)
  library(ggplot2)
  library(tidyr)
})

source(file.path(Sys.getenv("PDAC_CODE_ROOT"), "scripts/shared/runtime.R"))
root <- pdac_path("project_root")
out <- file.path(root, "results", "immune_salmon_11_20261012")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (folder in c("inputs", "tables", "figures", "qc")) dir.create(file.path(out, folder), showWarnings = FALSE)
patients <- pdac_cohort("primary_pdac")
samples <- as.vector(t(outer(patients, c("N", "T"), paste0)))
source <- file.path(root, "results/rnaseq_revise_11_archive/rnaseq_revised_11_20260916/tables/tpm_matrix_11_pdac.csv")
mapping_file <- file.path(root, "reproducible_pipeline_share_20260713_with_outputs/raw_outputs/results/figures/categories/04_cross_cohort_concordance/rnaseq_ensembl_to_symbol.csv")
tpm <- read.csv(source, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(samples %in% names(tpm)))
expression <- as.matrix(tpm[, samples])
storage.mode(expression) <- "double"
stopifnot(all(is.finite(expression)), all(expression >= 0))
ids <- sub("\\..*$", "", tpm[[1]])
expression <- rowsum(expression, ids, reorder = FALSE)
mapping <- read.csv(mapping_file, check.names = FALSE, stringsAsFactors = FALSE)
names(mapping)[1:2] <- c("ensembl", "symbol")
mapping$ensembl <- sub("\\..*$", "", mapping$ensembl)
stopifnot(!anyDuplicated(mapping$ensembl))
symbols <- mapping$symbol[match(rownames(expression), mapping$ensembl)]
mapped <- !is.na(symbols) & nzchar(trimws(symbols))
audit <- data.frame(ensembl = rownames(expression), symbol = symbols, mapped = mapped)
write.csv(audit, file.path(out, "qc/gene_mapping_audit.csv"), row.names = FALSE)
gene_tpm <- rowsum(expression[mapped, , drop = FALSE], symbols[mapped], reorder = FALSE)
write.csv(gene_tpm, file.path(out, "inputs/salmon_gene_symbol_tpm_11_pairs.csv"))
metadata <- data.frame(sample = samples, patient = sub("[NT]$", "", samples), condition = ifelse(grepl("T$", samples), "Tumour", "Normal"))
write.csv(metadata, file.path(out, "inputs/sample_metadata.csv"), row.names = FALSE)
write.csv(data.frame(sample = samples, original_tpm_total = colSums(expression), mapped_tpm_total = colSums(gene_tpm)), file.path(out, "qc/mapped_tpm_totals.csv"), row.names = FALSE)
signature_file <- system.file("extdata", "TIL10_signature.txt", package = "quantiseqr", mustWork = TRUE)
signature <- read.delim(signature_file, row.names = 1, check.names = FALSE)
removed <- read.delim(system.file("extdata", "TIL10_rmgenes.txt", package = "quantiseqr", mustWork = TRUE), header = FALSE)[[1]]
aberrant <- read.delim(system.file("extdata", "TIL10_TCGA_aberrant_immune_genes.txt", package = "quantiseqr", mustWork = TRUE), header = FALSE)[[1]]
effective <- setdiff(rownames(signature), union(removed, aberrant))
coverage <- data.frame(gene = rownames(signature), retained_in_signature = rownames(signature) %in% effective, present_in_tpm = rownames(signature) %in% rownames(gene_tpm))
write.csv(coverage, file.path(out, "qc/quantiseq_signature_coverage.csv"), row.names = FALSE)
coverage_fraction <- mean(effective %in% rownames(gene_tpm))
stopifnot(coverage_fraction >= 0.80)
results <- run_quantiseq(gene_tpm, signature_matrix = "TIL10", is_arraydata = FALSE, is_tumordata = TRUE, scale_mRNA = TRUE, method = "lsei", rm_genes = "default")
stopifnot(setequal(results$Sample, samples))
results <- results[match(samples, results$Sample), ]
fractions <- as.matrix(results[, -1])
stopifnot(all(is.finite(fractions)), all(fractions >= -1e-8), all(fractions <= 1 + 1e-8), all(abs(rowSums(fractions) - 1) < 1e-6))
write.csv(results, file.path(out, "tables/quantiseq_cell_fractions_11_pairs.csv"), row.names = FALSE)
write.csv(data.frame(sample = samples, fraction_sum = rowSums(fractions)), file.path(out, "qc/quantiseq_fraction_checks.csv"), row.names = FALSE)
cell_types <- names(results)[-1]
paired <- do.call(rbind, lapply(cell_types, function(cell_type) {
  normal <- results[match(paste0(patients, "N"), results$Sample), cell_type]
  tumour <- results[match(paste0(patients, "T"), results$Sample), cell_type]
  delta <- tumour - normal
  pvalue <- if (all(delta == 0)) 1 else wilcox.test(tumour, normal, paired = TRUE, alternative = "two.sided", exact = FALSE)$p.value
  data.frame(cell_type = cell_type, n_pairs = length(patients), normal_mean_percent = mean(normal) * 100, tumour_mean_percent = mean(tumour) * 100, median_delta_percentage_points = median(delta) * 100, p_value = pvalue)
}))
paired$FDR <- p.adjust(paired$p_value, method = "BH")
write.csv(paired, file.path(out, "tables/quantiseq_paired_comparison_11_pdac.csv"), row.names = FALSE)
long <- pivot_longer(results, cols = -Sample, names_to = "cell_type", values_to = "fraction")
long$patient <- factor(sub("[NT]$", "", long$Sample), levels = patients)
long$condition <- factor(ifelse(grepl("T$", long$Sample), "Tumour", "Normal"), levels = c("Normal", "Tumour"))
palette <- setNames(c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB", "#332288", "#88CCEE", "#44AA99", "#DDCC77"), cell_types)
stacked <- ggplot(long, aes(patient, fraction * 100, fill = cell_type)) +
  geom_col(width = 0.8) + facet_grid(. ~ condition) + scale_fill_manual(values = palette) +
  labs(x = "Patient", y = "Estimated cell fraction (%)", fill = "Cell type", title = "quanTIseq immune-cell composition: 11 paired PDAC samples", subtitle = "Salmon gene-level TPM; residual Other fraction retained") +
  theme_bw(base_size = 11) + theme(panel.grid.major.x = element_blank(), legend.position = "bottom")
ggsave(file.path(out, "figures/Figure_immune_composition_quantiseq_11_pdac.png"), stacked, width = 12, height = 6, dpi = 300)
paired_plot <- ggplot(long, aes(condition, fraction * 100, group = patient)) +
  geom_line(color = "#999999", alpha = 0.7) + geom_point(aes(color = condition), size = 1.6) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 4) + scale_color_manual(values = c(Normal = "#4477AA", Tumour = "#EE6677")) +
  labs(x = NULL, y = "Estimated cell fraction (%)", title = "Paired tumour-normal quanTIseq estimates: 11 PDAC patients") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out, "figures/Figure_immune_paired_changes_quantiseq_11_pdac.png"), paired_plot, width = 12, height = 8, dpi = 300)
capture.output(sessionInfo(), file = file.path(out, "qc/R_sessionInfo.txt"))
write.csv(data.frame(file = c(source, mapping_file, signature_file), md5 = unname(tools::md5sum(c(source, mapping_file, signature_file)))), file.path(out, "qc/input_checksums.csv"), row.names = FALSE)
writeLines(c("Primary cohort: 11 PDAC matched pairs; periampullary cases excluded.", "Input: Salmon gene-level TPM matrix with matched tumour and normal samples.", "Duplicate stripped Ensembl IDs and duplicate mapped symbols summed; missing symbols excluded and audited.", "Input is unlogged TPM. All 22 samples use is_arraydata=FALSE, is_tumordata=TRUE, scale_mRNA=TRUE, method=lsei, rm_genes=default so the same signature is used in both conditions.", "Tumour-mode removes aberrantly expressed signature genes in both conditions for primary paired comparability; method-dependent uncertainty remains.", "Other is a residual fraction, not a measured tumour-cell fraction or a distinct immune population.", "Paired tests: two-sided Wilcoxon with normal approximation (exact=FALSE, ties/zeros retained); BH across the 11 quanTIseq output types including Other.", "Cell fractions are estimates, not direct cell counts or validated immune infiltration measurements.", paste0("Effective signature coverage: ", sum(effective %in% rownames(gene_tpm)), "/", length(effective), " (", round(100 * coverage_fraction, 2), "%)."), "Hallmark programme scores and TME groups are calculated by generate_figureD_independent_programmes_11_pdac.py."), file.path(out, "METHODS_AND_LIMITATIONS.txt"))
cat("Samples:", nrow(results), "effective signature coverage:", sum(effective %in% rownames(gene_tpm)), "/", length(effective), "\n")
print(paired)
