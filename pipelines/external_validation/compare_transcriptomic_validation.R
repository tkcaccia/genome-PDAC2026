#!/usr/bin/env Rscript

# Compare paired tumour-normal effects across internal methods and independent
# cohorts. The script writes results only to the caller-supplied output directory.

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(
    "Usage: compare_transcriptomic_validation.R <internal_deseq.tsv> ",
    "<internal_edger.tsv> <internal_limma.tsv> <external_manifest.tsv> ",
    "<programme.gmt> <outdir>",
    call. = FALSE
  )
}

deseq_file <- normalizePath(args[[1]], mustWork = TRUE)
edger_file <- normalizePath(args[[2]], mustWork = TRUE)
limma_file <- normalizePath(args[[3]], mustWork = TRUE)
manifest_file <- normalizePath(args[[4]], mustWork = TRUE)
gmt_file <- normalizePath(args[[5]], mustWork = TRUE)
outdir <- args[[6]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required <- function(table, columns, label) {
  missing <- setdiff(columns, names(table))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

safe_cor <- function(x, y, method = "pearson") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3) return(NA_real_)
  unname(cor(x[keep], y[keep], method = method))
}

direction_percent <- function(x, y) {
  keep <- is.finite(x) & is.finite(y) & x != 0 & y != 0
  if (!any(keep)) return(NA_real_)
  100 * mean(sign(x[keep]) == sign(y[keep]))
}

read_gmt <- function(path) {
  records <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  records <- records[vapply(records, length, integer(1)) >= 3]
  sets <- lapply(records, function(x) unique(toupper(trimws(x[-c(1, 2)]))))
  names(sets) <- vapply(records, `[[`, character(1), 1)
  sets
}

run_gsea <- function(statistic, pathways, cohort) {
  statistic <- statistic[is.finite(statistic) & !is.na(names(statistic)) & names(statistic) != ""]
  statistic <- sort(statistic[!duplicated(names(statistic))], decreasing = TRUE)
  set.seed(42)
  result <- as.data.table(fgseaMultilevel(
    pathways = pathways,
    stats = statistic,
    minSize = 3,
    maxSize = 500,
    eps = 0
  ))
  result[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
  result[, cohort := cohort]
  setcolorder(result, c("cohort", setdiff(names(result), "cohort")))
  result
}

deseq <- fread(deseq_file)
edger <- fread(edger_file)
limma_result <- fread(limma_file)
required(deseq, c("gene_id", "gene_symbol", "DESeq2_log2FoldChange", "stat", "DESeq2_padj"), "DESeq2 table")
required(edger, c("gene_id", "logFC", "FDR"), "edgeR table")
required(limma_result, c("gene_id", "logFC", "adj.P.Val"), "limma-voom table")

# A shared merge verifies that comparisons use the same gene universe rather
# than silently comparing different filtering outcomes.
internal <- Reduce(
  function(x, y) merge(x, y, by = "gene_id", all = FALSE),
  list(
    deseq[, .(
      gene_id,
      gene_symbol,
      deseq_log2FC = DESeq2_log2FoldChange,
      deseq_stat = stat,
      deseq_FDR = DESeq2_padj
    )],
    edger[, .(gene_id, edger_log2FC = logFC, edger_FDR = FDR)],
    limma_result[, .(gene_id, limma_log2FC = logFC, limma_FDR = adj.P.Val)]
  )
)

method_pairs <- list(
  c("DESeq2", "edgeR quasi-likelihood", "deseq_log2FC", "edger_log2FC"),
  c("DESeq2", "limma-voom", "deseq_log2FC", "limma_log2FC"),
  c("edgeR quasi-likelihood", "limma-voom", "edger_log2FC", "limma_log2FC")
)
method_concordance <- rbindlist(lapply(method_pairs, function(pair) {
  x <- internal[[pair[[3]]]]
  y <- internal[[pair[[4]]]]
  data.table(
    method_1 = pair[[1]],
    method_2 = pair[[2]],
    common_genes = sum(is.finite(x) & is.finite(y)),
    pearson_r = safe_cor(x, y),
    spearman_rho = safe_cor(x, y, "spearman"),
    same_direction_percent = direction_percent(x, y)
  )
}))

method_significance <- data.table(
  method = c("DESeq2", "edgeR quasi-likelihood", "limma-voom"),
  genes_in_common_universe = nrow(internal),
  genes_FDR_0_05_abs_log2FC_1 = c(
    sum(internal$deseq_FDR < 0.05 & abs(internal$deseq_log2FC) > 1, na.rm = TRUE),
    sum(internal$edger_FDR < 0.05 & abs(internal$edger_log2FC) > 1, na.rm = TRUE),
    sum(internal$limma_FDR < 0.05 & abs(internal$limma_log2FC) > 1, na.rm = TRUE)
  )
)
fwrite(method_concordance, file.path(outdir, "internal_method_concordance.tsv"), sep = "\t")
fwrite(method_significance, file.path(outdir, "internal_method_significance.tsv"), sep = "\t")

manifest <- fread(manifest_file)
required(manifest, c("cohort", "result_file"), "External manifest")
manifest_dir <- dirname(manifest_file)
manifest[, result_file := vapply(result_file, function(path) {
  candidate <- if (grepl("^/", path)) path else file.path(manifest_dir, path)
  normalizePath(candidate, mustWork = TRUE)
}, character(1))]

internal_symbol <- copy(internal)
internal_symbol[, gene_symbol_join := toupper(trimws(gene_symbol))]
internal_symbol <- internal_symbol[!is.na(gene_symbol_join) & gene_symbol_join != ""]
setorder(internal_symbol, deseq_FDR, na.last = TRUE)
internal_symbol <- unique(internal_symbol, by = "gene_symbol_join")

external_tables <- list()
validation_rows <- list()
for (i in seq_len(nrow(manifest))) {
  cohort <- as.character(manifest$cohort[[i]])
  external <- fread(manifest$result_file[[i]])
  symbol_columns <- intersect(c("gene_symbol", "gene", "SYMBOL"), names(external))
  if (!length(symbol_columns)) {
    stop(cohort, " table has no gene-symbol column (expected gene_symbol, gene or SYMBOL)")
  }
  symbol_column <- symbol_columns[[1]]
  required(external, c(symbol_column, "logFC", "t", "adj.P.Val"), paste0(cohort, " table"))
  external[, gene_symbol_join := toupper(trimws(as.character(get(symbol_column))))]
  external <- external[!is.na(gene_symbol_join) & gene_symbol_join != ""]
  setorder(external, adj.P.Val, na.last = TRUE)
  external <- unique(external, by = "gene_symbol_join")
  external_tables[[cohort]] <- external

  joined <- merge(
    internal_symbol[, .(
      gene_symbol = gene_symbol_join,
      internal_log2FC = deseq_log2FC,
      internal_FDR = deseq_FDR
    )],
    external[, .(
      gene_symbol = gene_symbol_join,
      external_log2FC = logFC,
      external_FDR = adj.P.Val
    )],
    by = "gene_symbol"
  )
  internal_sig <- joined$internal_FDR < 0.05 & abs(joined$internal_log2FC) > 1
  external_sig <- joined$external_FDR < 0.05 & abs(joined$external_log2FC) > 1
  same <- sign(joined$internal_log2FC) == sign(joined$external_log2FC)
  validation_rows[[cohort]] <- data.table(
    external_cohort = cohort,
    common_genes = nrow(joined),
    pearson_r_log2FC = safe_cor(joined$internal_log2FC, joined$external_log2FC),
    spearman_rho_log2FC = safe_cor(joined$internal_log2FC, joined$external_log2FC, "spearman"),
    same_direction_percent = 100 * mean(same, na.rm = TRUE),
    internal_significant_in_overlap = sum(internal_sig, na.rm = TRUE),
    internal_significant_same_direction_percent = 100 * mean(same[internal_sig], na.rm = TRUE),
    both_significant = sum(internal_sig & external_sig, na.rm = TRUE),
    both_significant_same_direction = sum(internal_sig & external_sig & same, na.rm = TRUE)
  )
}
fwrite(rbindlist(validation_rows), file.path(outdir, "internal_external_concordance.tsv"), sep = "\t")

# If two or more external cohorts are supplied, quantify the intersection shared
# by the internal result and every external cohort without selecting genes first.
if (length(external_tables) >= 2) {
  three_way <- internal_symbol[, .(
    gene_symbol = gene_symbol_join,
    internal_log2FC = deseq_log2FC,
    internal_FDR = deseq_FDR
  )]
  external_names <- names(external_tables)
  for (cohort in external_names) {
    table <- external_tables[[cohort]]
    part <- table[, .(
      gene_symbol = gene_symbol_join,
      log2FC = logFC,
      FDR = adj.P.Val
    )]
    setnames(part, c("log2FC", "FDR"), paste0(c("log2FC_", "FDR_"), cohort))
    three_way <- merge(three_way, part, by = "gene_symbol")
  }
  effect_columns <- c("internal_log2FC", paste0("log2FC_", external_names))
  fdr_columns <- c("internal_FDR", paste0("FDR_", external_names))
  effect_matrix <- as.matrix(three_way[, ..effect_columns])
  fdr_matrix <- as.matrix(three_way[, ..fdr_columns])
  significant_all <- apply(fdr_matrix < 0.05 & abs(effect_matrix) > 1, 1, all)
  same_all <- apply(sign(effect_matrix), 1, function(x) length(unique(x[x != 0])) == 1)
  summary <- data.table(
    cohorts_compared = paste(c("internal", external_names), collapse = ";"),
    common_genes = nrow(three_way),
    significant_in_all = sum(significant_all, na.rm = TRUE),
    significant_in_all_same_direction = sum(significant_all & same_all, na.rm = TRUE)
  )
  fwrite(summary, file.path(outdir, "all_cohort_overlap_summary.tsv"), sep = "\t")
}

pathways <- read_gmt(gmt_file)
gsea_results <- list()
internal_rank <- internal_symbol$deseq_stat
names(internal_rank) <- internal_symbol$gene_symbol_join
gsea_results[["internal"]] <- run_gsea(internal_rank, pathways, "internal")
for (cohort in names(external_tables)) {
  external <- external_tables[[cohort]]
  rank <- external$t
  names(rank) <- external$gene_symbol_join
  gsea_results[[cohort]] <- run_gsea(rank, pathways, cohort)
}
gsea <- rbindlist(gsea_results, use.names = TRUE, fill = TRUE)
fwrite(gsea, file.path(outdir, "exact_programme_gsea.tsv"), sep = "\t")

gsea_summary <- gsea[, .(
  programmes_tested = .N,
  programmes_FDR_0_05 = sum(padj < 0.05, na.rm = TRUE),
  programmes_positive_NES = sum(NES > 0, na.rm = TRUE)
), by = cohort]
fwrite(gsea_summary, file.path(outdir, "exact_programme_gsea_summary.tsv"), sep = "\t")

writeLines(
  c(
    paste("DESeq2 input MD5:", unname(tools::md5sum(deseq_file))),
    paste("edgeR input MD5:", unname(tools::md5sum(edger_file))),
    paste("limma input MD5:", unname(tools::md5sum(limma_file))),
    paste("external manifest MD5:", unname(tools::md5sum(manifest_file))),
    paste("programme GMT MD5:", unname(tools::md5sum(gmt_file))),
    paste("R version:", R.version.string),
    "Interpretation: effect agreement, FDR overlap and pathway agreement are distinct validation endpoints."
  ),
  file.path(outdir, "validation_provenance.txt")
)

message("Validation outputs written to: ", normalizePath(outdir))
