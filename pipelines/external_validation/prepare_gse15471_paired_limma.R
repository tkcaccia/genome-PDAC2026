#!/usr/bin/env Rscript

# Download GSE15471 from NCBI GEO and produce a gene-level paired limma result.
# Technical replicate arrays are averaged within patient and condition before
# modelling, so replicate arrays are never treated as independent samples.

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(Biobase)
  library(GEOquery)
  library(hgu133plus2.db)
  library(limma)
  library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: prepare_gse15471_paired_limma.R <outdir>", call. = FALSE)
}

outdir <- args[[1]]
raw_dir <- file.path(outdir, "raw")
result_dir <- file.path(outdir, "results")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 600)

select_one_probe_per_gene <- function(expression, annotation_package) {
  mapping <- AnnotationDbi::select(
    annotation_package,
    keys = rownames(expression),
    keytype = "PROBEID",
    columns = "SYMBOL"
  )
  mapping <- unique(mapping[!is.na(mapping$SYMBOL), c("PROBEID", "SYMBOL")])

  # A probe assigned to multiple symbols is removed rather than resolved
  # arbitrarily. For genes represented by several unambiguous probes, the
  # highest-variance probe is retained.
  symbol_count <- table(mapping$PROBEID)
  ambiguous <- names(symbol_count)[symbol_count > 1]
  mapping <- mapping[!(mapping$PROBEID %in% ambiguous), ]
  probe_expression <- expression[mapping$PROBEID, , drop = FALSE]
  candidates <- data.frame(
    PROBEID = mapping$PROBEID,
    SYMBOL = mapping$SYMBOL,
    variance = matrixStats::rowVars(probe_expression),
    stringsAsFactors = FALSE
  )
  candidates <- candidates[order(candidates$SYMBOL, -candidates$variance), ]
  selected <- candidates[!duplicated(candidates$SYMBOL), ]
  gene_expression <- probe_expression[selected$PROBEID, , drop = FALSE]
  rownames(gene_expression) <- selected$SYMBOL
  list(expression = gene_expression, selected = selected, ambiguous = ambiguous)
}

series <- GEOquery::getGEO(
  "GSE15471",
  destdir = raw_dir,
  GSEMatrix = TRUE,
  getGPL = FALSE
)
eset <- series[[1]]
phenotype <- Biobase::pData(eset)
probe_expression <- Biobase::exprs(eset)

required_metadata <- c("geo_accession", "patient:ch1", "sample:ch1")
missing_metadata <- setdiff(required_metadata, names(phenotype))
if (length(missing_metadata)) {
  stop("GSE15471 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}

metadata <- data.frame(
  sample_id = as.character(phenotype$geo_accession),
  patient_id = as.character(phenotype[["patient:ch1"]]),
  condition = ifelse(phenotype[["sample:ch1"]] == "tumor", "tumor", "normal"),
  stringsAsFactors = FALSE
)
metadata$group_id <- paste(metadata$patient_id, metadata$condition, sep = "__")

group_ids <- unique(metadata$group_id)
collapsed <- vapply(group_ids, function(group_id) {
  sample_ids <- metadata$sample_id[metadata$group_id == group_id]
  if (length(sample_ids) == 1) {
    probe_expression[, sample_ids]
  } else {
    rowMeans(probe_expression[, sample_ids, drop = FALSE])
  }
}, numeric(nrow(probe_expression)))
rownames(collapsed) <- rownames(probe_expression)
colnames(collapsed) <- group_ids

collapsed_metadata <- unique(metadata[, c("patient_id", "condition", "group_id")])
pair_counts <- table(collapsed_metadata$patient_id)
complete_patients <- names(pair_counts)[pair_counts == 2]
collapsed_metadata <- collapsed_metadata[
  collapsed_metadata$patient_id %in% complete_patients,
  ,
  drop = FALSE
]
collapsed <- collapsed[, collapsed_metadata$group_id, drop = FALSE]
if (length(complete_patients) != 36) {
  warning("Expected 36 complete GSE15471 pairs; found ", length(complete_patients))
}

mapped <- select_one_probe_per_gene(collapsed, hgu133plus2.db)
gene_expression <- mapped$expression
metadata_model <- collapsed_metadata[match(colnames(gene_expression), collapsed_metadata$group_id), ]
stopifnot(identical(metadata_model$group_id, colnames(gene_expression)))
metadata_model$patient_id <- factor(metadata_model$patient_id)
metadata_model$condition <- factor(metadata_model$condition, levels = c("normal", "tumor"))

design <- model.matrix(~ patient_id + condition, data = metadata_model)
fit <- limma::eBayes(limma::lmFit(gene_expression, design))
result <- limma::topTable(
  fit,
  coef = "conditiontumor",
  number = Inf,
  sort.by = "none"
)
result$gene <- rownames(result)
result <- result[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]

write.table(
  result,
  file.path(result_dir, "GSE15471_limma_all_genes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  collapsed_metadata,
  file.path(result_dir, "GSE15471_paired_metadata.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  mapped$selected,
  file.path(result_dir, "GSE15471_selected_probe_per_gene.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
writeLines(
  c(
    paste("complete_pairs:", length(complete_patients)),
    paste("genes_tested:", nrow(result)),
    paste("significant_FDR_0.05_abs_logFC_1:", sum(result$adj.P.Val < 0.05 & abs(result$logFC) > 1)),
    paste("R_version:", R.version.string),
    paste("GEOquery_version:", as.character(packageVersion("GEOquery"))),
    paste("limma_version:", as.character(packageVersion("limma"))),
    paste("hgu133plus2.db_version:", as.character(packageVersion("hgu133plus2.db")))
  ),
  file.path(result_dir, "GSE15471_provenance.txt")
)

message("GSE15471 paired limma outputs written to: ", normalizePath(result_dir))
