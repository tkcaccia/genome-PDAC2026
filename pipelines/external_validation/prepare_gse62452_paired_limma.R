#!/usr/bin/env Rscript

# Download GSE62452 from NCBI GEO and produce a gene-level paired limma result.
# GEO does not provide a patient-ID field for this series. Only the E-scheme
# titles with an unambiguous shared E-number are retained as matched pairs.

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(Biobase)
  library(GEOquery)
  library(hugene10sttranscriptcluster.db)
  library(limma)
  library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: prepare_gse62452_paired_limma.R <outdir>", call. = FALSE)
}

outdir <- args[[1]]
raw_dir <- file.path(outdir, "raw")
result_dir <- file.path(outdir, "results")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 600)

normalize_e_id <- function(value) {
  raw_id <- sub("^(E[0-9]+).*$", "\\1", value)
  paste0("E", as.integer(sub("^E0*", "", raw_id)))
}

select_one_probe_per_gene <- function(expression, annotation_package) {
  mapping <- AnnotationDbi::select(
    annotation_package,
    keys = rownames(expression),
    keytype = "PROBEID",
    columns = "SYMBOL"
  )
  mapping <- unique(mapping[!is.na(mapping$SYMBOL), c("PROBEID", "SYMBOL")])
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
  "GSE62452",
  destdir = raw_dir,
  GSEMatrix = TRUE,
  getGPL = FALSE
)
eset <- series[[1]]
phenotype <- Biobase::pData(eset)
probe_expression <- Biobase::exprs(eset)

required_metadata <- c("geo_accession", "title", "tissue:ch1")
missing_metadata <- setdiff(required_metadata, names(phenotype))
if (length(missing_metadata)) {
  stop("GSE62452 metadata is missing: ", paste(missing_metadata, collapse = ", "))
}

titles <- as.character(phenotype$title)
is_e_scheme <- grepl("^E[0-9]+-?T", titles)
metadata <- data.frame(
  sample_id = as.character(phenotype$geo_accession),
  title = titles,
  condition = ifelse(phenotype[["tissue:ch1"]] == "Pancreatic tumor", "tumor", "normal"),
  patient_id = NA_character_,
  stringsAsFactors = FALSE
)
metadata$patient_id[is_e_scheme] <- vapply(titles[is_e_scheme], normalize_e_id, character(1))

# A defensible pair must contain exactly one tumour and one normal sample for
# the same normalized E-number. Numeric-title samples are excluded because no
# patient identifier can be recovered from their GEO metadata.
candidate <- metadata[is_e_scheme, , drop = FALSE]
pair_ok <- vapply(split(candidate$condition, candidate$patient_id), function(condition) {
  length(condition) == 2 && setequal(condition, c("normal", "tumor"))
}, logical(1))
complete_patients <- names(pair_ok)[pair_ok]
paired_metadata <- candidate[candidate$patient_id %in% complete_patients, , drop = FALSE]
paired_metadata <- paired_metadata[order(paired_metadata$patient_id, paired_metadata$condition), ]
if (length(complete_patients) != 45) {
  warning("Expected 45 complete GSE62452 E-scheme pairs; found ", length(complete_patients))
}

paired_expression <- probe_expression[, paired_metadata$sample_id, drop = FALSE]
mapped <- select_one_probe_per_gene(
  paired_expression,
  hugene10sttranscriptcluster.db
)
gene_expression <- mapped$expression
metadata_model <- paired_metadata[match(colnames(gene_expression), paired_metadata$sample_id), ]
stopifnot(identical(metadata_model$sample_id, colnames(gene_expression)))
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
  file.path(result_dir, "GSE62452_limma_all_genes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  paired_metadata,
  file.path(result_dir, "GSE62452_paired_metadata.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  mapped$selected,
  file.path(result_dir, "GSE62452_selected_probe_per_gene.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
writeLines(
  c(
    paste("complete_pairs:", length(complete_patients)),
    paste("excluded_without_recoverable_pair:", nrow(metadata) - nrow(paired_metadata)),
    paste("genes_tested:", nrow(result)),
    paste("significant_FDR_0.05_abs_logFC_1:", sum(result$adj.P.Val < 0.05 & abs(result$logFC) > 1)),
    paste("R_version:", R.version.string),
    paste("GEOquery_version:", as.character(packageVersion("GEOquery"))),
    paste("limma_version:", as.character(packageVersion("limma"))),
    paste("hugene10sttranscriptcluster.db_version:", as.character(packageVersion("hugene10sttranscriptcluster.db")))
  ),
  file.path(result_dir, "GSE62452_provenance.txt")
)

message("GSE62452 paired limma outputs written to: ", normalizePath(result_dir))
