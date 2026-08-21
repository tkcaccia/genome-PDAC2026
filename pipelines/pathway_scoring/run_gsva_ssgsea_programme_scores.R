#!/usr/bin/env Rscript

# Patient-data-safe GSVA/ssGSEA programme scoring template.
#
# This script reads a gene-by-sample expression matrix and a GMT gene-set file,
# calculates ssGSEA scores with GSVA when available, optionally calculates GSVA
# scores, and converts programme scores into cohort-relative
# low/intermediate/high classes.

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    expression = NULL,
    gmt = NULL,
    metadata = NULL,
    gene_column = "gene",
    sample_column = "sample_id",
    condition_column = "condition",
    tumour_label = "Tumour",
    min_size = 5,
    max_size = 500,
    out_prefix = NULL
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- sub("^--", "", key)
    if (!name %in% names(out)) stop("Unknown argument: --", name)
    if (i == length(args)) stop("Missing value for argument: ", key)
    out[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  out$min_size <- as.integer(out$min_size)
  out$max_size <- as.integer(out$max_size)
  if (is.null(out$expression)) stop("--expression is required")
  if (is.null(out$gmt)) stop("--gmt is required")
  if (is.null(out$out_prefix)) stop("--out-prefix is required")
  out
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  gene_sets <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    genes <- unique(fields[-c(1, 2)])
    genes[genes != ""]
  })
  names(gene_sets) <- vapply(strsplit(lines, "\t", fixed = TRUE), `[`, character(1), 1)
  gene_sets
}

filter_gene_sets <- function(gene_sets, genes, min_size, max_size) {
  filtered <- lapply(gene_sets, intersect, genes)
  sizes <- vapply(filtered, length, integer(1))
  filtered[sizes >= min_size & sizes <= max_size]
}

read_expression <- function(path, gene_column) {
  expr_dt <- fread(path)
  if (!gene_column %in% names(expr_dt)) {
    stop("Gene column not found: ", gene_column)
  }
  genes <- expr_dt[[gene_column]]
  expr_dt[[gene_column]] <- NULL
  expr <- as.matrix(expr_dt)
  storage.mode(expr) <- "numeric"
  rownames(expr) <- make.unique(as.character(genes))
  expr
}

select_samples <- function(expr, metadata_path, sample_column, condition_column, tumour_label) {
  if (is.null(metadata_path)) return(expr)
  meta <- fread(metadata_path)
  if (!sample_column %in% names(meta)) stop("Sample column not found in metadata: ", sample_column)
  if (!condition_column %in% names(meta)) stop("Condition column not found in metadata: ", condition_column)
  tumour_samples <- meta[get(condition_column) == tumour_label, get(sample_column)]
  keep <- intersect(colnames(expr), tumour_samples)
  if (length(keep) == 0) stop("No tumour samples from metadata were found in expression matrix columns.")
  expr[, keep, drop = FALSE]
}

run_ssgsea <- function(expr, gene_sets) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("The GSVA R package is required for ssGSEA scoring.")
  }
  if (exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    param <- GSVA::ssgseaParam(expr, gene_sets, normalize = TRUE)
    return(GSVA::gsva(param, verbose = FALSE))
  }
  GSVA::gsva(expr, gene_sets, method = "ssgsea", kcdf = "Gaussian", abs.ranking = FALSE, verbose = FALSE)
}

run_gsva <- function(expr, gene_sets) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("The GSVA R package is required for GSVA scoring.")
  }
  if (exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    param <- GSVA::gsvaParam(expr, gene_sets, kcdf = "Gaussian")
    return(GSVA::gsva(param, verbose = FALSE))
  }
  GSVA::gsva(expr, gene_sets, method = "gsva", kcdf = "Gaussian", verbose = FALSE)
}

write_matrix <- function(mat, path) {
  out <- data.table(programme = rownames(mat), as.data.frame(mat, check.names = FALSE))
  fwrite(out, path, sep = "\t")
}

make_classes <- function(score_mat) {
  class_list <- lapply(seq_len(nrow(score_mat)), function(i) {
    scores <- as.numeric(score_mat[i, ])
    q <- stats::quantile(scores, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
    class <- ifelse(scores <= q[1], "low", ifelse(scores >= q[2], "high", "intermediate"))
    data.table(
      programme = rownames(score_mat)[i],
      sample_id = colnames(score_mat),
      score = scores,
      programme_class = class
    )
  })
  rbindlist(class_list)
}

make_summary <- function(classes) {
  classes[, .(
    median_score = median(score, na.rm = TRUE),
    min_score = min(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    low_n = sum(programme_class == "low", na.rm = TRUE),
    intermediate_n = sum(programme_class == "intermediate", na.rm = TRUE),
    high_n = sum(programme_class == "high", na.rm = TRUE)
  ), by = programme]
}

main <- function() {
  opt <- parse_args()
  dir.create(dirname(opt$out_prefix), recursive = TRUE, showWarnings = FALSE)

  expr <- read_expression(opt$expression, opt$gene_column)
  expr <- select_samples(expr, opt$metadata, opt$sample_column, opt$condition_column, opt$tumour_label)
  gene_sets <- read_gmt(opt$gmt)
  gene_sets <- filter_gene_sets(gene_sets, rownames(expr), opt$min_size, opt$max_size)
  if (length(gene_sets) == 0) stop("No gene sets passed min/max size filters after matching to expression genes.")

  notes <- c(
    paste("expression_file:", opt$expression),
    paste("gmt_file:", opt$gmt),
    paste("samples_scored:", ncol(expr)),
    paste("genes_in_expression:", nrow(expr)),
    paste("gene_sets_scored:", length(gene_sets)),
    paste("GSVA_installed:", requireNamespace("GSVA", quietly = TRUE))
  )

  ssgsea <- run_ssgsea(expr, gene_sets)
  write_matrix(ssgsea, paste0(opt$out_prefix, ".ssgsea_scores.tsv"))

  gsva_result <- tryCatch(run_gsva(expr, gene_sets), error = function(e) e)
  if (inherits(gsva_result, "error")) {
    notes <- c(notes, paste("GSVA_mode:", paste("skipped", gsva_result$message)))
  } else {
    write_matrix(gsva_result, paste0(opt$out_prefix, ".gsva_scores.tsv"))
    notes <- c(notes, "GSVA_mode: completed")
  }

  classes <- make_classes(ssgsea)
  fwrite(classes, paste0(opt$out_prefix, ".programme_classes.tsv"), sep = "\t")
  fwrite(make_summary(classes), paste0(opt$out_prefix, ".programme_summary.tsv"), sep = "\t")

  notes <- c(notes, paste("R_version:", R.version.string))
  if (requireNamespace("GSVA", quietly = TRUE)) {
    notes <- c(notes, paste("GSVA_version:", as.character(utils::packageVersion("GSVA"))))
  }
  writeLines(notes, paste0(opt$out_prefix, ".method_notes.txt"))
}

main()

