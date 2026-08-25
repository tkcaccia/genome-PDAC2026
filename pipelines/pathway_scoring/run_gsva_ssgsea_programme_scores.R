#!/usr/bin/env Rscript

# Reproducible GSVA/ssGSEA programme scoring from normalized RNA expression.
# The script supports Ensembl identifiers, records gene-set coverage, and writes
# the exact transformed expression matrix used for scoring.

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    expression = NULL,
    gmt = NULL,
    metadata = NULL,
    gtf = NULL,
    gene_column = "gene_id",
    sample_column = "sample_id",
    condition_column = "condition",
    tumour_label = "Tumour",
    transform = "auto",
    min_size = 5,
    max_size = 500,
    out_prefix = NULL
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- gsub("-", "_", sub("^--", "", key))
    if (!name %in% names(out)) stop("Unknown argument: ", key)
    if (i == length(args)) stop("Missing value for argument: ", key)
    out[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  out$min_size <- as.integer(out$min_size)
  out$max_size <- as.integer(out$max_size)
  if (is.null(out$expression)) stop("--expression is required")
  if (is.null(out$gmt)) stop("--gmt is required")
  if (is.null(out$out_prefix)) stop("--out-prefix is required")
  if (!out$transform %in% c("auto", "none", "log2p1")) {
    stop("--transform must be auto, none or log2p1")
  }
  out
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  invalid <- vapply(fields, length, integer(1)) < 3
  if (any(invalid)) stop("Invalid GMT line(s): ", paste(which(invalid), collapse = ", "))
  gene_sets <- lapply(fields, function(x) unique(toupper(x[-c(1, 2)][x[-c(1, 2)] != ""])))
  names(gene_sets) <- vapply(fields, `[`, character(1), 1)
  descriptions <- vapply(fields, `[`, character(1), 2)
  list(gene_sets = gene_sets, descriptions = descriptions)
}

parse_gtf_gene_map <- function(path) {
  # Select only the feature and attributes columns. This avoids retaining the
  # complete decompressed GTF as millions of R character strings.
  gtf <- fread(
    path,
    sep = "\t",
    header = FALSE,
    quote = "",
    select = c(3, 9),
    col.names = c("feature", "attributes"),
    showProgress = FALSE
  )
  gtf <- gtf[feature == "gene"]
  gene_ids <- sub('.*gene_id "([^"]+)".*', "\\1", gtf$attributes)
  gene_names <- sub('.*gene_name "([^"]+)".*', "\\1", gtf$attributes)
  valid <- gene_ids != gtf$attributes & gene_names != gtf$attributes & nzchar(gene_ids) & nzchar(gene_names)
  map <- data.table(
    gene_id = sub("\\..*$", "", gene_ids[valid]),
    gene_symbol = toupper(gene_names[valid])
  )
  unique(map[gene_id != "" & gene_symbol != ""], by = "gene_id")
}

read_expression <- function(path, gene_column, gtf_path, transform) {
  expr_dt <- fread(path)
  if (!gene_column %in% names(expr_dt)) stop("Gene column not found: ", gene_column)

  genes <- sub("\\..*$", "", as.character(expr_dt[[gene_column]]))
  expr_dt[[gene_column]] <- NULL
  sample_names <- names(expr_dt)
  expr <- as.matrix(expr_dt)
  storage.mode(expr) <- "numeric"
  if (any(!is.finite(expr))) stop("Expression matrix contains non-finite values")
  if (any(expr < 0) && transform == "log2p1") stop("log2p1 cannot be applied to negative values")

  ensembl_fraction <- mean(startsWith(genes, "ENSG"))
  mapped_fraction <- 0
  if (ensembl_fraction > 0.5) {
    if (is.null(gtf_path)) stop("Expression uses Ensembl IDs; provide --gtf for gene-symbol mapping")
    gene_map <- parse_gtf_gene_map(gtf_path)
    symbols <- gene_map$gene_symbol[match(genes, gene_map$gene_id)]
    mapped_fraction <- mean(!is.na(symbols))
    genes <- symbols
  } else {
    genes <- toupper(genes)
    mapped_fraction <- 1
  }

  keep <- !is.na(genes) & genes != ""
  expr <- expr[keep, , drop = FALSE]
  genes <- genes[keep]

  # Multiple Ensembl IDs can map to one symbol. Their normalized expression is
  # summed before log transformation so one biological gene is scored once.
  collapse_dt <- data.table(gene = genes, expr)
  collapsed <- collapse_dt[, lapply(.SD, sum, na.rm = TRUE), by = gene]
  collapsed_genes <- collapsed$gene
  collapsed[, gene := NULL]
  expr <- as.matrix(collapsed)
  rownames(expr) <- collapsed_genes
  colnames(expr) <- sample_names

  transform_used <- transform
  if (transform == "auto") {
    transform_used <- if (all(expr >= 0) && stats::quantile(expr, 0.99, na.rm = TRUE) > 50) {
      "log2p1"
    } else {
      "none"
    }
  }
  if (transform_used == "log2p1") expr <- log2(expr + 1)

  list(
    matrix = expr,
    ensembl_fraction = ensembl_fraction,
    mapped_fraction = mapped_fraction,
    transform_used = transform_used
  )
}

select_samples <- function(expr, metadata_path, sample_column, condition_column, tumour_label) {
  if (is.null(metadata_path)) return(expr)
  meta <- fread(metadata_path)
  required <- c(sample_column, condition_column)
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0) stop("Metadata columns missing: ", paste(missing, collapse = ", "))
  tumour_samples <- meta[get(condition_column) == tumour_label, get(sample_column)]
  keep <- intersect(colnames(expr), tumour_samples)
  if (length(keep) == 0) stop("No tumour samples from metadata were found in the expression matrix")
  expr[, keep, drop = FALSE]
}

filter_gene_sets <- function(gene_sets, genes, min_size, max_size) {
  matched <- lapply(gene_sets, intersect, genes)
  sizes <- vapply(matched, length, integer(1))
  list(sets = matched[sizes >= min_size & sizes <= max_size], sizes = sizes)
}

run_ssgsea <- function(expr, gene_sets) {
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("The GSVA R package is required")
  if (exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    param <- GSVA::ssgseaParam(expr, gene_sets, normalize = TRUE)
    return(GSVA::gsva(param, verbose = FALSE))
  }
  GSVA::gsva(expr, gene_sets, method = "ssgsea", kcdf = "Gaussian", abs.ranking = FALSE, verbose = FALSE)
}

run_gsva <- function(expr, gene_sets) {
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("The GSVA R package is required")
  if (exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    param <- GSVA::gsvaParam(expr, gene_sets, kcdf = "Gaussian")
    return(GSVA::gsva(param, verbose = FALSE))
  }
  GSVA::gsva(expr, gene_sets, method = "gsva", kcdf = "Gaussian", verbose = FALSE)
}

write_matrix <- function(mat, path, first_column) {
  out <- data.table(value = rownames(mat), as.data.frame(mat, check.names = FALSE))
  setnames(out, "value", first_column)
  fwrite(out, path, sep = "\t")
}

make_classes <- function(score_mat) {
  rbindlist(lapply(seq_len(nrow(score_mat)), function(i) {
    scores <- as.numeric(score_mat[i, ])
    q <- stats::quantile(scores, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE)
    data.table(
      programme = rownames(score_mat)[i],
      sample_id = colnames(score_mat),
      score = scores,
      programme_class = ifelse(scores <= q[1], "low", ifelse(scores >= q[2], "high", "intermediate"))
    )
  }))
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

  expression_info <- read_expression(opt$expression, opt$gene_column, opt$gtf, opt$transform)
  expr <- select_samples(
    expression_info$matrix,
    opt$metadata,
    opt$sample_column,
    opt$condition_column,
    opt$tumour_label
  )
  gmt <- read_gmt(opt$gmt)
  filtered <- filter_gene_sets(gmt$gene_sets, rownames(expr), opt$min_size, opt$max_size)
  if (length(filtered$sets) == 0) stop("No gene sets passed coverage and size filters")

  coverage <- data.table(
    programme = names(gmt$gene_sets),
    description = unname(gmt$descriptions),
    genes_in_gmt = vapply(gmt$gene_sets, length, integer(1)),
    genes_matched = unname(filtered$sizes)
  )
  coverage[, retained := programme %in% names(filtered$sets)]
  fwrite(coverage, paste0(opt$out_prefix, ".gene_set_coverage.tsv"), sep = "\t")
  write_matrix(expr, paste0(opt$out_prefix, ".expression_used.tsv"), "gene_symbol")

  ssgsea <- run_ssgsea(expr, filtered$sets)
  write_matrix(ssgsea, paste0(opt$out_prefix, ".ssgsea_scores.tsv"), "programme")

  gsva_result <- tryCatch(run_gsva(expr, filtered$sets), error = function(e) e)
  mode_notes <- c()
  if (inherits(gsva_result, "error")) {
    mode_notes <- c(mode_notes, paste("GSVA_mode: skipped", gsva_result$message))
  } else {
    write_matrix(gsva_result, paste0(opt$out_prefix, ".gsva_scores.tsv"), "programme")
    mode_notes <- c(mode_notes, "GSVA_mode: completed")
  }

  classes <- make_classes(ssgsea)
  fwrite(classes, paste0(opt$out_prefix, ".programme_classes.tsv"), sep = "\t")
  fwrite(make_summary(classes), paste0(opt$out_prefix, ".programme_summary.tsv"), sep = "\t")

  notes <- c(
    paste("expression_file:", normalizePath(opt$expression)),
    paste("expression_md5:", unname(tools::md5sum(opt$expression))),
    paste("gmt_file:", normalizePath(opt$gmt)),
    paste("gmt_md5:", unname(tools::md5sum(opt$gmt))),
    paste("gtf_file:", if (is.null(opt$gtf)) "not_used" else normalizePath(opt$gtf)),
    paste("gtf_md5:", if (is.null(opt$gtf)) "not_used" else unname(tools::md5sum(opt$gtf))),
    paste("samples_scored:", ncol(expr)),
    paste("genes_scored:", nrow(expr)),
    paste("input_ensembl_fraction:", sprintf("%.4f", expression_info$ensembl_fraction)),
    paste("mapped_gene_fraction:", sprintf("%.4f", expression_info$mapped_fraction)),
    paste("transform_used:", expression_info$transform_used),
    paste("gene_sets_retained:", length(filtered$sets)),
    paste("R_version:", R.version.string),
    paste("GSVA_version:", as.character(utils::packageVersion("GSVA"))),
    mode_notes
  )
  writeLines(notes, paste0(opt$out_prefix, ".method_notes.txt"))
}

main()
