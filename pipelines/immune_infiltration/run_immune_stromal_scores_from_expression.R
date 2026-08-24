#!/usr/bin/env Rscript

# Calculate immune and stromal scores from a patient-data-free expression matrix.
#
# immunedeconv expects HGNC gene symbols in rows, samples in columns, and
# non-log expression values. TPM is preferred. A log2(TPM + 1) matrix can be
# converted back to TPM by selecting --input-scale log2_plus_one. Do not use
# that option for VST, rlog, z-scores, or an unknown transformation.

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    expression = NULL,
    gene_column = "gene",
    input_scale = NULL,
    methods = "estimate,mcp_counter,epic,xcell,quantiseq",
    tumour_mode = "true",
    cibersort_script = NULL,
    cibersort_lm22 = NULL,
    cibersort_permutations = "100",
    out_dir = NULL
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

  required <- c("expression", "input_scale", "out_dir")
  missing <- required[vapply(out[required], is.null, logical(1))]
  if (length(missing) > 0) {
    stop("Missing required arguments: --", paste(gsub("_", "-", missing), collapse = ", --"))
  }
  if (!out$input_scale %in% c("linear", "log2_plus_one")) {
    stop("--input-scale must be 'linear' or 'log2_plus_one'.")
  }
  tumour_mode_value <- tolower(out$tumour_mode)
  if (!tumour_mode_value %in% c("true", "yes", "1", "false", "no", "0")) {
    stop("--tumour-mode must be true or false.")
  }
  out$tumour_mode <- tumour_mode_value %in% c("true", "yes", "1")
  out$cibersort_permutations <- as.integer(out$cibersort_permutations)
  out$methods <- trimws(strsplit(out$methods, ",", fixed = TRUE)[[1]])
  out
}

read_expression <- function(path, gene_column, input_scale) {
  dt <- fread(path)
  if (!gene_column %in% names(dt)) stop("Gene column not found: ", gene_column)

  genes <- trimws(as.character(dt[[gene_column]]))
  value_columns <- setdiff(names(dt), gene_column)
  if (length(value_columns) == 0) stop("The expression file contains no sample columns.")

  for (column in value_columns) set(dt, j = column, value = as.numeric(dt[[column]]))
  dt[, gene_symbol := genes]
  dt <- dt[!is.na(gene_symbol) & gene_symbol != ""]

  # Gene-set and deconvolution tools require one row per HGNC symbol.
  dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol, .SDcols = value_columns]
  expr <- as.matrix(dt[, ..value_columns])
  rownames(expr) <- dt$gene_symbol
  storage.mode(expr) <- "numeric"

  if (any(!is.finite(expr))) stop("Expression matrix contains non-finite values after parsing.")
  if (input_scale == "log2_plus_one") expr <- pmax(2^expr - 1, 0)
  if (any(expr < 0)) stop("Deconvolution input contains negative values; use non-negative linear TPM-like data.")
  expr
}

write_score_matrix <- function(result, method, out_dir) {
  result <- as.data.frame(result, check.names = FALSE)
  if (ncol(result) < 2) stop("Method returned no sample score columns: ", method)
  names(result)[1] <- "feature"
  fwrite(result, file.path(out_dir, paste0(method, "_scores.tsv")), sep = "\t")
}

run_one <- function(expr, method, opt) {
  if (method == "cibersort") {
    if (is.null(opt$cibersort_script) || is.null(opt$cibersort_lm22)) {
      stop("CIBERSORT requires --cibersort-script CIBERSORT.R and --cibersort-lm22 LM22.txt.")
    }
    immunedeconv::set_cibersort_binary(opt$cibersort_script)
    immunedeconv::set_cibersort_mat(opt$cibersort_lm22)
    return(immunedeconv::deconvolute(
      expr,
      method = "cibersort",
      tumor = opt$tumour_mode,
      arrays = FALSE,
      perm = opt$cibersort_permutations
    ))
  }

  immunedeconv::deconvolute(
    expr,
    method = method,
    tumor = opt$tumour_mode,
    arrays = FALSE
  )
}

main <- function() {
  opt <- parse_args()
  if (!requireNamespace("immunedeconv", quietly = TRUE)) {
    stop("The immunedeconv R package is required. See this folder's README for installation notes.")
  }

  allowed <- c("estimate", "mcp_counter", "epic", "xcell", "quantiseq", "cibersort")
  unknown <- setdiff(opt$methods, allowed)
  if (length(unknown) > 0) stop("Unsupported methods: ", paste(unknown, collapse = ", "))

  dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
  expr <- read_expression(opt$expression, opt$gene_column, opt$input_scale)

  status <- rbindlist(lapply(opt$methods, function(method) {
    message("Running ", method, "...")
    result <- tryCatch(run_one(expr, method, opt), error = function(e) e)
    if (inherits(result, "error")) {
      return(data.table(method = method, completed = FALSE, message = result$message))
    }
    write_score_matrix(result, method, opt$out_dir)
    data.table(method = method, completed = TRUE, message = "completed")
  }))

  fwrite(status, file.path(opt$out_dir, "method_status.tsv"), sep = "\t")
  writeLines(
    c(
      paste("expression_file:", opt$expression),
      paste("input_scale:", opt$input_scale),
      paste("genes_used:", nrow(expr)),
      paste("samples_scored:", ncol(expr)),
      paste("tumour_mode_for_EPIC_quanTIseq:", opt$tumour_mode),
      paste("R_version:", R.version.string),
      paste("immunedeconv_version:", as.character(utils::packageVersion("immunedeconv")))
    ),
    file.path(opt$out_dir, "method_notes.txt")
  )

  if (!all(status$completed)) {
    warning("One or more methods failed. See method_status.tsv; successful independent methods were retained.")
  }
  message("Immune/stromal scoring complete: ", normalizePath(opt$out_dir))
}

main()
