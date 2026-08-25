#!/usr/bin/env Rscript

# Convert a gene-by-sample count matrix to gene-level TPM using union exon
# lengths from a declared GTF. This is a reproducible fallback when transcript-
# level Salmon TPM is unavailable; Salmon TPM remains preferable when complete.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: make_gene_tpm_from_counts.R <counts.tsv> <annotation.gtf[.gz]> <out_prefix>")
}

counts_file <- args[[1]]
gtf_file <- args[[2]]
out_prefix <- args[[3]]
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

message("Reading count matrix...")
counts_dt <- fread(counts_file)
gene_column <- names(counts_dt)[1]
gene_ids <- sub("\\..*$", "", as.character(counts_dt[[gene_column]]))
counts_dt[[gene_column]] <- NULL
counts <- as.matrix(counts_dt)
storage.mode(counts) <- "numeric"
if (any(!is.finite(counts)) || any(counts < 0)) stop("Counts must be finite and non-negative")

message("Reading GTF exons...")
gtf <- fread(
  gtf_file,
  sep = "\t",
  header = FALSE,
  quote = "",
  select = c(1, 3, 4, 5, 9),
  col.names = c("chromosome", "feature", "start", "end", "attributes")
)
gtf <- gtf[feature == "exon"]
gtf[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attributes)]
gtf[, gene_symbol := sub('.*gene_name "([^"]+)".*', "\\1", attributes)]
gtf[, gene_id := sub("\\..*$", "", gene_id)]
gtf <- gtf[gene_id != attributes & gene_symbol != attributes]

message("Calculating union exon lengths...")
exons <- GRanges(
  seqnames = gtf$chromosome,
  ranges = IRanges(start = gtf$start, end = gtf$end),
  gene_id = gtf$gene_id
)
exons_by_gene <- split(exons, exons$gene_id)
union_lengths <- vapply(reduce(exons_by_gene), function(x) sum(width(x)), numeric(1))

symbol_map <- unique(gtf[, .(gene_id, gene_symbol)], by = "gene_id")
length_index <- match(gene_ids, names(union_lengths))
symbol_index <- match(gene_ids, symbol_map$gene_id)
keep <- !is.na(length_index) & !is.na(symbol_index) & union_lengths[length_index] > 0
if (mean(keep) < 0.8) warning("Fewer than 80% of count-matrix genes matched GTF exon lengths")

counts <- counts[keep, , drop = FALSE]
gene_ids <- gene_ids[keep]
gene_symbols <- toupper(symbol_map$gene_symbol[symbol_index[keep]])
length_bp <- union_lengths[length_index[keep]]

rpk <- counts / (length_bp / 1000)
denominators <- colSums(rpk)
if (any(denominators <= 0)) stop("At least one sample has no positive expression")
tpm <- sweep(rpk, 2, denominators, "/") * 1e6

# Sum genes that share the same HGNC symbol after Ensembl mapping.
tpm_dt <- data.table(gene_symbol = gene_symbols, tpm)
tpm_by_symbol <- tpm_dt[gene_symbol != "", lapply(.SD, sum), by = gene_symbol]
setnames(tpm_by_symbol, names(tpm_by_symbol)[-1], colnames(counts))
fwrite(tpm_by_symbol, paste0(out_prefix, ".gene_tpm.tsv"), sep = "\t")

mapping <- data.table(
  gene_id = gene_ids,
  gene_symbol = gene_symbols,
  union_exon_length_bp = as.numeric(length_bp)
)
fwrite(mapping, paste0(out_prefix, ".gene_length_mapping.tsv"), sep = "\t")

writeLines(
  c(
    paste("counts_file:", normalizePath(counts_file)),
    paste("counts_md5:", unname(tools::md5sum(counts_file))),
    paste("gtf_file:", normalizePath(gtf_file)),
    paste("gtf_md5:", unname(tools::md5sum(gtf_file))),
    paste("input_genes:", length(keep)),
    paste("genes_with_length_and_symbol:", sum(keep)),
    paste("mapping_fraction:", sprintf("%.4f", mean(keep))),
    paste("output_symbols:", nrow(tpm_by_symbol)),
    paste("samples:", ncol(counts)),
    "method: counts divided by union exon length, scaled to one million RPK per sample",
    "limitation: gene-level TPM fallback; complete Salmon transcript TPM is preferable",
    paste("R_version:", R.version.string),
    paste("GenomicRanges_version:", as.character(packageVersion("GenomicRanges")))
  ),
  paste0(out_prefix, ".method_notes.txt")
)
