#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(
    "Usage: run_standard_de_from_star_readspergene.R <star_analysis_dir> <outdir> [phenotype_assignment.tsv|-] [annotation.gtf.gz]\n",
    "star_analysis_dir should contain sample subfolders with *ReadsPerGene.out.tab files.",
    call. = FALSE
  )
}

star_dir <- normalizePath(args[[1]], mustWork = TRUE)
outdir <- args[[2]]
phenotype_file <- if (length(args) >= 3 && args[[3]] != "-") args[[3]] else NA_character_
gtf_file <- if (length(args) >= 4) normalizePath(args[[4]], mustWork = TRUE) else NA_character_
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_star_counts <- function(path) {
  # STAR writes four count columns. The manuscript used the unstranded column
  # because the project-level count matrix was built from unstranded counts.
  dt <- fread(path, header = FALSE, col.names = c("gene_id", "unstranded", "sense", "antisense"))
  dt <- dt[!grepl("^N_", gene_id)]
  dt[, .(gene_id, count = as.integer(round(unstranded)))]
}

files <- list.files(star_dir, pattern = "ReadsPerGene\\.out\\.tab$", recursive = TRUE, full.names = TRUE)
if (length(files) < 4) {
  stop("Too few STAR ReadsPerGene files found in: ", star_dir, call. = FALSE)
}

sample_from_file <- function(path) {
  sample_dir <- basename(dirname(path))
  if (grepl("^[0-9]+[NT]$", sample_dir)) {
    return(paste0(sample_dir, "_RNA"))
  }
  base <- basename(path)
  sub("^PC0*([0-9]+[NT]).*$", "\\1_RNA", base)
}

normalize_phenotype_group <- function(x) {
  x <- as.character(x)
  key <- tolower(gsub("[^A-Za-z0-9]", "", x))
  x[key %in% c("immunehighstromallow", "immunehighstromalow")] <- "ImmuneHigh_StromalLow"
  x[key %in% c(
    "stromalhighemthighimmunelow",
    "stromalemthighimmunelow",
    "stromahighemthighimmunelow",
    "stromalhighenthighimmunelow",
    "stromahighenthighimmunelow"
  )] <- "StromalHigh_EMTHigh_ImmuneLow"
  x[key == "intermediate"] <- "Intermediate"
  x
}

read_gene_symbol_map <- function(path) {
  if (is.na(path)) return(data.table(gene_id = character(), gene_symbol = character()))
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
  gtf[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attributes)]
  gtf[, gene_symbol := sub('.*gene_name "([^"]+)".*', "\\1", attributes)]
  gtf[, gene_id := sub("\\..*$", "", gene_id)]
  unique(gtf[gene_id != attributes & gene_symbol != attributes, .(gene_id, gene_symbol)], by = "gene_id")
}

add_gene_symbols <- function(result, gene_map) {
  result[, gene_symbol := gene_map$gene_symbol[match(sub("\\..*$", "", gene_id), gene_map$gene_id)]]
  setcolorder(result, c("gene_id", "gene_symbol", setdiff(names(result), c("gene_id", "gene_symbol"))))
  result
}

count_list <- lapply(files, read_star_counts)
sample_ids <- vapply(files, sample_from_file, character(1))
names(count_list) <- sample_ids
count_list <- Map(function(dt, sample_id) {
  setnames(dt, "count", sample_id)
  dt
}, count_list, sample_ids)

count_dt <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_list)
for (j in setdiff(names(count_dt), "gene_id")) {
  set(count_dt, which(is.na(count_dt[[j]])), j, 0L)
}
count_dt <- count_dt[order(gene_id)]

counts <- as.matrix(count_dt[, -1, with = FALSE])
rownames(counts) <- count_dt$gene_id
storage.mode(counts) <- "integer"

metadata <- data.table(sample_id = colnames(counts))
metadata[, patient_id := sub("^([0-9]+)[NT]_RNA$", "\\1", sample_id)]
metadata[, condition := ifelse(grepl("T_RNA$", sample_id), "Tumour", "Normal")]
metadata[, condition := factor(condition, levels = c("Normal", "Tumour"))]
metadata[, patient_id := factor(patient_id)]
setorder(metadata, patient_id, condition)
counts <- counts[, metadata$sample_id, drop = FALSE]

if (!is.na(phenotype_file) && file.exists(phenotype_file)) {
  phenotype <- fread(phenotype_file)
  phenotype[, patient_id := as.character(patient_id)]
  phenotype[, phenotype_group := normalize_phenotype_group(phenotype_group)]
  metadata[, patient_id_chr := as.character(patient_id)]
  metadata <- merge(
    metadata,
    phenotype[, .(patient_id, phenotype_group)],
    by.x = "patient_id_chr",
    by.y = "patient_id",
    all.x = TRUE,
    sort = FALSE
  )
  metadata[, patient_id_chr := NULL]
} else {
  metadata[, phenotype_group := NA_character_]
}

fwrite(count_dt, file.path(outdir, "star_unstranded_gene_counts_matrix.tsv"), sep = "\t")
fwrite(metadata, file.path(outdir, "rnaseq_metadata_for_standard_DE.tsv"), sep = "\t")

keep <- filterByExpr(counts, group = metadata$condition)
filtered_counts <- counts[keep, , drop = FALSE]
gene_map <- read_gene_symbol_map(gtf_file)

# Paired tumour-normal model used by all three count-analysis frameworks:
#   ~ patient_id + condition
# The patient term absorbs baseline differences between matched individuals.
# The conditionTumour coefficient then estimates tumour versus matched normal.
dge <- DGEList(counts = filtered_counts)
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ patient_id + condition, data = metadata)
voom_fit <- voom(dge, design, plot = FALSE)
fit <- eBayes(lmFit(voom_fit, design))

limma_res <- topTable(fit, coef = "conditionTumour", number = Inf, sort.by = "P")
limma_res <- data.table(gene_id = rownames(limma_res), limma_res)
limma_res <- add_gene_symbols(limma_res, gene_map)
limma_res[, significant_FDR_0_05_logFC_1 := adj.P.Val < 0.05 & abs(logFC) > 1]
fwrite(limma_res, file.path(outdir, "DE_tumour_vs_normal_paired_limma_voom.tsv"), sep = "\t")
fwrite(data.table(gene_id = rownames(voom_fit$E), voom_fit$E), file.path(outdir, "limma_voom_logCPM.tsv"), sep = "\t")

# edgeR quasi-likelihood provides an independent negative-binomial count model
# with robust dispersion estimation. It uses exactly the same filtered counts
# and paired design as DESeq2 and limma-voom.
dge_ql <- estimateDisp(dge, design, robust = TRUE)
edge_fit <- glmQLFit(dge_ql, design, robust = TRUE)
edge_test <- glmQLFTest(edge_fit, coef = "conditionTumour")
edge_res <- as.data.table(topTags(edge_test, n = Inf, sort.by = "PValue")$table, keep.rownames = "gene_id")
edge_res <- add_gene_symbols(edge_res, gene_map)
edge_res[, significant_FDR_0_05_logFC_1 := FDR < 0.05 & abs(logFC) > 1]
fwrite(edge_res, file.path(outdir, "DE_tumour_vs_normal_paired_edgeR_QL.tsv"), sep = "\t")

# DESeq2 is the prespecified primary count model. edgeR quasi-likelihood and
# limma-voom above are sensitivity analyses fitted to the same filtered
# integer-count matrix and paired design. The result should only be described
# as DESeq2 when this section completes successfully.
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = as.data.frame(metadata),
  design = ~ patient_id + condition
)
dds <- DESeq(dds, quiet = TRUE)
deseq_res <- as.data.table(results(dds, contrast = c("condition", "Tumour", "Normal")), keep.rownames = "gene_id")
setnames(deseq_res, old = c("log2FoldChange", "pvalue", "padj"), new = c("DESeq2_log2FoldChange", "DESeq2_pvalue", "DESeq2_padj"))
deseq_res <- add_gene_symbols(deseq_res, gene_map)
deseq_res[, significant_FDR_0_05_logFC_1 := DESeq2_padj < 0.05 & abs(DESeq2_log2FoldChange) > 1]
fwrite(deseq_res, file.path(outdir, "DE_tumour_vs_normal_paired_DESeq2.tsv"), sep = "\t")

# A conservative consensus table makes the method sensitivity explicit. A
# gene is called consensus-significant only when both negative-binomial count
# models pass FDR and effect-size thresholds with the same effect direction.
consensus <- merge(
  deseq_res[, .(
    gene_id,
    gene_symbol,
    DESeq2_log2FoldChange,
    DESeq2_pvalue,
    DESeq2_padj
  )],
  edge_res[, .(
    gene_id,
    edgeR_log2FoldChange = logFC,
    edgeR_PValue = PValue,
    edgeR_FDR = FDR
  )],
  by = "gene_id",
  all = FALSE
)
consensus <- merge(
  consensus,
  limma_res[, .(
    gene_id,
    limma_voom_log2FoldChange = logFC,
    limma_voom_PValue = P.Value,
    limma_voom_FDR = adj.P.Val
  )],
  by = "gene_id",
  all = FALSE
)
consensus[, effect_direction_concordant :=
  sign(DESeq2_log2FoldChange) == sign(edgeR_log2FoldChange) &
  sign(DESeq2_log2FoldChange) == sign(limma_voom_log2FoldChange)]
consensus[, significant_DESeq2_edgeR_consensus :=
  effect_direction_concordant &
  DESeq2_padj < 0.05 & abs(DESeq2_log2FoldChange) > 1 &
  edgeR_FDR < 0.05 & abs(edgeR_log2FoldChange) > 1]
setorder(consensus, DESeq2_padj, edgeR_FDR)
fwrite(consensus, file.path(outdir, "DE_tumour_vs_normal_paired_count_model_consensus.tsv"), sep = "\t")

norm_counts <- counts(dds, normalized = TRUE)
fwrite(data.table(gene_id = rownames(norm_counts), norm_counts), file.path(outdir, "DESeq2_normalized_counts.tsv"), sep = "\t")

extreme_groups <- c("ImmuneHigh_StromalLow", "StromalHigh_EMTHigh_ImmuneLow")
tumour_meta <- metadata[condition == "Tumour" & phenotype_group %in% extreme_groups]
if (nrow(tumour_meta) >= 4 && length(unique(tumour_meta$phenotype_group)) == 2) {
  # This tumour-only contrast starts after phenotype labels have already been
  # assigned from upstream immune/stromal/EMT score integration. It does not
  # calculate those scores itself.
  tumour_meta[, phenotype_group := factor(phenotype_group, levels = extreme_groups)]
  tumour_counts <- counts[, tumour_meta$sample_id, drop = FALSE]
  tumour_keep <- filterByExpr(tumour_counts, group = tumour_meta$phenotype_group)
  tumour_dge <- DGEList(counts = tumour_counts[tumour_keep, , drop = FALSE])
  tumour_dge <- calcNormFactors(tumour_dge, method = "TMM")
  tumour_design <- model.matrix(~ phenotype_group, data = tumour_meta)
  tumour_voom <- voom(tumour_dge, tumour_design, plot = FALSE)
  tumour_fit <- eBayes(lmFit(tumour_voom, tumour_design))
  tumour_res <- topTable(tumour_fit, coef = "phenotype_groupStromalHigh_EMTHigh_ImmuneLow", number = Inf, sort.by = "P")
  tumour_res <- data.table(gene_id = rownames(tumour_res), tumour_res)
  tumour_res <- add_gene_symbols(tumour_res, gene_map)
  tumour_res[, significant_FDR_0_05_logFC_1 := adj.P.Val < 0.05 & abs(logFC) > 1]
  fwrite(tumour_res, file.path(outdir, "DE_tumour_phenotype_stromal_EMT_high_vs_immune_high_limma_voom.tsv"), sep = "\t")
} else {
  fwrite(
    data.table(reason = "Skipped: phenotype assignment missing or fewer than two extreme groups."),
    file.path(outdir, "DE_tumour_phenotype_stromal_EMT_high_vs_immune_high_SKIPPED.tsv"),
    sep = "\t"
  )
}

summary <- data.table(
  metric = c(
    "star_readspergene_files",
    "samples",
    "patients",
    "input_genes",
    "filtered_genes",
    "limma_significant_FDR_0_05_abs_logFC_1",
    "edger_ql_significant_FDR_0_05_abs_logFC_1",
    "deseq2_significant_FDR_0_05_abs_logFC_1",
    "deseq2_edger_consensus_significant_FDR_0_05_abs_logFC_1"
  ),
  value = c(
    length(files),
    ncol(counts),
    length(unique(metadata$patient_id)),
    nrow(counts),
    nrow(filtered_counts),
    sum(limma_res$significant_FDR_0_05_logFC_1, na.rm = TRUE),
    sum(edge_res$significant_FDR_0_05_logFC_1, na.rm = TRUE),
    sum(deseq_res$significant_FDR_0_05_logFC_1, na.rm = TRUE),
    sum(consensus$significant_DESeq2_edgeR_consensus, na.rm = TRUE)
  )
)
fwrite(summary, file.path(outdir, "standard_DE_summary.tsv"), sep = "\t")

sink(file.path(outdir, "standard_DE_session_info.txt"))
cat("Command:", paste(commandArgs(), collapse = " "), "\n\n")
cat("STAR input directory:", star_dir, "\n")
cat("Output directory:", normalizePath(outdir), "\n")
cat("GENCODE annotation:", if (is.na(gtf_file)) "not supplied" else gtf_file, "\n")
cat("GENCODE annotation MD5:", if (is.na(gtf_file)) "not supplied" else unname(tools::md5sum(gtf_file)), "\n")
cat("Primary model: DESeq2 with design ~ patient_id + condition\n")
cat("Sensitivity models: edgeR quasi-likelihood and limma-voom with the same design\n")
cat("Consensus rule: DESeq2 and edgeR FDR < 0.05, absolute log2FC > 1, concordant direction\n\n")
print(sessionInfo())
sink()

message("Standard DE complete: ", normalizePath(outdir))
