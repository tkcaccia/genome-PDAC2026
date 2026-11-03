#!/usr/bin/env Rscript
# Current collaborator code export. Expression-to-deconvolution scores, paired comparisons and tumour-only Hallmark/TME groups.
# Source: Original Source and Code/scripts/run_salmon_tme_methods_11_pdac.R
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
suppressPackageStartupMessages({library(immunedeconv); library(ggplot2); library(tidyr)})
source(file.path(Sys.getenv("PDAC_CODE_ROOT"), "scripts/shared/runtime.R"))
root <- pdac_path("project_root")
out <- file.path(root, "results/immune_salmon_11_20261012")
input <- file.path(out, "inputs/salmon_gene_symbol_tpm_11_pairs.csv")
tpm <- as.matrix(read.csv(input, row.names = 1, check.names = FALSE))
storage.mode(tpm) <- "double"
metadata <- read.csv(file.path(out, "inputs/sample_metadata.csv"))
patients <- unique(metadata$patient)
samples <- metadata$sample
stopifnot(identical(colnames(tpm), samples), ncol(tpm) == 22, all(is.finite(tpm)), all(tpm >= 0))
log_tpm <- log2(tpm + 1)
reference_dir <- file.path(out, "inputs/references")
dir.create(reference_dir, showWarnings = FALSE)
marker_file <- file.path(reference_dir, "MCPcounter_genes.txt")
marker_url <- "https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"
if (!file.exists(marker_file)) download.file(marker_url, marker_file, quiet = TRUE, mode = "wb")
markers <- read.delim(marker_file, check.names = FALSE, colClasses = "character")
stopifnot(all(c("HUGO symbols", "Cell population") %in% names(markers)))
coverage <- list()
record_coverage <- function(method, feature, genes) {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  data.frame(method = method, feature = feature, signature_genes = length(genes), mapped_genes = sum(genes %in% rownames(tpm)), coverage_percent = 100 * mean(genes %in% rownames(tpm)))
}
coverage[["MCPcounter"]] <- do.call(rbind, lapply(unique(markers[["Cell population"]]), function(cell) record_coverage("MCPcounter", cell, markers[markers[["Cell population"]] == cell, "HUGO symbols"])))
epic_signature <- EPIC::TRef
coverage[["EPIC"]] <- record_coverage("EPIC", "TRef signature genes", epic_signature$sigGenes)
estimate_reference <- readRDS(system.file("extdata/estimate/estimate_genes.rds", package = "immunedeconv", mustWork = TRUE))
estimate_sets <- estimate_reference$signature_geneset
coverage[["ESTIMATE"]] <- do.call(rbind, lapply(rownames(estimate_sets), function(cell) record_coverage("ESTIMATE", cell, as.character(unlist(estimate_sets[cell, -1])))))
xcell_sets <- GSEABase::geneIds(xCell::xCell.data$signatures)
coverage[["xCell"]] <- do.call(rbind, lapply(names(xcell_sets), function(cell) record_coverage("xCell", cell, xcell_sets[[cell]])))
write.csv(do.call(rbind, coverage), file.path(out, "qc/multimethod_signature_coverage.csv"), row.names = FALSE)

gsva_compat <- function(expr, geneSets, method, ssgsea.norm, parallel.sz, ...) {
  stopifnot(method == "ssgsea")
  parameters <- GSVA::ssgseaParam(expr, geneSets, minSize = 1, maxSize = Inf, alpha = 0.25, normalize = ssgsea.norm)
  GSVA::gsva(parameters, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
}
replace_gsva_call <- function(expression) {
  if (!is.call(expression)) return(expression)
  if (identical(expression[[1]], quote(GSVA::gsva))) expression[[1]] <- as.name("gsva_compat")
  for (index in seq_along(expression)[-1]) expression[[index]] <- replace_gsva_call(expression[[index]])
  expression
}
# xCell 1.3 uses the retired GSVA call; adapt only that call in this R session.
raw_function <- get("rawEnrichmentAnalysis", asNamespace("xCell"))
body(raw_function) <- replace_gsva_call(body(raw_function))
adapter_env <- new.env(parent = environment(raw_function))
adapter_env$gsva_compat <- gsva_compat
environment(raw_function) <- adapter_env
assignInNamespace("rawEnrichmentAnalysis", raw_function, ns = "xCell")

results <- list()
status <- list()
for (method in c("estimate", "mcp_counter", "epic", "xcell")) {
  cat("Running", method, "\n")
  result <- tryCatch({
    if (method == "estimate") scores <- immunedeconv::deconvolute_estimate(log_tpm)
    if (method == "mcp_counter") scores <- immunedeconv::deconvolute_mcp_counter(log_tpm, genes = markers)
    if (method == "epic") {
      epic_fit <- EPIC::EPIC(tpm, reference = "TRef")
      write.csv(epic_fit$fit.gof, file.path(out, "qc/epic_goodness_of_fit.csv"), row.names = TRUE)
      scores <- t(epic_fit$cellFractions)
    }
    if (method == "xcell") scores <- immunedeconv::deconvolute_xcell(log_tpm, arrays = FALSE, parallel.sz = 1)
    if (method == "estimate") scores <- scores[rownames(scores) != "TumorPurity", , drop = FALSE]
    scores <- as.matrix(scores[, samples, drop = FALSE])
    stopifnot(ncol(scores) == 22, all(is.finite(scores)), !anyDuplicated(rownames(scores)))
    if (method == "epic") stopifnot(all(scores >= -1e-6), all(scores <= 1 + 1e-6), all(abs(colSums(scores) - 1) < 1e-5))
    scores
  }, error = function(error) error)
  if (inherits(result, "error")) {
    status[[method]] <- data.frame(method = method, status = "failed", detail = conditionMessage(result))
    cat("FAILED", method, conditionMessage(result), "\n")
    next
  }
  results[[method]] <- result
  write.csv(result, file.path(out, "tables", paste0(method, "_native_scores_11_pairs.csv")))
  status[[method]] <- data.frame(method = method, status = "completed", detail = paste(nrow(result), "features; 22 samples"))
  if (method == "epic" && any(epic_fit$fit.gof$convergeCode != 0)) {
    status[[method]]$status <- "completed_with_convergence_warning"
    status[[method]]$detail <- paste("Provisional full-cohort results; nonconverged:", paste(rownames(epic_fit$fit.gof)[epic_fit$fit.gof$convergeCode != 0], collapse = "; "))
  }
}
write.csv(do.call(rbind, status), file.path(out, "qc/multimethod_run_status.csv"), row.names = FALSE)
paired_tables <- list()
for (method in names(results)) {
  scores <- results[[method]]
  paired <- do.call(rbind, lapply(rownames(scores), function(cell) {
    normal <- scores[cell, paste0(patients, "N")]
    tumour <- scores[cell, paste0(patients, "T")]
    delta <- tumour - normal
    pvalue <- if (all(delta == 0)) 1 else wilcox.test(tumour, normal, paired = TRUE, exact = FALSE)$p.value
    data.frame(method = method, feature = cell, n_pairs = 11, normal_mean = mean(normal), tumour_mean = mean(tumour), median_paired_delta = median(delta), p_value = pvalue)
  }))
  paired$FDR <- p.adjust(paired$p_value, method = "BH")
  paired$qc_status <- status[[method]]$status
  write.csv(paired, file.path(out, "tables", paste0(method, "_paired_comparison_11_pdac.csv")), row.names = FALSE)
  paired_tables[[method]] <- paired
  if (method == "epic") {
    failed_samples <- rownames(epic_fit$fit.gof)[epic_fit$fit.gof$convergeCode != 0]
    passed_patients <- setdiff(patients, sub("[NT]$", "", failed_samples))
    sensitivity <- do.call(rbind, lapply(rownames(scores), function(cell) {
      normal <- scores[cell, paste0(passed_patients, "N")]
      tumour <- scores[cell, paste0(passed_patients, "T")]
      delta <- tumour - normal
      pvalue <- if (all(delta == 0)) 1 else wilcox.test(tumour, normal, paired = TRUE, exact = FALSE)$p.value
      data.frame(method = method, feature = cell, n_pairs = length(passed_patients), normal_mean = mean(normal), tumour_mean = mean(tumour), median_paired_delta = median(delta), p_value = pvalue)
    }))
    sensitivity$FDR <- p.adjust(sensitivity$p_value, method = "BH")
    write.csv(sensitivity, file.path(out, "tables/epic_converged_pairs_sensitivity.csv"), row.names = FALSE)
  }
  selected <- if (method == "xcell") head(paired[order(paired$FDR, paired$p_value), "feature"], 12) else rownames(scores)
  long <- data.frame(feature = rownames(scores[selected, , drop = FALSE]), scores[selected, , drop = FALSE], check.names = FALSE)
  long <- pivot_longer(long, cols = -feature, names_to = "sample", values_to = "native_score")
  long$patient <- sub("[NT]$", "", long$sample)
  long$condition <- factor(ifelse(grepl("T$", long$sample), "Tumour", "Normal"), levels = c("Normal", "Tumour"))
  plot <- ggplot(long, aes(condition, native_score, group = patient)) + geom_line(color = "#AAAAAA") +
    geom_point(aes(color = condition), size = 1.6) + facet_wrap(~feature, scales = "free_y", ncol = 4) +
    scale_color_manual(values = c(Normal = "#4477AA", Tumour = "#EE6677")) +
    labs(x = NULL, y = if (method == "epic") "Estimated cell fraction" else "Native method score", title = paste(method, "paired tumour-normal inference: 11 PDAC pairs"), subtitle = if (method == "xcell") "12 lowest-FDR features shown; all scores/tests retained in CSV" else if (method == "epic") status[[method]]$detail else NULL) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(out, "figures", paste0("Figure_TME_paired_", method, "_11_pdac.png")), plot, width = 12, height = if (method == "estimate") 4.5 else 8, dpi = 300)
}
write.csv(do.call(rbind, paired_tables), file.path(out, "tables/TME_multimethod_paired_comparison_11_pdac.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "qc/R_sessionInfo_multimethod.txt"))
packages <- c("immunedeconv", "EPIC", "MCPcounter", "xCell", "GSVA")
write.csv(data.frame(package = packages, version = vapply(packages, function(package) as.character(packageVersion(package)), character(1))), file.path(out, "qc/multimethod_package_versions.csv"), row.names = FALSE)
write.csv(data.frame(file = marker_file, source = marker_url, md5 = unname(tools::md5sum(marker_file))), file.path(out, "qc/mcp_reference_checksum.csv"), row.names = FALSE)
writeLines(c("All methods use the same summed gene-symbol Salmon TPM matrix from the quanTIseq input audit.", "EPIC: unlogged TPM; TRef (tumor=TRUE) for both tumour and adjacent-normal samples to ensure a common basis; mRNA scaling enabled.", "MCP-counter: log2(TPM+1); official HTTPS marker list cached; scores are not fractions and are comparable between samples for a given cell population.", "xCell: log2(TPM+1), rnaseq=TRUE, one worker; within-sample ranks, unnormalised ssGSEA, transformation and spillover retained. Session-local adapter uses GSVA ssgseaParam(alpha=0.25,normalize=FALSE), because installed xCell uses retired GSVA syntax.", "ESTIMATE: log2(TPM+1), native rank scoring; ImmuneScore, StromalScore and ESTIMATEScore retained. Microarray-calibrated TumorPurity is excluded from inference.", "Each method tested separately using two-sided paired Wilcoxon (exact=FALSE) and BH within method. Do not combine native scores/fractions across methods.", "These estimates do not imply cell counts or histological verification.", "GSVA/ssGSEA for 19 custom programmes and the 12-feature TME phenotype classification are not performed: exact prespecified definitions unavailable."), file.path(out, "MULTIMETHOD_METHODS_AND_LIMITATIONS.txt"))
print(do.call(rbind, status))
for (method in names(paired_tables)) cat(method, "FDR<0.05:", sum(paired_tables[[method]]$FDR < 0.05), "/", nrow(paired_tables[[method]]), "\n")
if (any(vapply(status, function(row) row$status == "failed", logical(1)))) stop("One or more TME methods failed; inspect qc/multimethod_run_status.csv")