#!/usr/bin/env Rscript

cran_packages <- c("data.table", "ggplot2", "matrixStats", "patchwork", "circlize", "remotes")
bioconductor_packages <- c(
  "AnnotationDbi", "Biobase", "DESeq2", "edgeR", "fgsea", "GEOquery",
  "GSVA", "hgu133plus2.db", "hugene10sttranscriptcluster.db", "limma",
  "GenomicRanges", "ComplexHeatmap", "MCPcounter"
)

missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
missing_bioc <- setdiff(bioconductor_packages, rownames(installed.packages()))
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

if (!requireNamespace("immunedeconv", quietly = TRUE)) {
  remotes::install_github("omnideconv/immunedeconv")
}
