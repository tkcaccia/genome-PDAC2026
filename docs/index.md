---
layout: default
title: genome PDAC2026
---

# genome PDAC2026

This website documents the PDAC2026 pancreatic cancer multi-omics analyses so
readers can trace each method from its required inputs to its outputs and the
figures or tables it supports.

The public repository contains source code, configuration templates and method
documentation only. Patient data, results, figures and manuscript files remain
private.

## Analytical Scope

| Analysis | Purpose |
| --- | --- |
| Whole-exome sequencing | Summarise retained somatic variants, KRAS evidence and computational microsatellite instability |
| RNA sequencing | Compare matched tumour-normal gene expression using paired DESeq2 and characterise expression programmes |
| Immune and stromal analysis | Estimate method-specific microenvironment features, compare matched tissues and describe tumour expression phenotypes |
| Proteomics | Assemble quantified protein abundances and document protein-level inputs |
| Multi-omics integration | Integrate genomic, transcriptomic and protein features using MOFA |
| Cohort comparisons | Describe TCGA PAAD context and analyse periampullary cases separately |

## Start Here

- [Figure and table index](figures-and-tables.md)
- [Analysis workflow and scoring definitions](workflow.md)
- [Private input requirements](inputs.md)
- [Reproducibility and implementation coverage](source-coverage.md)
- [Code files and checksums](source_manifest.json)
- [Browse the code on GitHub](https://github.com/tkcaccia/genome-PDAC2026)

Use the figure/table index to locate a script, the workflow page to understand
the calculation, and the input guide to identify the required file format.
Implementation coverage is stated explicitly: rendering an existing figure is
distinguished from reproducing the analysis behind it.
