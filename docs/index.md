# genome-PDAC2026

## Methods Scope In Plain Language

This repository explains how pancreatic-cancer DNA and RNA sequencing outputs can be processed into mutation, copy-number, structural-variant, fusion, expression, pathway and tumour-microenvironment summaries. It focuses on reusable code and the order in which inputs move through the analysis. Patient-level findings, real result tables, clinical metadata, figures and manuscript text are intentionally excluded from this public site.

The RNA workflow starts from gene expression, calculates immune/stromal and biological-programme scores, combines those scores into exploratory tumour phenotypes, and only then tests gene-expression differences between phenotype groups. The DNA workflow independently produces germline/somatic variants, copy-number and structural-variant evidence; those variables are integrated with RNA phenotypes only after both upstream analyses are complete.

## Repository Purpose

This GitHub Pages site documents the pipelines and code used for the analysis in a patient-data-safe way. It includes reusable scripts, configuration files and detailed method notes, but excludes all patient-level result tables, raw sequencing data, clinical metadata, figures, reports, manuscript drafts and reviewer comments.

## Analysis Components

- Bulk RNA-seq quantification with `nf-core/rnaseq 3.24.0`
- Tumour-normal and germline WES analysis with `nf-core/sarek 3.8.1`
- On-target WES depth and SAGE/PURPLE tumour-normal call-quality auditing
- Copy-number and structural-variant review from ASCAT and Manta outputs
- RNA fusion calling with pinned `nf-core/rnafusion 4.1.0`, four configured
  methods and a bounded sequential workstation queue
- COSMIC somatic annotation using local COSMIC database tables
- SBS/DBS/ID, CNV and SV signature workflows using SigProfilerAssignment
- Immune/stromal deconvolution using ESTIMATE, MCP-counter, EPIC, xCell and
  quanTIseq, with optional licensed CIBERSORT support documented but excluded
  from the audited rerun
- KRAS-focused computational evidence review
- MSI/MMR/hypermutation review using MSIsensor-pro and strict somatic filtering
- Integrated patient-level genotype/phenotype interpretation
- Direct transcriptomic validation against independently processed public paired cohorts, with exact-GMT pathway comparison
- Pancreatic-head-only sensitivity checks for expression and tumour-microenvironment effects

For students starting with a normalized RNA expression matrix, the complete expression-to-phenotype walkthrough is in the existing [pipeline details](pipelines.md) and the READMEs under `pipelines/immune_infiltration`, `pipelines/pathway_scoring` and `pipelines/phenotype_assignment`. The pathway README also explains that a GMT is a separate gene-set definition file rather than a value calculated from patient expression.

## Audit Boundary

The current code repository does not imply that all historical primary outputs remain recoverable. The August 2026 storage audit re-ran and verified the RNA downstream layer. In September 2026, the provider DRAGEN WES enrichment package was recovered for all 28 specimens and Oncoanalyser SAGE/PURPLE outputs were recovered for 13 of 14 pairs. The WES audit therefore now reports capture-target depth, matched callable territory and caller-filtered evidence directly from primary QC/call files. One pair still lacks its final published somatic VCF because the historical Oncoanalyser run failed during downstream LINX plotting. Historical rnafusion, detailed Sarek ASCAT/Manta and SigProfiler trees remain incomplete unless separately documented.

## Documentation Pages

- [Pipeline details](pipelines.md)
