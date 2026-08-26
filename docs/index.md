# genome-PDAC2026

## Methods Scope In Plain Language

This repository explains how pancreatic-cancer DNA and RNA sequencing outputs can be processed into mutation, copy-number, structural-variant, fusion, expression, pathway and tumour-microenvironment summaries. It focuses on reusable code and the order in which inputs move through the analysis. Patient-level findings, real result tables, clinical metadata, figures and manuscript text are intentionally excluded from this public site.

The RNA workflow starts from gene expression, calculates immune/stromal and biological-programme scores, combines those scores into exploratory tumour phenotypes, and only then tests gene-expression differences between phenotype groups. The DNA workflow independently produces germline/somatic variants, copy-number and structural-variant evidence; those variables are integrated with RNA phenotypes only after both upstream analyses are complete.

## Repository Purpose

This GitHub Pages site documents the pipelines and code used for the analysis in a patient-data-safe way. It includes reusable scripts, configuration files and detailed method notes, but excludes all patient-level result tables, raw sequencing data, clinical metadata, figures, reports, manuscript drafts and reviewer comments.

## Analysis Components

- Bulk RNA-seq quantification with `nf-core/rnaseq 3.24.0`
- Tumour-normal and germline WES analysis with `nf-core/sarek 3.8.1`
- Copy-number and structural-variant review from ASCAT and Manta outputs
- RNA fusion calling with pinned `nf-core/rnafusion 4.1.0`, four configured
  methods and a bounded sequential workstation queue
- COSMIC somatic annotation using local COSMIC database tables
- SBS/DBS/ID, CNV and SV signature workflows using SigProfilerAssignment
- Immune/stromal deconvolution using ESTIMATE, MCP-counter, CIBERSORT, EPIC, xCell and quanTIseq
- KRAS-focused computational evidence review
- MSI/MMR/hypermutation review using MSIsensor-pro and strict somatic filtering
- Integrated patient-level genotype/phenotype interpretation

For students starting with a normalized RNA expression matrix, the complete expression-to-phenotype walkthrough is in the existing [pipeline details](pipelines.md) and the READMEs under `pipelines/immune_infiltration`, `pipelines/pathway_scoring` and `pipelines/phenotype_assignment`. The pathway README also explains that a GMT is a separate gene-set definition file rather than a value calculated from patient expression.

## Audit Boundary

The current code repository does not imply that all historical primary outputs remain recoverable. The August 2026 storage audit re-ran and verified the RNA downstream layer. It found only recovered patient-summary fields for several genomic analyses; primary Sarek, historical rnafusion and SigProfiler output trees were absent from the connected drives. Because the RNA FASTQs remain available, this repository now documents how to regenerate RNA-fusion outputs safely. A regenerated sample is complete only after its execution report, trace and requested caller outputs pass the documented checks. The absent WES-derived primary evidence still requires source recovery before event-level re-audit.

## Documentation Pages

- [Pipeline details](pipelines.md)
