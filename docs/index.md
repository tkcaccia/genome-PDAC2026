# genome-PDAC2026

## Project Message in Plain Language

This project studied pancreatic cancers from 14 African patients using several layers of information: DNA mutations, germline-relevant variants, copy-number changes, structural variants, RNA expression, RNA fusion calls, immune/stromal signals and integrated patient-level summaries.

The central analysis message is cautious but important: this small African PDAC cohort contains patient-level noncanonical molecular findings that should be interpreted in the context of published PDAC literature, not as a formal ancestry comparison. In particular, one tumour emerged as a strong computational candidate for a KRAS-wild-type, MSI-high/MMR-deficient pancreatic cancer. That finding is biologically interesting because most pancreatic ductal adenocarcinomas are KRAS-driven, whereas MSI/MMR-deficient tumours can reflect a different DNA-repair biology. The result is not presented as clinically confirmed; it should be validated with orthogonal pathology or molecular testing.

The second message is that the tumour microenvironment varied across patients. Some tumours looked stromal-rich, fibroblast/CAF-rich and EMT-high, while others looked more immune-high and stromal-low. These expression phenotypes are exploratory because the cohort is small and immune-deconvolution tools use different score scales and assumptions, but they provide a useful genotype-phenotype framework for future validation.

The analysis should not overclaim broad cohort-wide hypermutation. Earlier broader mutation-burden signals became weaker after stricter somatic filtering. The stronger conclusion is one high-priority MSI/MMR-deficient candidate and several KRAS-wild-type or non-clean-KRAS tumours needing focused validation.

## Repository Purpose

This GitHub Pages site documents the pipelines and code used for the analysis in a patient-data-safe way. It includes reusable scripts, configuration files and detailed method notes, but excludes all patient-level result tables, raw sequencing data, clinical metadata, figures, reports, manuscript drafts and reviewer comments.

## Analysis Components

- Bulk RNA-seq quantification with `nf-core/rnaseq 3.24.0`
- Tumour-normal and germline WES analysis with `nf-core/sarek 3.8.1`
- Copy-number and structural-variant review from ASCAT and Manta outputs
- RNA fusion review with `nf-core/rnafusion 4.1.0`
- COSMIC somatic annotation using local COSMIC database tables
- SBS/DBS/ID, CNV and SV signature workflows using SigProfilerAssignment
- Immune/stromal deconvolution using ESTIMATE, MCP-counter, CIBERSORT, EPIC, xCell and quanTIseq
- KRAS-focused algorithmic review
- MSI/MMR/hypermutation review using MSIsensor-pro and strict somatic filtering
- Integrated patient-level genotype/phenotype interpretation

## Documentation Pages

- [Pipeline details](pipelines.md)
- [Differential expression clarification](differential_expression.md)
