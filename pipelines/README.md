# PDAC2026 Pipeline Index

This folder contains reusable, patient-data-safe code for the pancreatic and
periampullary multi-omics study. It preserves scientific methods, not private
inputs, result tables, machine administration or credentials.

## Current Audit Status

The August 2026 audit verified the following analyses from current source
files:

- paired tumour-normal RNA differential expression with DESeq2, edgeR and
  limma-voom;
- documented gene-level TPM construction from STAR counts and GENCODE v46;
- GSVA and ssGSEA scoring of 19 curated programmes;
- ESTIMATE, MCP-counter, EPIC, xCell and quanTIseq;
- paired tumour-normal immune/stromal statistics;
- cohort-relative tumour-microenvironment phenotype assignment;
- anonymized manuscript plots and technical-validation summaries.

A September 2026 recovery found provider DRAGEN enrichment outputs for all 28
WES specimens and SAGE/PURPLE outputs for 13 of 14 matched pairs. The
`wes_qc_calls` workflow now audits capture-target depth, paired callable
territory and caller-filtered somatic evidence. One Oncoanalyser run completed
variant calling but failed during LINX plotting before final files were
published; that pair remains explicitly missing rather than being reconstructed
by subtracting separate germline VCFs. Detailed historical ASCAT/Manta,
rnafusion and SigProfiler trees remain incomplete unless separately documented.

## Analysis Folders

| Folder | Purpose | Main input |
| --- | --- | --- |
| `rnaseq` | nf-core/rnaseq launch and paired DE | RNA FASTQ or STAR gene counts |
| `pathway_scoring` | GSVA/ssGSEA programme activity | normalized/log expression plus GMT |
| `immune_infiltration` | immune/stromal inference and paired tests | linear TPM-like expression |
| `phenotype_assignment` | RNA score integration and TME labels | method-specific score tables |
| `driver_mutation_review` | safe TP53 filtering example | VEP-annotated somatic variants |
| `sarek_germline` | nf-core/sarek germline WES | normal WES samplesheet |
| `sarek_tumor_normal` | nf-core/sarek tumour-normal WES | matched WES samplesheet |
| `wes_qc_calls` | target depth, paired callability and genomic-call audit | DRAGEN enrichment and SAGE/PURPLE outputs |
| `sarek_sv_cna` | WES-compatible ASCAT/Manta summaries and circos plots | Sarek CNA/SV outputs |
| `rnafusion` | nf-core/rnafusion wrapper | RNA FASTQ samplesheet |
| `cosmic_annotation` | licensed local COSMIC annotation | annotated somatic calls |
| `cosmic_signatures` | SBS/DBS/ID assignment | PASS somatic VCFs |
| `cosmic_sv_cna_signatures` | exploratory CNV48/SV32 assignment | ASCAT/Manta outputs |
| `integrated_analysis` | anonymized multi-omics figures | harmonized safe summaries |
| `validation` | KRAS/MSI/TMB summary figure | harmonized safe summaries |

## Interpretation Boundaries

- The cohort contains pancreatic-head and ampullary/periampullary cancers and
  should not be described as a histologically uniform classic-PDAC cohort.
- WES-derived CNA and SV calls have lower and less uniform structural
  resolution than whole-genome sequencing.
- RNA deconvolution outputs use method-specific scales and are not all
  percentages.
- The 3/3/8 TME labels are cohort-relative exploratory groups created from RNA
  scores, not clinical diagnoses or independent validation.
- CIBERSORT requires its licensed script and LM22 matrix. It is optional and
  was excluded from the audited rerun because those files were unavailable.
- Signature, fusion and detailed event-level claims require restoration of the
  primary output files and execution logs.

All runtime paths must be supplied through command-line arguments,
configuration or environment variables. No retired workstation path is a
valid default.
