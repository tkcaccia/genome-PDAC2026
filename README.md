# genome-PDAC2026

Reusable, patient-data-safe pipeline code and methods documentation for the PDAC2026 African pancreatic cancer multi-omics manuscript.

This repository is intended to document how the analyses were performed, not to distribute private data or patient-level results. Raw FASTQ/BAM/CRAM/VCF/MAF/CNV/SV files, expression matrices, reports, figures and real result tables are deliberately excluded.

## Website

The GitHub Pages documentation lives in `docs/` and is designed to be published from the `main` branch using GitHub Pages.

Main pages:

- `docs/index.md`: manuscript story and pipeline overview
- `docs/pipelines.md`: exhaustive pipeline and method map
- `docs/differential_expression.md`: student-facing differential-expression clarification and recommended analysis code
- `docs/manuscript_reviewer_comments.md`: reviewer-style comments for Journal of Translational Medicine

## Code Layout

- `pipelines/rnaseq/`: nf-core/rnaseq launcher/configuration and differential-expression templates
- `pipelines/sarek_tumor_normal/`: nf-core/sarek tumour-normal WES calling
- `pipelines/sarek_germline/`: nf-core/sarek germline WES calling
- `pipelines/sarek_sv_cna/`: SV/CNA recovery, interchromosomal translocation summaries and circos-style plotting
- `pipelines/rnafusion/`: nf-core/rnafusion launcher/configuration
- `pipelines/cosmic_signatures/`: SigProfilerAssignment SBS/DBS/ID signature assignment
- `pipelines/cosmic_sv_cna_signatures/`: SigProfilerAssignment CNV/SV signature assignment
- `pipelines/cosmic_annotation/`: local COSMIC database annotation of somatic calls
- `pipelines/cohort_autodraft/`: samplesheet/autodraft helper code
- `pipelines/shared_runtime/`: shared shell and validation helpers
- `pipelines/monitoring_and_recovery/`: remote watchdog/recovery scripts

## Differential Expression Clarification

One historical downstream file name used during exploratory analysis implied `DESeq2_or_limma`, but one fallback calculation used custom Python log2CPM and paired statistics rather than DESeq2, edgeR or limma. This repository addresses that concern by documenting the distinction and providing a reproducible `limma-voom` paired workflow in:

- `pipelines/rnaseq/differential_expression_limma_template.R`
- `pipelines/rnaseq/phenotype_group_comparison_limma_template.R`
- `docs/differential_expression.md`

For manuscript language, the defensible phrasing is: RNA-seq quantification was performed with nf-core/rnaseq; downstream tumour-normal and phenotype comparisons should be reported as limma-voom/edgeR-based when generated with the scripts in this repository. Historical fallback outputs should be labelled exploratory and not described as DESeq2.

## Safety

This repository should remain code-only. Before committing, run:

```bash
git status --short
find . -type f \( -name "*.tsv" -o -name "*.csv" -o -name "*.vcf*" -o -name "*.maf" -o -name "*.bam" -o -name "*.cram" -o -name "*.fastq*" \)
```

No patient-data files should be staged.
