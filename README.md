# genome-PDAC2026

Reusable, patient-data-safe pipeline code and methods documentation for the PDAC2026 African pancreatic cancer multi-omics analysis.

This repository is intended to document how the analyses were performed, not to distribute private data or patient-level results. Raw FASTQ/BAM/CRAM/VCF/MAF/CNV/SV files, expression matrices, reports, figures and real result tables are deliberately excluded.

## Website

The GitHub Pages documentation lives in `docs/` and is designed to be published from the `main` branch using GitHub Pages.

Main pages:

- `docs/index.md`: project overview and pipeline index
- `docs/pipelines.md`: exhaustive pipeline and method map

## Code Layout

- `pipelines/rnaseq/`: nf-core/rnaseq launcher/configuration and differential-expression templates
- `pipelines/immune_infiltration/`: paired tumour-normal immune/stromal score comparison template
- `pipelines/driver_mutation_review/`: safe example showing how TP53 mutation evidence was extracted from VEP-annotated Sarek somatic calls
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

## RNA-seq Downstream Analysis

RNA-seq downstream analysis has two different layers that should not be mixed up. First, expression matrices are used to calculate immune, stromal, EMT, CAF and pathway-related scores. Those scores are then integrated to assign broad tumour phenotype groups. Only after those groups exist does the limma-voom phenotype comparison test which genes differ between groups.

The reproducible tumour-normal and phenotype-comparison code is in:

- `pipelines/rnaseq/differential_expression_limma_template.R`
- `pipelines/rnaseq/phenotype_group_comparison_limma_template.R`
- `pipelines/rnaseq/run_standard_de_from_star_readspergene.R`
- `pipelines/immune_infiltration/paired_tumour_normal_immune_comparison.R`
- `docs/pipelines.md`

For reporting, the defensible phrasing is: RNA-seq quantification was performed with nf-core/rnaseq; downstream tumour-normal and phenotype comparisons were performed with limma-voom/edgeR-based scripts in this repository, with DESeq2 used as a sensitivity analysis where explicitly stated. Historical fallback outputs should be labelled exploratory and not described as DESeq2.

Important distinction for students: the phenotype-group limma script does not calculate immune, stromal or EMT scores. Those scores are upstream inputs used to assign tumours to groups such as `StromalHigh_EMTHigh_ImmuneLow` or `ImmuneHigh_StromalLow`. The limma script only compares expression between groups after that assignment exists in the metadata.

Tumour-normal immune/stromal comparisons are a separate analysis. The immune script takes deconvolution score matrices from ESTIMATE, MCP-counter, CIBERSORT LM22, EPIC, xCell or quanTIseq, matches tumour and normal samples by patient, calculates tumour-minus-normal deltas, and performs paired Wilcoxon tests with FDR correction across features within each method.

## Safety

This repository should remain code-only. Before committing, run:

```bash
git status --short
find . -type f \( -name "*.tsv" -o -name "*.csv" -o -name "*.vcf*" -o -name "*.maf" -o -name "*.bam" -o -name "*.cram" -o -name "*.fastq*" \)
```

No patient-data files should be staged.
