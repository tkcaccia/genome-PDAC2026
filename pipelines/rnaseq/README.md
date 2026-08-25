# RNA-seq

Reusable RNA expression workflow and audited downstream analysis based on
`nf-core/rnaseq 3.24.0` outputs.

Files:

- `run_rnaseq.sh`: launcher wrapper
- `rnaseq.config`: Nextflow process/resource configuration
- `launch_example.sh`: reusable example for launching or resuming nf-core/rnaseq
- `_common.sh` and `validate_samplesheet.py`: copied shared helpers required by the wrapper
- `run_standard_de_from_star_readspergene.R`: constructs STAR gene-count matrices, runs primary paired DESeq2, and runs paired edgeR quasi-likelihood and limma-voom sensitivity analyses
- `phenotype_group_comparison_limma_template.R`: compares expression between already assigned tumour phenotype groups; it does not calculate immune/stromal/EMT scores
- `de_diagnostics_sensitivity.R`: checks sample QC, pairing, PCA, limma-vs-DESeq2 effect-size concordance, significant-gene overlap and discordant genes
- `make_gene_tpm_from_counts.R`: creates a documented gene-level TPM fallback from gene counts and union exon lengths when complete Salmon TPM is unavailable

## What Each RNA-seq Script Answers

| Script | Biological question | Main required inputs | Main outputs |
| --- | --- | --- | --- |
| `run_rnaseq.sh` | Process RNA-seq reads and create expression outputs | nf-core/rnaseq samplesheet and references | STAR/Salmon expression and QC outputs |
| `run_standard_de_from_star_readspergene.R` | Which genes differ between matched tumour and normal samples? | STAR output directory; optional phenotype assignment table | Count matrix, metadata, primary paired DESeq2, edgeR and limma-voom tables, consensus diagnostics and normalized/log expression |
| `de_three_method_diagnostics.R` | Do DESeq2, edgeR and limma-voom agree in effect size and significance? | The three method-specific result tables | Effect-size correlations, significance overlap and method-sensitivity summary |
| `make_gene_tpm_from_counts.R` | Can immune fraction methods be rerun when complete Salmon TPM is unavailable? | Gene-count matrix plus the matching GTF | Gene-symbol TPM, union exon lengths and checksummed method notes |
| `phenotype_group_comparison_limma_template.R` | Among tumours only, which genes differ between pre-defined phenotype groups? | Count matrix plus metadata containing `phenotype_group` | Tumour-normal limma table, tumour phenotype-group limma table or skipped-file explanation |

The immune-infiltration calculations are not performed in this RNA-seq folder.
Starting from TPM-like expression, use
`../immune_infiltration/run_immune_stromal_scores_from_expression.R`; then use
`../immune_infiltration/paired_tumour_normal_immune_comparison.R` for the
matched tumour-normal statistics.

## Gene-Level TPM Fallback

Complete Salmon transcript TPM is preferred for EPIC, quanTIseq and CIBERSORT. If only gene counts remain, generate a documented gene-level TPM sensitivity input using union exon lengths from the same annotation release:

```bash
Rscript make_gene_tpm_from_counts.R \
  results/star_unstranded_gene_counts_matrix.tsv \
  references/gencode.v46.primary_assembly.annotation.gtf.gz \
  results/rna_expression
```

The output is suitable for sensitivity analysis, not a claim that gene-level TPM is numerically identical to Salmon transcript-level quantification.

## Important: Scoring Is Upstream of Group Comparison

The tumour phenotype comparison uses labels such as `StromalHigh_EMTHigh_ImmuneLow` and `ImmuneHigh_StromalLow`. These labels come from upstream immune/stromal/EMT scoring and integrated interpretation. The limma script only compares expression between groups after the labels already exist.

Where the labels come from:

- The starting point is the normalized RNA-seq expression matrix generated from the nf-core/rnaseq/STAR gene-count outputs.
- Immune/stromal context was estimated with tools such as ESTIMATE, MCP-counter, CIBERSORT and immunedeconv methods where available.
- CAF, EMT, hypoxia, angiogenesis and pathway programmes were scored from curated gene sets using the expression matrix. A patient-data-safe GSVA/ssGSEA implementation is provided in `../pathway_scoring/run_gsva_ssgsea_programme_scores.R`.
- The GMT gene-set input is `../../templates/pdac_programme_gene_sets_example.gmt`. It is a separate curated gene-list file, not an output generated from the expression matrix.
- `../phenotype_assignment/assemble_tme_score_table.py` merged the method-specific matrices by RNA sample ID, and `../phenotype_assignment/assign_tme_phenotype_groups.py` standardized selected features and calculated the labels.
- Tumours with high stromal/CAF/EMT signal and relatively low immune signal were labelled `StromalHigh_EMTHigh_ImmuneLow`; tumours with stronger immune signal and low stromal signal were labelled `ImmuneHigh_StromalLow`; unclear cases were labelled `Intermediate`.

The labels are therefore derived phenotype summaries used for an exploratory tumour-only comparison. They are not raw percentages, and they are not calculated inside the limma script itself.

In other words:

```text
expression matrix -> immune/stromal/EMT scores -> phenotype_group labels -> limma group comparison
```

If a metadata file has no `phenotype_group` column, the phenotype-group comparison cannot be run.

A safe example phenotype-assignment file is available at:

```text
templates/phenotype_assignment_template.tsv
```
