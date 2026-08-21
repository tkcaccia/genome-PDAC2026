# GSVA/ssGSEA Pathway and Programme Scoring

This folder contains patient-data-safe code documenting how RNA-seq expression matrices can be converted into pathway and programme activity scores such as:

- PDAC classical/progenitor, basal/squamous, epithelial, mesenchymal and stromal-rich programmes
- CAF and extracellular matrix (ECM) programmes
- epithelial-mesenchymal transition (EMT), invasion, hypoxia and angiogenesis programmes
- metabolic programmes
- cytolytic, antigen-presentation, checkpoint and interferon-gamma immune programmes

In the PDAC2026 analysis, these scores were used as tumour-level features for interpretation and phenotype assignment. They are not raw expression values, not percentages and not clinical diagnostic classes.

## Inputs

The example script expects:

1. A gene-by-sample expression matrix in TSV format.
2. A GMT gene-set file containing curated programme gene sets.
3. Optionally, a metadata table identifying tumour samples.

The expression matrix should look like this:

| gene | sample_1 | sample_2 | sample_3 |
| --- | --- | --- | --- |
| KRT19 | 10.2 | 8.8 | 14.1 |
| COL1A1 | 6.3 | 11.5 | 4.2 |

The first column must contain gene symbols unless you change `--gene-column`.

## Run

```bash
Rscript run_gsva_ssgsea_programme_scores.R \
  --expression path/to/log_expression_matrix.tsv \
  --gmt ../../templates/pdac_programme_gene_sets_example.gmt \
  --metadata path/to/sample_metadata.tsv \
  --sample-column sample_id \
  --condition-column condition \
  --tumour-label Tumour \
  --out-prefix results/pdac_programme_scores
```

If metadata are supplied, only tumour samples are scored by default. If no metadata are supplied, all expression-matrix columns are scored.

## Outputs

The script writes:

- `<out-prefix>.ssgsea_scores.tsv`: ssGSEA score matrix, programme by sample.
- `<out-prefix>.gsva_scores.tsv`: GSVA score matrix, programme by sample, if GSVA mode succeeds.
- `<out-prefix>.programme_classes.tsv`: per-sample low/intermediate/high calls by cohort tertiles.
- `<out-prefix>.programme_summary.tsv`: median, min, max and high/intermediate/low counts per programme.
- `<out-prefix>.method_notes.txt`: package versions and notes about which mode ran.

## Create A Manuscript-Style Table And Heatmap

After programme scores have been integrated into a tumour-level table, use:

```bash
Rscript create_programme_score_table_figure.R \
  --score-table path/to/integrated_tumour_programme_scores.tsv \
  --sample-column sample_id \
  --phenotype-column phenotype_group \
  --out-dir results/programme_score_summary
```

This writes:

- `programme_score_summary.tsv`
- `programme_scores_long.tsv`
- `programme_score_heatmap.png`
- `programme_score_heatmap.pdf`

This is the safe template corresponding to the manuscript's GSVA/ssGSEA programme-score summary table and heatmap.

## Important Interpretation Notes

- Scores are relative to the expression matrix and gene sets used.
- Low/intermediate/high classes are cohort-relative tertiles, not clinical cutoffs.
- ssGSEA/GSVA values should not be compared directly with ESTIMATE, MCP-counter, CIBERSORT, EPIC, xCell or quanTIseq scores without method-specific interpretation.
- In small cohorts, programme scores are best used for descriptive biology, heatmaps and hypothesis generation.

## Relationship To Phenotype Assignment

The downstream tumour phenotype labels were assigned after these programme scores were integrated with immune/stromal deconvolution outputs:

```text
expression matrix
  -> GSVA/ssGSEA programme scores
  -> immune/stromal deconvolution scores
  -> immune/stromal/EMT phenotype assignment
  -> limma-voom phenotype-group comparison
```

See also:

- `../phenotype_assignment/assign_tme_phenotype_groups.py`
- `../immune_infiltration/paired_tumour_normal_immune_comparison.R`
- `../rnaseq/phenotype_group_comparison_limma_template.R`
