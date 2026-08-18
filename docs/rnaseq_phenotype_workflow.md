# RNA-seq Phenotype Workflow: Scores, Groups and Comparisons

This page explains the part of the analysis that caused confusion: the difference between **immune/stromal/EMT scoring**, **phenotype-group assignment**, and **differential-expression comparison between groups**.

The short version is:

1. RNA-seq expression values were generated first.
2. Immune, stromal, CAF, EMT and related scores were calculated from expression matrices.
3. Tumours were assigned to broad phenotype groups using those scores.
4. A limma-voom script then compared expression between the already assigned groups.

The phenotype-group comparison script does **not** calculate immune, stromal or EMT scores by itself. It expects those groups to already exist in the metadata.

## Why the Confusion Happens

A script named like `run_phenotype_group_comparison.R` sounds as if it should create the phenotype groups. In this analysis, it does not. It only asks:

```text
Given tumour samples that have already been labelled as immune-high/stromal-low or stromal-high/EMT-high/immune-low, which genes differ between those two groups?
```

That means the input metadata must already contain a column such as:

```text
phenotype_group
```

The group labels come from upstream immune/stromal/EMT scoring and integrated interpretation, not from the limma comparison itself.

## Step 1: RNA-seq Quantification

RNA-seq was processed using `nf-core/rnaseq 3.24.0`. STAR produced `ReadsPerGene.out.tab` files. These were merged into a gene-count matrix.

Reusable code:

```text
pipelines/rnaseq/run_standard_de_from_star_readspergene.R
```

Important outputs from this step include:

```text
star_unstranded_gene_counts_matrix.tsv
rnaseq_metadata_for_standard_DE.tsv
DE_tumour_vs_normal_paired_limma_voom.tsv
DE_tumour_vs_normal_paired_DESeq2.tsv
limma_voom_logCPM.tsv
DESeq2_normalized_counts.tsv
```

## Step 2: Immune, Stromal and EMT Scoring

This is a biological scoring step, not a differential-expression test.

The analysis used normalized/log-transformed RNA expression to calculate multiple types of scores:

- Immune and stromal admixture using ESTIMATE.
- Immune and stromal cell-population signals using MCP-counter.
- Immune-cell fraction estimates using CIBERSORT LM22.
- Additional deconvolution using EPIC, xCell and quanTIseq.
- CAF, EMT, hypoxia, angiogenesis, pathway and PDAC-specific programmes from curated gene sets.

These methods use different mathematical scales. For example, ESTIMATE scores are not percentages, MCP-counter scores are abundance-like scores, xCell gives enrichment-like scores, and CIBERSORT/EPIC/quanTIseq produce model-derived fractions under their own assumptions.

The important point for the phenotype comparison is that these scores were used to label tumours into broad groups.

## Step 3: Phenotype-Group Assignment

After scoring, tumours were assigned to broad groups such as:

```text
StromalHigh_EMTHigh_ImmuneLow
ImmuneHigh_StromalLow
Intermediate
```

Older notes or scripts may use equivalent labels with additional underscores:

```text
Stromal_EMT_high_Immune_low
Immune_high_Stromal_low
Intermediate
```

The current reusable scripts normalize these aliases internally, but students should still check the metadata file before running the comparison.

Minimum phenotype metadata:

```text
patient_id    phenotype_group
23            StromalHigh_EMTHigh_ImmuneLow
35            ImmuneHigh_StromalLow
```

A safe example template is provided at:

```text
templates/phenotype_assignment_template.tsv
```

For full differential-expression scripts, metadata normally also includes:

```text
sample_id
patient_id
condition
phenotype_group
```

where `condition` is `Tumour` or `Normal`.

## Step 4: Tumour-Normal Differential Expression

This comparison asks:

```text
Which genes differ between tumour and matched normal tissue?
```

Recommended model:

```r
design <- model.matrix(~ patient_id + condition, data = metadata)
```

Interpretation:

- `patient_id` controls for matching.
- `conditionTumour` estimates tumour versus normal expression.
- This is the main paired tumour-normal comparison.

Reusable scripts:

```text
pipelines/rnaseq/differential_expression_limma_template.R
pipelines/rnaseq/run_standard_de_from_star_readspergene.R
```

## Step 5: Tumour Phenotype-Group Comparison

This comparison asks a different question:

```text
Among tumours only, which genes differ between stromal-high/EMT-high/immune-low tumours and immune-high/stromal-low tumours?
```

Recommended model:

```r
design <- model.matrix(~ phenotype_group, data = tumour_metadata)
```

Reusable script:

```text
pipelines/rnaseq/phenotype_group_comparison_limma_template.R
```

This script:

- reads a count matrix,
- reads metadata with a `phenotype_group` column,
- runs a paired all-patient tumour-normal limma-voom comparison,
- runs a tumour-only comparison between the two extreme phenotype groups if both groups are present,
- writes a skipped file if the phenotype comparison cannot be run.

This script does **not**:

- calculate ESTIMATE scores,
- calculate MCP-counter scores,
- calculate CIBERSORT fractions,
- calculate EMT scores,
- decide which tumours are immune-high or stromal-high.

Those decisions must already be present in the metadata.

## Which Output Should Be Sent to a Student?

For the corrected tumour-normal analysis, send:

```text
DE_tumour_vs_normal_paired_limma_voom.tsv
DE_tumour_vs_normal_paired_DESeq2.tsv
standard_DE_summary.tsv
de_diagnostics_summary.tsv
de_diagnostics_limma_vs_deseq2_logfc.pdf
```

For phenotype-group comparison, send:

```text
DE_tumour_phenotype_stromal_EMT_high_vs_immune_high_limma_voom.tsv
```

or, for the template script:

```text
DE_tumour_stromal_emt_high_vs_immune_high_limma.tsv
```

Also send the metadata file used to create the result, because without `phenotype_group` the comparison cannot be understood.

## How to Explain This in Plain Words

The tumour-normal analysis compares cancer tissue with matched normal tissue. It answers the question: what changes in the tumour?

The phenotype-group analysis compares one type of tumour with another type of tumour. It answers the question: among tumours, do stromal/EMT-high cancers have different expression from immune-high/stromal-low cancers?

Immune/stromal/EMT scoring is the step that defines the tumour categories. The limma script is only the final statistical comparison after those categories have already been defined.

## Common Mistakes

- Do not say the phenotype comparison calculated immune or EMT scores. It did not.
- Do not run the phenotype comparison without checking that `phenotype_group` exists.
- Do not interpret ESTIMATE, MCP-counter, xCell, CIBERSORT, EPIC and quanTIseq scores as the same kind of number.
- Do not treat the phenotype comparison as definitive; the extreme groups were small, so the result is exploratory.
- Do not describe historical custom log2CPM/t-test fallback outputs as DESeq2, edgeR or limma.
