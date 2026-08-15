# Differential Expression Clarification for Students

## Short Answer

Your concern is correct: a historical downstream file name suggested `DESeq2_or_limma`, but one exploratory fallback analysis used a custom Python calculation of log2CPM and paired statistics. That fallback is not the same as DESeq2, edgeR or limma and should not be described as such.

The corrected, reproducible analysis supplied in this repository uses `edgeR` normalization plus `limma-voom`, which is a standard RNA-seq differential-expression workflow.

## Why DESeq2, edgeR or limma-voom Are Preferred

RNA-seq counts are discrete and mean-variance dependent. Standard RNA-seq differential-expression tools model this structure better than a simple t-test on transformed counts:

- DESeq2 models counts using negative binomial generalized linear models and estimates dispersion.
- edgeR also models count data using negative binomial methods and TMM normalization.
- limma-voom estimates the mean-variance trend after count normalization and then applies precision weights in a linear modelling framework.

For this project, `limma-voom` is especially practical because it handles paired designs and flexible phenotype-group models cleanly.

## Main Paired Tumour-Normal Model

For matched tumour-normal samples, the recommended model is:

```r
design <- model.matrix(~ patient_id + condition, data = metadata)
```

This means:

- `patient_id` controls for baseline expression differences between patients.
- `condition` estimates the tumour-versus-normal effect.
- The coefficient `conditionTumour` is interpreted as tumour expression relative to normal expression after accounting for pairing.

The provided script is:

```text
pipelines/rnaseq/differential_expression_limma_template.R
```

Example:

```bash
Rscript pipelines/rnaseq/differential_expression_limma_template.R \
  counts.tsv \
  metadata.tsv \
  results/differential_expression
```

Required metadata columns:

- `sample_id`
- `patient_id`
- `condition`, with values `Tumour` and `Normal`

## Phenotype Group Comparisons

The student question also mentions immune/stromal/EMT tumour phenotype group comparisons. These are not the same as the tumour-normal comparison.

The phenotype analysis asks:

```text
Among tumour samples only, do stromal-high/EMT-high/immune-low tumours differ from immune-high/stromal-low tumours?
```

Recommended tumour-only model:

```r
design <- model.matrix(~ phenotype_group, data = tumour_metadata)
```

The provided reusable script is:

```text
pipelines/rnaseq/phenotype_group_comparison_limma_template.R
```

Required metadata columns:

- `sample_id`
- `patient_id`
- `condition`
- `phenotype_group`

Recommended group labels:

- `Stromal_EMT_high_Immune_low`
- `Immune_high_Stromal_low`
- `Intermediate`

## Optional Two-Way/Interaction-Style Question

A two-factor question can be framed as:

```text
Does the tumour-normal expression change differ between phenotype groups?
```

This is an interaction question. A clean model is possible only if there are enough matched tumour and normal samples in each group:

```r
design <- model.matrix(~ patient_id + condition * phenotype_group, data = metadata)
```

But with only 14 patients and only three tumours in each extreme phenotype group, this is statistically fragile. These results should be called exploratory.

## How to Explain This in Methods

Recommended wording:

```text
Tumour-normal differential expression was assessed using a paired design where possible. Reproducible downstream scripts use TMM normalization and limma-voom linear modelling with patient identity included as a blocking factor. Phenotype-group comparisons were performed as exploratory tumour-only and paired-delta analyses because the extreme immune/stromal groups were small.
```

Avoid saying:

```text
DESeq2 was used
```

unless the DESeq2 script was actually run and its outputs are the ones being reported.

Avoid saying:

```text
All genes were significantly different
```

Instead say:

```text
A large number of genes showed tumour-normal expression differences, but statistical significance should be interpreted using FDR correction, effect-size thresholds and the model used to generate the result.
```

## Practical Interpretation

The tumour-normal analysis identifies genes whose expression differs between pancreatic tumour and matched normal tissue. This includes tumour-cell-intrinsic biology, but also differences in tissue composition, tumour purity, fibrosis, immune infiltration and stromal admixture.

The phenotype-group analysis identifies genes associated with tumour microenvironment states. For example, stromal-high/EMT-high/immune-low tumours may show extracellular matrix, fibroblast and invasion-related genes, whereas immune-high/stromal-low tumours may show stronger immune-cell or inflammatory expression programmes.

Because the cohort is small, these analyses should be treated as hypothesis-generating rather than definitive biomarkers.
