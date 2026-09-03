# External Transcriptomic Validation and Anatomical Sensitivity

This folder documents the validation added after the internal paired tumour-normal
RNA-seq analysis. It contains code only. No internal expression matrix, clinical
metadata, external-cohort result table or patient identifier is stored here.

## Why this analysis was added

The internal analysis used three paired differential-expression frameworks on the
same filtered genes and the same design (`~ patient_id + condition`). DESeq2 was the
prespecified primary model; edgeR quasi-likelihood and limma-voom were sensitivity
models. A difference in the number of FDR-significant genes is therefore reported as
method sensitivity, while agreement in log2 fold-change direction and magnitude is
reported separately.

Two public tumour-normal pancreatic datasets were then used to ask whether the
internal effects reproduced outside the study:

- GSE15471: 36 usable matched tumour-normal pairs after technical replicates were
  averaged within patient and condition.
- GSE62452: 45 confidently reconstructed matched tumour-normal pairs from the
  unambiguous paired naming scheme.

Each external cohort is analysed independently with paired limma using
`~ patient_id + condition`. Probe identifiers are mapped to HGNC gene symbols;
probes mapping to more than one symbol are removed, and the highest-variance
unambiguous probe is retained per symbol. The external results are then compared
directly with the complete internal DESeq2 table. The exact version-controlled
19-programme GMT used internally is applied to every cohort's ranked statistic so
pathway validation is not confused with gene-level replication.

## Scripts

### Public-cohort preparation

The two preparation scripts start from the public GEO series matrices and create the
external paired-limma tables used by the comparison script:

```bash
Rscript pipelines/external_validation/prepare_gse15471_paired_limma.R \
  /restricted/public_validation/GSE15471

Rscript pipelines/external_validation/prepare_gse62452_paired_limma.R \
  /restricted/public_validation/GSE62452
```

`prepare_gse15471_paired_limma.R` reads explicit patient and sample-type fields from
GEO and averages technical replicate arrays within patient and condition before
modelling. `prepare_gse62452_paired_limma.R` retains only titles following the
unambiguous E-number scheme, normalizes inconsistent leading zeros and requires
exactly one tumour and one normal sample per reconstructed patient. Numeric-title
samples are excluded because their patient pairing cannot be recovered defensibly.

Both scripts remove probes mapping to multiple gene symbols, select the
highest-variance unambiguous probe per symbol, fit the paired limma model and write
the complete unthresholded result table plus metadata and provenance. The GEO
download is public, but it requires network access and the platform-specific
Bioconductor annotation packages listed in `env/requirements_R.txt`.

### `compare_transcriptomic_validation.R`

Inputs:

1. Complete internal DESeq2, edgeR and limma-voom result tables.
2. A two-column external manifest with `cohort` and `result_file`.
3. The exact programme GMT used for the manuscript.
4. An output directory outside this public repository.

The internal tables are expected to use the column names produced by
`pipelines/rnaseq/run_standard_de_from_star_readspergene.R`. External limma tables
must contain a gene-symbol column (`gene` or `gene_symbol`), `logFC`, `t` and
`adj.P.Val`.

Example manifest, created outside the repository:

```text
cohort	result_file
GSE15471	/path/to/GSE15471/results/GSE15471_limma_all_genes.tsv
GSE62452	/path/to/GSE62452/results/GSE62452_limma_all_genes.tsv
```

Example command:

```bash
Rscript pipelines/external_validation/compare_transcriptomic_validation.R \
  /restricted/de/DE_tumour_vs_normal_paired_DESeq2.tsv \
  /restricted/de/DE_tumour_vs_normal_paired_edgeR_QL.tsv \
  /restricted/de/DE_tumour_vs_normal_paired_limma_voom.tsv \
  /restricted/external_manifest.tsv \
  templates/pdac_programme_gene_sets_example.gmt \
  /restricted/results/external_validation
```

The script reports internal method concordance, direct internal-to-external
gene-level concordance, three-cohort overlap when at least two external cohorts are
provided, and exact-GMT preranked GSEA. A positive pathway normalized enrichment
score means enrichment toward genes increased in tumour relative to matched normal.

### `compare_anatomical_sensitivity.R`

This script compares a complete primary DE table with a DE table produced after
restricting metadata to one anatomical group, such as pancreatic-head tumours. Both
tables must have been generated independently with the same filtering rule, design
and software. It also accepts optional long-format paired feature-delta tables to
test whether immune/TME directions persist after the restriction.

```bash
Rscript pipelines/external_validation/compare_anatomical_sensitivity.R \
  /restricted/full/DE_tumour_vs_normal_paired_DESeq2.tsv \
  /restricted/restricted/DE_tumour_vs_normal_paired_DESeq2.tsv \
  /restricted/full/paired_feature_deltas.tsv \
  /restricted/restricted/paired_feature_deltas.tsv \
  /restricted/results/anatomical_sensitivity
```

Use `-` for both feature-delta arguments when only DE results are available. The
feature tables must contain `method`, `feature`, `median_delta`, `p_value` and
`FDR`. Directional retention is descriptive; an unchanged direction does not imply
statistical significance.

## Interpretation boundary

External datasets differ in platform, preprocessing and cohort composition.
Genome-wide correlations can therefore be modest even when robust biological
programmes replicate. Report gene-level overlap, direction and pathway-level
agreement separately. The anatomical analysis is a sensitivity analysis, not a
formal comparison between tumour sites, and small groups should not support ancestry
or histology claims.
