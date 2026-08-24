# Immune And Stromal Infiltration

This folder documents both parts of the RNA-derived immune/stromal workflow:

1. Calculate method-specific immune/stromal scores from a normalized expression matrix.
2. Compare the resulting scores between matched tumour and normal samples.

These analyses use RNA expression, not WES variants.

## Preparing Expression Input

All methods expect genes in rows, samples in columns and HGNC gene symbols as row names/first-column values. However, the required numerical scale is not identical for every downstream task.

| Available expression matrix | GSVA/ssGSEA | ESTIMATE/MCP-counter/xCell | EPIC/quanTIseq/CIBERSORT |
| --- | --- | --- | --- |
| Linear TPM-like expression | May be log2-transformed first | Suitable | Preferred input |
| `log2(TPM + 1)` | Suitable as supplied | Can be back-transformed | Back-transform with `2^x - 1` |
| limma log2-CPM | Suitable | Recover TPM-like linear expression for this runner | Recover TPM-like linear expression |
| DESeq2 VST/rlog | Suitable | Recover TPM-like linear expression for this runner | Do not back-transform as if it were TPM |
| Gene-wise z-scores | Not suitable as the starting matrix | Not suitable | Not suitable |

The official `immunedeconv` interface recommends TPM-normalized, non-log data. Therefore, if the only file is a generic normalized/log matrix, first determine how it was transformed. Only select `log2_plus_one` below when the values truly equal `log2(linear abundance + 1)`. If the matrix is VST, rlog or an unknown transformation, use it for GSVA/ssGSEA but recover TPM-like expression from the nf-core/rnaseq Salmon outputs for fraction-oriented deconvolution.

Install the required R packages in the analysis environment if they are absent:

```r
install.packages(c("data.table", "remotes"))
remotes::install_github("omnideconv/immunedeconv")
```

## Calculate Scores From Expression

The script `run_immune_stromal_scores_from_expression.R` provides one documented entry point for ESTIMATE, MCP-counter, EPIC, xCell and quanTIseq through `immunedeconv`:

```bash
Rscript run_immune_stromal_scores_from_expression.R \
  --expression path/to/tpm_expression.tsv \
  --gene-column gene \
  --input-scale linear \
  --methods estimate,mcp_counter,epic,xcell,quantiseq \
  --tumour-mode true \
  --out-dir results/immune_scores
```

For a matrix known to be `log2(TPM + 1)`, use:

```bash
Rscript run_immune_stromal_scores_from_expression.R \
  --expression path/to/log2_tpm_plus_1.tsv \
  --gene-column gene \
  --input-scale log2_plus_one \
  --methods estimate,mcp_counter,epic,xcell,quantiseq \
  --tumour-mode true \
  --out-dir results/immune_scores
```

The runner averages duplicate gene-symbol rows, converts a declared `log2(x + 1)` matrix back to the linear scale, runs each method independently, retains successful outputs if another method fails, and records method status and package version.

`--tumour-mode true` applies the tumour-optimized EPIC/quanTIseq procedure consistently to every supplied sample, including matched normals, so samples remain directly comparable within the run. Record this choice in the methods.

## CIBERSORT LM22

CIBERSORT is an academic-licensed special case. Its R source and LM22 signature cannot be redistributed in this public repository. After obtaining `CIBERSORT.R` and `LM22.txt` through the official CIBERSORT access route, add it explicitly:

```bash
Rscript run_immune_stromal_scores_from_expression.R \
  --expression path/to/tpm_expression.tsv \
  --gene-column gene \
  --input-scale linear \
  --methods cibersort \
  --cibersort-script /secure/path/CIBERSORT.R \
  --cibersort-lm22 /secure/path/LM22.txt \
  --cibersort-permutations 100 \
  --out-dir results/immune_scores
```

Do not commit the licensed CIBERSORT files or patient expression data to GitHub.

## Method-Specific Outputs And Scales

The scoring runner writes feature-by-sample matrices such as:

- `estimate_scores.tsv`
- `mcp_counter_scores.tsv`
- `epic_scores.tsv`
- `xcell_scores.tsv`
- `quantiseq_scores.tsv`
- `cibersort_scores.tsv`, when licensed inputs are configured
- `method_status.tsv` and `method_notes.txt`

ESTIMATE, MCP-counter and xCell primarily support comparisons of the same feature between samples. Relative-mode CIBERSORT supports comparison of cell types within a sample. EPIC and quanTIseq provide fraction-like estimates that support both types of comparison. The methods use different models and scales; their raw values must not be pooled as though all were percentages.

## Paired Tumour-Normal Comparison

Use `paired_tumour_normal_immune_comparison.R` separately on each method's score matrix and a metadata table containing `sample_id`, `patient_id` and `condition`:

```bash
Rscript paired_tumour_normal_immune_comparison.R \
  results/immune_scores/mcp_counter_scores.tsv \
  path/to/sample_metadata.tsv \
  results/mcp_counter_paired \
  feature
```

The script matches tumour and normal samples within patient, calculates tumour-minus-normal deltas, runs paired Wilcoxon signed-rank and paired t-tests feature by feature, and applies Benjamini-Hochberg false-discovery-rate (FDR) correction within each method.

## Connection To Phenotype Assignment

The method-specific matrices are merged with the GSVA/ssGSEA programme matrix by `../phenotype_assignment/assemble_tme_score_table.py`. The merged tumour-level table is then supplied to `../phenotype_assignment/assign_tme_phenotype_groups.py`. Neither the paired comparison script nor the later limma phenotype comparison calculates these upstream scores.

Official software documentation: [immunedeconv input and method guide](https://omnideconv.org/immunedeconv/articles/immunedeconv.html) and [CIBERSORT path configuration in immunedeconv](https://omnideconv.org/immunedeconv/reference/).
