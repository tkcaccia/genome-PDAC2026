# Tumour Phenotype Assignment: Immune/Stromal/EMT Groups

This folder documents how tumour RNA-seq downstream scores were integrated into broad tumour-microenvironment (TME) expression phenotypes. The classification used RNA expression-derived features, not WES variants.

The manuscript-facing cohort labels were:

- `StromalHigh_EMTHigh_ImmuneLow`: 3 tumours
- `ImmuneHigh_StromalLow`: 3 tumours
- `Intermediate_or_mixed`: 8 tumours

These labels are cohort-relative exploratory summaries. They are not percentages, clinical diagnostic classes or formal PDAC transcriptional subtypes.

## Complete Data Flow From An Expression Matrix

```text
normalized RNA expression
  |-- linear TPM-like matrix -> ESTIMATE/MCP-counter/EPIC/xCell/quanTIseq
  |                              (optional CIBERSORT only with licensed inputs)
  |-- normalized log matrix + curated GMT -> GSVA/ssGSEA programme scores
  -> merge every score table by RNA sample ID
  -> z-standardize selected features across tumours
  -> calculate immune, stromal and EMT meta-scores
  -> assign 3/3/8 phenotype groups
  -> add phenotype_group to tumour metadata
  -> run the tumour-only limma-voom phenotype comparison
```

The WES-derived KRAS, TP53, CDKN2A, SMAD4, copy-number and DNA-damage-repair variables were compared with these phenotypes later. They did not determine the phenotype labels.

## Which Script Does What

| Script/file | Starts from | Produces |
| --- | --- | --- |
| `../immune_infiltration/run_immune_stromal_scores_from_expression.R` | Linear TPM-like expression, or known `log2(TPM + 1)` | ESTIMATE, MCP-counter, EPIC, xCell, quanTIseq and optional licensed CIBERSORT score matrices |
| `../pathway_scoring/run_gsva_ssgsea_programme_scores.R` | Normalized/log expression plus a GMT file | Programme-by-sample GSVA/ssGSEA score matrices |
| `../../templates/pdac_programme_gene_sets_example.gmt` | Version-controlled biological marker definitions | The gene membership read by the GSVA/ssGSEA script; this file is not generated from expression |
| `assemble_tme_score_table.py` | Method-specific score matrices | One tumour-level table containing all prefixed score columns |
| `assign_tme_phenotype_groups.py` | The merged tumour-level score table | Immune, stromal and EMT meta-scores plus phenotype labels |
| `../rnaseq/phenotype_group_comparison_limma_template.R` | Raw counts plus metadata containing the assigned labels | Gene-level expression comparison between the two extreme phenotype groups |

## Step 1: Calculate Immune/Stromal Scores

Follow `../immune_infiltration/README.md`. The basic command for linear TPM-like expression is:

```bash
Rscript ../immune_infiltration/run_immune_stromal_scores_from_expression.R \
  --expression path/to/tpm_expression.tsv \
  --gene-column gene \
  --input-scale linear \
  --methods estimate,mcp_counter,epic,xcell,quantiseq \
  --out-dir results/immune_scores
```

If the available matrix is exactly `log2(TPM + 1)`, use `--input-scale log2_plus_one`. Do not use that conversion for DESeq2 VST/rlog or an unknown log transformation.

The score layer supported by the reusable code includes:

- ESTIMATE ImmuneScore, StromalScore and combined ESTIMATEScore.
- MCP-counter immune/stromal population scores, including fibroblasts.
- optional CIBERSORT LM22 immune-cell estimates when the licensed files are available.
- EPIC, xCell and quanTIseq outputs through `immunedeconv`.

These methods use different scales and assumptions. Their raw output values were not treated as a common percentage scale.

In the audited PDAC2026 rerun, ESTIMATE, MCP-counter, EPIC, xCell and quanTIseq
all completed. CIBERSORT did not run and did not contribute to classification.
quanTIseq contributed to the separate tumour-normal immune comparison, but the
final phenotype meta-score used the prespecified ESTIMATE, MCP-counter, EPIC,
xCell and ssGSEA columns listed in the feature manifest.

## Step 2: Calculate CAF/ECM/EMT And Other Programme Scores

Follow `../pathway_scoring/README.md` and run:

```bash
Rscript ../pathway_scoring/run_gsva_ssgsea_programme_scores.R \
  --expression path/to/log_expression_matrix.tsv \
  --gmt ../../templates/pdac_programme_gene_sets_example.gmt \
  --metadata path/to/sample_metadata.tsv \
  --sample-column sample_id \
  --condition-column condition \
  --tumour-label Tumour \
  --out-prefix results/pdac_programme_scores
```

The GMT is the text file containing the gene members of each programme. It is not another R script and it is not inferred from this cohort. `run_gsva_ssgsea_programme_scores.R` is the R script that reads it and calculates the sample scores.

Starting from only a normalized/log matrix requires checking its scale. GSVA
and ssGSEA can use log2-CPM, DESeq2 VST, or `log2(TPM + 1)` directly. Immune
deconvolution instead requires non-negative linear TPM-like abundance. Only a
matrix known to equal `log2(TPM + 1)` may be inverted with `2^x - 1`. A VST,
rlog, generic log-CPM, or unknown transform must not be inverted as TPM; recover
Salmon TPM or create the documented gene-level TPM fallback from raw counts and
the matching GTF with `../rnaseq/make_gene_tpm_from_counts.R`.

## Step 3: Merge Score Matrices By Sample

`assemble_tme_score_table.py` accepts either feature-by-sample matrices from the reusable runners or archived sample-by-feature tables. It detects the latter when a `sample_id` column is present, gives every feature a method prefix, merges samples by exact ID and optionally retains tumour samples from metadata:

```bash
python assemble_tme_score_table.py \
  --score-table estimate=results/immune_scores/estimate_scores.tsv \
  --score-table mcp_counter=results/immune_scores/mcp_counter_scores.tsv \
  --score-table epic=results/immune_scores/epic_scores.tsv \
  --score-table xcell=results/immune_scores/xcell_scores.tsv \
  --score-table quantiseq=results/immune_scores/quantiseq_scores.tsv \
  --programme-table results/pdac_programme_scores.ssgsea_scores.tsv \
  --metadata path/to/sample_metadata.tsv \
  --sample-column sample_id \
  --condition-column condition \
  --tumour-label Tumour \
  --output results/tumour_score_table.tsv
```

If licensed CIBERSORT output is unavailable, omit that one argument and continue with the independent methods. Non-numeric columns such as `patient_id` and `condition` in archived sample-by-feature tables are not treated as scores. The output has one row per tumour and columns such as `estimate_immune_score`, `mcp_counter_fibroblasts`, `epic_cancer_associated_fibroblast` and `programme_emt_invasion`. Inspect the actual header before selecting columns because package versions can change feature labels.

## Step 4: Calculate Meta-Scores And Assign Groups

`assign_tme_phenotype_groups.py` performs the classification:

1. It z-standardizes every selected feature across tumours. This puts method-specific values on comparable relative scales without calling them percentages.
2. It averages selected immune features into `immune_meta_score`.
3. It averages stromal, fibroblast, CAF and ECM features into `stromal_meta_score`.
4. It averages EMT/invasion features into `emt_meta_score`.
5. It calculates `immune_meta_score - stromal_meta_score` for the immune-high/stromal-low contrast.
6. It calculates `mean(stromal_meta_score, emt_meta_score) - immune_meta_score` for the stromal/EMT-high/immune-low contrast.
7. It assigns the highest-ranking non-overlapping tumours to the two extremes and leaves the others intermediate.

An illustrative command is:

```bash
python assign_tme_phenotype_groups.py \
  --scores results/tumour_score_table.tsv \
  --sample-column sample_id \
  --immune-columns estimate_immune_score,mcp_counter_t_cells,xcell_immune_score,programme_cytolytic_activity,programme_ifng_response \
  --stromal-columns estimate_stromal_score,mcp_counter_fibroblasts,epic_cancer_associated_fibroblast,xcell_stroma_score,programme_caf_ecm \
  --emt-columns programme_emt_invasion,programme_pdac_mesenchymal \
  --target-per-extreme 3 \
  --out-prefix results/tme_phenotype_assignment
```

The column list above reflects the biological feature families used in the audited workflow but remains illustrative because package feature labels can change. Inspect the exact machine-safe names in `tumour_score_table.tsv`. The project-facing 3/3/8 split corresponds to `--target-per-extreme 3`; the script prevents the same tumour from entering both extreme groups.

The outputs are:

- `tme_phenotype_assignment.assignments.tsv`: sample-level meta-scores, contrast scores and phenotype labels.
- `tme_phenotype_assignment.group_counts.tsv`: group sizes.

## Robustness Audit

The 3/3/8 rule deliberately selects three ranked tumours at each extreme; it
does not discover three natural clusters. Run the supplied audit to confirm
exact reconstruction and quantify sensitivity to the 75th-percentile rule,
leaving out one feature, and leaving out one method family:

```bash
python audit_tme_phenotype_robustness.py \
  --scores results/tumour_score_table.tsv \
  --assignments results/tme_phenotype_assignment.assignments.tsv \
  --feature-manifest results/tme_phenotype_assignment.feature_manifest.tsv \
  --target-per-extreme 3 \
  --upper-quantile 0.75 \
  --out-prefix results/tme_phenotype_robustness
```

The sample-level robustness table is a private result and must not be committed
to this public repository.

In the corrected PDAC2026 audit, the recorded ranked labels reproduced exactly
with no missing selected features. All 12 leave-one-feature-out tests retained
the labels. Excluding ESTIMATE, EPIC or xCell also retained every label;
excluding MCP-counter changed two labels. The 75th-percentile alternative
produced a 4/4/6 split and agreed with 12/14 ranked labels. These results support
the descriptive extremes while confirming that the exact 3/3/8 split is a
prespecified ranking rule rather than a naturally inferred cluster solution.

## Step 5: Compare Expression Between Phenotype Groups

Add `phenotype_group` from the assignment output to the tumour metadata, then run `../rnaseq/phenotype_group_comparison_limma_template.R` using raw counts. This last script answers which genes differ between `StromalHigh_EMTHigh_ImmuneLow` and `ImmuneHigh_StromalLow` tumours. It does not calculate immune, stromal or EMT scores.

## Why This Analysis Was Included

Bulk PDAC RNA-seq contains malignant epithelial, fibroblast, endothelial and immune-cell signals. The integrated labels summarize whether a tumour lies near a stromal/CAF/EMT-rich extreme, an immune-rich/stromal-low extreme, or neither. This supports the manuscript's genotype-phenotype integration, but the small 3-versus-3 expression comparison remains exploratory and requires cautious multiple-testing interpretation.
