# Tumour Phenotype Assignment: Immune/Stromal/EMT Groups

This folder documents, with patient-data-safe example code, how tumours were classified into broad tumour-microenvironment expression phenotypes for exploratory downstream comparisons.

The final manuscript-facing cohort labels were:

- `StromalHigh_EMTHigh_ImmuneLow`: 3 tumours
- `ImmuneHigh_StromalLow`: 3 tumours
- `Intermediate`: 8 tumours

These are cohort-relative expression phenotype groups. They are not clinical diagnostic classes, not percentages and not formal PDAC transcriptional subtype labels.

## Input Data Used In The Real Workflow

The classification used tumour RNA-seq downstream outputs, not WES variant calls. The score layer came from normalized/log expression matrices and immune/stromal/pathway scoring outputs, including:

- ESTIMATE immune, stromal and combined ESTIMATE scores.
- MCP-counter immune and stromal population scores, including fibroblast signal.
- CIBERSORT LM22 immune-cell estimates.
- EPIC, xCell and quanTIseq deconvolution outputs through immunedeconv.
- Curated gene-set scores for CAF, extracellular matrix (ECM), epithelial-mesenchymal transition (EMT), hypoxia, angiogenesis and pathway programmes.

Because these tools use different scales, the workflow never treated all outputs as percentages. Each feature was first standardized across tumours before integration.

## Classification Logic

The patient-data-safe implementation in `assign_tme_phenotype_groups.py` uses this transparent rule:

1. Read a tumour-level score table with one row per tumour.
2. Convert selected immune, stromal and EMT-related score columns to z-scores across tumours.
3. Calculate three meta-scores:
   - `immune_meta_score`: average z-score of immune-related features.
   - `stromal_meta_score`: average z-score of stromal/fibroblast/CAF/ECM features.
   - `emt_meta_score`: average z-score of EMT/invasion features.
4. Calculate two contrast scores:
   - `immune_high_stromal_low_score = immune_meta_score - stromal_meta_score`
   - `stromal_emt_high_immune_low_score = mean(stromal_meta_score, emt_meta_score) - immune_meta_score`
5. Assign phenotype groups using cohort-relative quantile thresholds:
   - `ImmuneHigh_StromalLow`: high immune-high/stromal-low contrast.
   - `StromalHigh_EMTHigh_ImmuneLow`: high stromal/EMT-high/immune-low contrast.
   - `Intermediate`: all remaining tumours.

For the PDAC2026 cohort, the final reviewed assignment was 3 immune-high/stromal-low, 3 stromal-high/EMT-high/immune-low and 8 intermediate/mixed tumours. The final labels were used as metadata for downstream limma-voom phenotype-group comparisons.

## Why This Was Done

PDAC bulk RNA-seq is strongly shaped by tumour microenvironment composition. A direct tumour-only comparison can miss the fact that some tumours are dominated by stromal/CAF/EMT signals, while others show stronger immune infiltration and lower stromal signal. These groups were therefore used to ask an exploratory question:

> Which tumour expression programmes differ between stromal-high/EMT-high/immune-low tumours and immune-high/stromal-low tumours?

This is coherent with the manuscript because the study aims to integrate genotype and phenotype, not only catalogue mutations. The phenotype labels provide a simple, transparent way to relate RNA expression, immune/stromal deconvolution, EMT biology and genomic features in a small cohort.

## Run The Example

```bash
python assign_tme_phenotype_groups.py \
  --scores path/to/tumour_score_table.tsv \
  --sample-column sample_id \
  --immune-columns estimate_immune_score,mcp_t_cells,cibersort_cd8_t_cells,quantiseq_cd8_t_cells \
  --stromal-columns estimate_stromal_score,mcp_fibroblasts,epic_caf,xcell_fibroblasts,caf_score,ecm_score \
  --emt-columns emt_score,invasion_score,hypoxia_score \
  --out-prefix results/tme_phenotype_assignment
```

The script writes:

- `tme_phenotype_assignment.assignments.tsv`: tumour-level meta-scores and assigned phenotype labels.
- `tme_phenotype_assignment.group_counts.tsv`: number of tumours assigned to each group.

The column names above are examples. Use the column names from your own deconvolution and gene-set score tables.

## Important Interpretation Notes

- The labels are relative to the analysed cohort and depend on the available score columns.
- The labels should be used for exploratory biology and visualization, not for clinical classification.
- The strongest comparison in the manuscript was between the two extremes: `StromalHigh_EMTHigh_ImmuneLow` versus `ImmuneHigh_StromalLow`.
- Tumours without a clear extreme signal were intentionally kept as `Intermediate` rather than forced into a binary class.

