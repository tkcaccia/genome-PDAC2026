---
layout: default
title: Figures and tables
---

# Figures and Tables

[Overview](index.md) | [Workflow](workflow.md) | [Inputs](inputs.md) | [Reproducibility](source-coverage.md)

This index connects each publication figure and table to its input files,
analysis script and rendering function. Panel identifiers below are the reference
for matching the publication to code filenames. Patient-level inputs and figures
are not distributed here.

## Main Figures

| Figure | Input to the plotting step | Rendering code | Computation and reproduction requirements |
| --- | --- | --- | --- |
| 1A mutation burden | `WES/wes_mutation_burden.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `wes_figures()` | [WES cohort summaries](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/reanalyze_wes_cohorts.py); rendering reads the derived burden table |
| 1B oncoprint | `WES/wes_mutation_matrix.blacklist_filtered.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `wes_figures()` | [Gene-filtered mutation matrix](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/complete_wes_outputs.py) |
| 1C gene summary | `WES/gene_summary_table.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `wes_figures()` | [Gene-level aggregation](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/complete_wes_outputs.py) |
| 2A RNA volcano | `RNAseq/native_deseq2_results_11_pdac.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `rnaseq_figures()` | [Native paired DESeq2](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/rnaseq/run_native_deseq2_revised_11.R) |
| 2B RNA response PCA | `RNAseq/tumour_normal_log2_delta_matrix_11_pdac.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `rnaseq_figures()`; [PCA calculation](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/rnaseq/regenerate_figure13_rna_pca_publication.py) | Patient-level tumour-normal responses, not the DESeq2 count input |
| 2C immune composition | `RNAseq/quantiseq_cell_fractions_11_pairs.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `rnaseq_figures()`; [immune renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_immune_vectors.py), `quantiseq_composition()` | [quanTIseq from linear TPM](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/run_salmon_quantiseq_11_pdac.R) |
| 3A-C programmes, TME scores and groups | Hallmark row z-scores, TME feature z-scores, composites, medians and assignments | [Immune renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_immune_vectors.py), `independent_programmes()` | [Programme and grouping calculation](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/generate_figureD_independent_programmes_11_pdac.py); renders the three panels |
| 4 proteomics volcano | Source PDF; matching numeric table is required for statistical reproduction | [PDF conversion only](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/convert_figures.py) | Exact differential-abundance and plotting implementations are not included |
| 5A mutation-DE overlap | Native DESeq2 results and filtered mutation matrix | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `mofa_figures()` | Intersection of mutated genes with the formal differential-expression gene set |
| 5B factor heatmap and 5C clusters | `MOFA/mofa_factors.csv`, `MOFA/mofa_patient_clusters.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `mofa_figures()` | [MOFA fitting and clustering](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/run_mofa_revised_11.py) |
| 5D candidate scores | `MOFA/mofa_top20_composite_score_11_pdac.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `mofa_figures()` | [Cross-modal candidate score](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/regenerate_figure11_composite_score.py) |
| 6A joint RNA PCA | `RNAseq/Figure_13_rna_pca_tcga_11_pdac_coordinates.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `comparison_figures()` | [Joint absolute-expression PCA](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tcga/generate_figure13_rna_pca_tcga_11.py) |
| 6B mutation-frequency contrast | `MOFA/revised_pdac_tcga_mutation_frequency.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `comparison_figures()` | [Mutation comparison](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tcga/compare_mutation_frequencies.py) writes `revised_11_pdac_tcga_mutation_frequency.csv`; stage that derived table under the rendering input filename |
| 6C pathway themes | `MOFA/revised_11_vs_tcga_pathway_theme_concordance.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `comparison_figures()` | [Recurrent-gene enrichment and theme aggregation](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tcga/regenerate_figure17_african_vs_tcga_genuine.py) |
| 7A KRAS status | `Ampulla/kras_variants_all_cohorts.tsv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `comparison_figures()` | [Variant-level KRAS oncoprint](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/periampullary/regenerate_pdac_ampulla_kras_oncoprint.py); no retained call must not be interpreted as proven biological wild type |
| 7B shared pathways | `Ampulla/pdac_ampulla_shared_pathways.csv` | [Study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `comparison_figures()` | [Separate exploratory group enrichment](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/periampullary/generate_pdac_ampulla_exploratory_figures.py) |

The vector renderer colours RNA volcano points using adjusted p-values alone.
The formal DEG table uses both FDR and effect-size thresholds. Point colouring
and membership of the statistical DEG list are therefore separate definitions.

## Supplementary Figures

| Figure | Code and input | Output filename and reproduction requirements |
| --- | --- | --- |
| S1 paired immune/stromal changes | [Paired feature summary](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/generate_figureC_immune_heatmap_11_pdac.py), then [immune renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_immune_vectors.py), `immune_stromal()` | Uses method-wise paired-test tables and the converged-pair EPIC subset |
| S2 normalization | [Normalization helper](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/proteomics/07_qc_normalize_inputs.py), then [study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `normalization_supplement()` | Plots protein integration matrices. Reproduction of precursor/technical-replicate panels additionally requires their numeric source data and replicate definitions |
| S3 CV; S4 peptides; S5 protein groups; S6 hydrophobicity; S7 isoelectric point | [PDF conversion only](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/convert_figures.py) | Extraction-method calculation and plotting implementations are not included |
| S8 peptide/protein comparison | No generating implementation included | Requires numeric source data, replicate definitions and plotting code |
| S9 larger-cohort proteomics volcano | [PDF conversion only](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/convert_figures.py) | Output stem `Supplementary/Figure S8 volcano plot 19`; exact differential-analysis table/code is not included |
| S10 MOFA reconstruction | [Native reconstruction calculations](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/generate_remaining_mofa_publication_figures_11.py), then [study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `mofa_supplements()` | Output stem `Supplementary/Figure S9 mofa native variance`; input `MOFA/mofa_native_variance_explained_11_pdac.csv` |
| S11 exclusion sensitivity | [MOFA exclusion calculation](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/run_mofa_revised_11.py), then [study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `mofa_supplements()` | Output stem `Supplementary/Figure S10 mofa cluster robustness`; requires factors, clusters and sensitivity summary. Exclusions recluster factor scores without refitting MOFA |
| S12 modality support | [Modality-support summary](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/generate_remaining_mofa_publication_figures_11.py), then [study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `mofa_supplements()` | Output stem `Supplementary/Figure S11 mofa modality`; input `MOFA/mofa_modality_support_top200_11_pdac.csv` |
| S13 integrated heatmap | [Integrated feature construction](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/generate_kras_integrated_figures_11.py), then [study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `mofa_supplements()` | Output stem `Supplementary/Figure S12 integrated heatmap`; input `MOFA/kras_integrated_feature_heatmap_11_pdac.csv` |
| S14 TCGA modality counts; S15 ancestry overlap | [Source-metadata summaries](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tcga/generate_revised_tcga_ancestry_figure.py), then [study renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_study_vectors.py), `ancestry_supplements()` | Output stems `Supplementary/Figure S13 tcga ancestry counts` and `Supplementary/Figure S14 ancestry overlap`; requires clinical/ancestry metadata in addition to expression downloads |

## Tables

| Table | Implementation | Scope |
| --- | --- | --- |
| 1 KRAS evidence | [KRAS extraction](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/extract_deepsomatic_kras_fields.py) and [table formatter](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tables/generate_revised_tables_kras_msi.py) | Reads variant-level data; does not embed patient results |
| 2 RNA pathway summary | [Direction-specific enrichment](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/rnaseq/generate_native_rna_pathway_table_11.py), [top-five formatter](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tables/generate_native_rna_pathway_top5_table_11.py), [PNG formatter](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tables/render_native_rna_pathway_top5_png.py) | Formatter ranks terms across its input; an upregulated-only summary requires an upregulated-only input table |
| S1 gene exclusion rules | [Filtering definitions](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/complete_wes_outputs.py), `BLACKLIST_PATTERNS`, `BLACKLIST_GENES`, `blacklist_reason()` | Rules are source code; removed-variant tables remain private |
| S2 clinical characteristics | No table-generation implementation included | Private clinical metadata and table-assembly code are required |
| S3 full-cohort MSI summary | [Paired MSI launcher](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/run_msisensor_pro_20260828.sh), [MSI aggregation](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/complete_wes_outputs.py), [subset formatter](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tables/generate_revised_tables_kras_msi.py) | Aggregation covers configured cohorts. The formatter prints the primary-cohort subset; full-cohort presentation requires an additional table-assembly step |
| S4 quanTIseq paired estimates | [quanTIseq calculations](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/run_salmon_quantiseq_11_pdac.R), then [immune renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_immune_vectors.py), `quantiseq_table()` | Reads paired-test results and formats means, deltas, p-values and BH FDR |
| S5-S6 proteomics search/modification parameters | No table-generation implementation included | Requires the mass-spectrometry search configuration and modification settings |
| S7-S10 proteomics pathway tables | Exact cohort/direction/feature-selection implementations are not included | Requires differential-protein tables, selected gene/protein lists, enrichment backgrounds and library releases |

The [code manifest](source_manifest.json) provides script identities and checksums.
See [reproducibility](source-coverage.md) for execution records and implementation
coverage.
