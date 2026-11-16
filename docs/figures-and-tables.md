---
layout: default
title: Figures and tables
---

# Figures and Tables

[Overview](index.md) | [Workflow](workflow.md) | [Inputs](inputs.md) | [Reproducibility](source-coverage.md)

This index follows the manuscript figure/table numbering. Source-script numbers
inside filenames are execution identifiers, not publication figure numbers.
Only code, method descriptions and empty configuration templates are public.
The input packages contain additional exploratory figures; those are not added
as separate manuscript analyses here.

## Main Figures

| Panel | Generating source | Private input and output |
| --- | --- | --- |
| 1A variant counts | [WES panels](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/generate_primary_wes_panels.py) | `wes/analysis/outputs/primary_pdac/wes_mutation_burden.csv`; before gene-family exclusions |
| 1B oncoprint | [WES panels](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/generate_primary_wes_panels.py) | `wes_mutation_matrix.blacklist_filtered.csv`; binary top-30 recurrence |
| 1C gene counts | [WES panels](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/generate_primary_wes_panels.py) | `gene_summary_table.csv`; post-exclusion top-20 variant counts |
| 2A RNA volcano | [Paired DESeq2](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/rnaseq/run_native_deseq2_revised_11.R) | Gene-level Salmon NumReads; formal count-based tumour-normal contrast |
| 2B RNA response PCA | [Response PCA](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/rnaseq/regenerate_figure13_rna_pca_publication.py) | Paired TPM responses; top-variable features and standardisation |
| 2C immune composition | [quanTIseq](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/run_salmon_quantiseq_11_pdac.R), [immune renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_immune_vectors.py) | Linear symbol TPM; paired cell-fraction matrix |
| 3A-C programmes and TME groups | [Programme calculation](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/generate_figureD_independent_programmes_11_pdac.py), [immune renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_immune_vectors.py) | Tumour-only Hallmark ES, native TME features, z-scores, composites and median rules |
| 4 protein volcano | Exact statistical and plotting producer not supplied | Requires paired protein model, full differential-abundance table and plotted statistic; PDF conversion is not statistical reproduction |
| 5A mutation-DE overlap | [Formal overlap](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/formal_de_mutation_overlap.py) | Formal DESeq2 genes and retained WES genes; excludes response-threshold substitutes |
| 5B factor heatmap | [MOFA fit/export](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/run_mofa_revised_11.py) | Saved factors and clusters |
| 5C clusters | [Cluster panel](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/manuscript_clustering.py) | Aligned factor scores and model cluster assignments |
| 5D candidates | [Composite scoring](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/regenerate_figure11_composite_score.py) | Mapped maximum absolute loadings; within-modality percentile ranks |
| 6A joint RNA PCA | [TCGA PCA](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tcga/generate_figure13_rna_pca_tcga_11.py) | Absolute study/TCGA tumour expression on a common log scale |
| 6B mutation frequencies | [Frequency comparison](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tcga/compare_mutation_frequencies.py) | Study retained calls and public MC3 gene matrix; source paths must match configured run |
| 6C broad pathway themes | [Theme comparison](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tcga/regenerate_figure17_african_vs_tcga_genuine.py) | Recurrent study/TCGA genes; all returned Enrichr terms grouped by keywords |
| 7A KRAS status | [Variant-level strip](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/periampullary/regenerate_pdac_ampulla_kras_oncoprint.py) | `kras_variants_all_cohorts.tsv`; no detected hotspot is not proof of biological wild type |
| 7B pathway sharing | [Group pathway comparison](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/periampullary/generate_pdac_ampulla_exploratory_figures.py) | WES recurrence and top absolute RNA/protein responses; returned-term sharing, not a between-group test |

## Supplementary Figures

| Figure | Source and scope |
| --- | --- |
| S1 paired TME changes | [Paired summary](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/generate_figureC_immune_heatmap_11_pdac.py), [immune renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_immune_vectors.py); method-wise paired tests and converged EPIC subset |
| S2 protein normalisation | [Integration/QC helper](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/proteomics/07_qc_normalize_inputs.py) covers prepared matrices; exact technical-replicate/precursor-panel producers are still required |
| S3-S8 extraction-method optimisation | Exact computation and plotting producers not supplied; numeric source tables and replicate definitions are required |
| S9 larger-cohort protein volcano | Exact paired differential-abundance/plotting producer not supplied |
| S10 reconstruction | [Corrected S10 producer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/generate_figure_s10_reconstruction.py); aligns stored model names and computes per-factor/full-model training R² |
| S11 exclusions | [Cluster stability](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/cluster_stability_panel.py); reclusters fixed factor scores without MOFA refitting |
| S12 modality support | [Ranked input](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/export_modality_support_input.py), then [MOFA summaries](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/generate_remaining_mofa_publication_figures_11.py), `modality_support()` |
| S13 integrated features | [Integrated heatmap](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/generate_kras_integrated_figures_11.py); mutation indicators, RNA responses and tumour protein abundances |
| S14-S15 TCGA modality/ancestry | [Ancestry summary](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tcga/generate_revised_tcga_ancestry_figure.py); public clinical ancestry labels and available modalities |

## Tables

| Table | Source and scope |
| --- | --- |
| 1 KRAS evidence | [WES panel/table producer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/generate_primary_wes_panels.py), [FORMAT extraction](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/extract_deepsomatic_kras_fields.py) |
| 2 RNA pathways | [Enrichment](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/rnaseq/generate_native_rna_pathway_table_11.py), [top-five formatting](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tables/generate_native_rna_pathway_top5_table_11.py), [table rendering](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tables/render_native_rna_pathway_top5_png.py); supply the correct direction-specific table |
| S1 exclusion rules | [Explicit gene rules](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/complete_wes_outputs.py), `blacklist_reason()` |
| S2 clinical characteristics | Exact table-generation source not supplied; real clinical metadata remain private |
| S3 MSI | [Scan](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/raw/26_msi_scan.sbatch), [paired scoring](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/raw/27_msi_paired.sbatch), [full-cohort summary](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/summarize_msi.py) |
| S4 quanTIseq | [Native calculation](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/run_salmon_quantiseq_11_pdac.R), [immune table renderer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/figures/generate_immune_vectors.py), `quantiseq_table()` |
| S5 modifications/search comparison | Exact table-generation/search producer not supplied |
| S6-S9 protein pathways | Exact measurable-background, protein-list selection and enrichment producers not supplied |

## Figure Assets and Execution

Keep figure SVG/PDF files in private `MOFA/figures/` and `WES/figures/` folders,
numeric source data in their private `source_data/` folders, and public generating
code in `scripts/mofa/` and `scripts/wes/`. Do not commit any real plotting input
or figure. Most supplied SVGs contain embedded PNG images; Figure S10 has a
native vector producer. File extension and vector editability are not equivalent.

The [Nextflow graph](https://github.com/tkcaccia/genome-PDAC2026/blob/main/workflow/nextflow/README.md) has site/private-package
dependencies. Its stub/output-audit evidence must not be described as a new full
scientific execution. See [source identities](source_manifest.json) for exact
script hashes and documented configuration adaptations.
