---
layout: default
title: Private inputs
---

# Private Inputs

[Overview](index.md) | [Workflow](workflow.md) | [Figure and table index](figures-and-tables.md)

No private inputs are distributed. The public example configuration contains
empty lists and placeholder locations, so it cannot accidentally reproduce a
patient-level analysis. Do not put a populated configuration in this repository.

## Configuration

[config/private.example.json](https://github.com/tkcaccia/genome-PDAC2026/blob/main/config/private.example.json) documents these keys:

| Key | Meaning |
| --- | --- |
| `paths.project_root` | Private project tree containing the relative input/output paths used by the analysis scripts |
| `paths.figures_out` | New private output directory for vector rendering/conversion |
| `paths.study_tables` | Private presentation input root with `WES/`, `RNAseq/`, `MOFA/`, `Ampulla/` and `Proteomics/` subdirectories |
| `paths.immune_tables` | Private folder containing native scores, paired comparisons, Hallmark scores, composites and group assignments |
| `paths.source_pdfs` | Private source-PDF tree used only by the conversion helper |
| `paths.salmon_quant_dir` | Directory of gene-level Salmon `*.quant.genes.sf` files |
| `paths.hallmark_gmt` | Cached public Enrichr Hallmark 2020 GMT, if already available |
| `cohorts.primary_pdac`, `cohorts.periampullary`, `cohorts.all` | Ordered patient identifiers, without tissue suffixes; preserve the order used in the source matrices |
| `cohorts.epic_converged` | Pairs with converged EPIC fits in both tissues, obtained from the private fit diagnostics |
| `cohorts.kras_hotspot_retained`, `cohorts.kras_hotspot_not_retained` | Variant-evidence-defined comparison sets, not inferred biological wild-type status |
| `values.mofa_exclusion_scenarios` | Scenario-to-sample-list mapping, with tissue suffixes matching factor-score row names; include `full: []` |
| `values.mofa_factor_annotations` | Optional result-specific annotations required by the reconstruction presentation helper; keys such as `factor_1` |
| `values.pca_label_offsets` | Optional patient-label display offsets; no biological grouping is encoded here |

Python uses [pdac_runtime.py](https://github.com/tkcaccia/genome-PDAC2026/blob/main/pdac_runtime.py). R uses
[scripts/shared/runtime.R](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/shared/runtime.R) and `jsonlite`.
Both read the same private JSON. EPIC convergence eligibility and KRAS comparison
sets are independent inclusion criteria.

Set the following variables in the analysis computer's shell:

```bash
export PDAC_CODE_ROOT=/path/to/genome-PDAC2026
export PDAC_PRIVATE_CONFIG=/path/to/private/pdac_config.json
export PYTHONPATH="$PDAC_CODE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
```

The standalone Slurm MSI launcher additionally requires `PDAC_MSI_WORKDIR`.
Its private working directory must contain the paired `sample_manifest.tsv`,
duplicate-marked BAMs/indexes, executable and microsatellite list for the run.

## Input Contracts

| Input | Expected structure |
| --- | --- |
| Salmon quantification | Tab-separated `Name`, `TPM`, `NumReads`; parsers expect numeric patient identifiers followed by `N`/`T` in source sample names |
| Gene TPM/count matrix | Genes in rows, tissue samples in columns; stable gene IDs or symbols with the matching mapping |
| Gene mapping | Ensembl identifier to gene-symbol CSV matching the quantification annotation release |
| Tumour-normal response | Gene by patient matrix of `log2((tumour TPM+1)/(normal TPM+1))`; not a count matrix |
| Variant calls | One VEP-annotated DeepSomatic VCF per pair; consequences, gene symbols and relevant FORMAT depth/allele fields |
| MSI raw summaries | MSIsensor-pro site-count/unstable-count/score files; no hardcoded result substitutes |
| Protein report | Quantified report with the `R.Condition`, `R.Label`, `PG.ProteinAccessions`, `PG.Quantity` and decoy fields used by the parser |
| MOFA model/exports | Private HDF5 model and aligned factor, loading, cluster and sensitivity tables |
| TME native scores | Method-specific score matrices with matching specimen columns |
| Paired TME tests | Feature/cell-type names, means/deltas, paired p-values, method-wise BH FDR and pair counts |
| Programme scores | Tumour-only Hallmark ES and row z-scores, selected GMT coverage, composite score definitions, medians and assignments |
| TCGA inputs | Public Xena expression/RPPA/mutation files, plus clinical metadata with genetic-ancestry labels when needed |

## Analysis and Presentation Paths

Analysis scripts use relative result names, including
dated directories such as `results/rnaseq_revised_11_native_deseq2_20260917/`,
`results/immune_salmon_11_20261012/` and `results/mofa_revised_11_20260916/`.
These names identify the paths expected by the scripts. Private hard-drive
addresses are configured separately.

The figure renderers use a flattened input tree. For example,
`study_tables/RNAseq/native_deseq2_results_11_pdac.csv` is an already derived table,
not a raw Salmon input. Stage the derived inputs while preserving file contents,
provenance and sample order. Rendering requires those intermediate files; it does
not reconstruct them from raw data.

Use a new private output tree for each run. Several scripts write
fixed filenames and can overwrite results if pointed at an existing analysis.
Do not import an analysis module just to discover its command-line arguments;
some modules perform directory creation or data loading at module scope.
