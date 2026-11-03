---
layout: default
title: Analysis workflow
---

# Analysis Workflow

[Overview](index.md) | [Figure and table index](figures-and-tables.md) | [Inputs](inputs.md)

This page describes the current collaborators' source, not the superseded
STAR/Sarek/custom-expression workflow. Reorganization did not execute these steps.

## WES and MSI

Previously aligned, duplicate-marked and DeepSomatic-called WES files are external
inputs. The full alignment/calling/VEP job definitions were not present in this
source selection, so they are not replaced with earlier Sarek commands.

[The VCF parser](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/wes_somatic_01_wes_build_mutation_matrix_FIXED.py)
reads per-patient VEP-annotated DeepSomatic calls and applies the supplied PASS,
impact and available depth/allele-support filters. [Cohort summaries](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/reanalyze_wes_cohorts.py)
split PDAC from periampullary cases. [The completion helper](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/complete_wes_outputs.py)
applies the explicit gene-family/large-gene sensitivity rules and produces binary
gene matrices and gene summaries. This exclusion list does not prove that every
removed variant is an artefact.

[The MSI launcher](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/wes/run_msisensor_pro_20260828.sh) uses the paired
`msisensor-pro msi` command, tumour/normal `.dedup.bam` files, the project
`tools/microsatellites.list`, minimum coverage 20 and four threads.
The aggregation code uses MSS <5%, intermediate 5 to <20%, and high >=20% as
exploratory categories. Exact executable version and construction/checksum of
the site list are not supplied by the launcher. Private fallback scores were
removed: missing original MSI files are not silently replaced by embedded values.

## RNA Expression and DESeq2

[Salmon matrix assembly](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/rnaseq/02_rnaseq_build_tpm_matrix.py) reads
`Name`, `TPM` and `NumReads` from gene-level `*.quant.genes.sf` files. TPM and
estimated counts are different outputs and must not be interchanged.

[Native R DESeq2](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/rnaseq/run_native_deseq2_revised_11.R) imports
`NumReads`, strips Ensembl version suffixes, merges repeated gene IDs, replaces
missing entries with zero, filters counts >=10 in at least the number of pairs,
and rounds estimated counts. It fits `~ patient + condition`, with normal as the
reference, uses a tumour-versus-normal Wald contrast, and applies BH correction.
The exported DEG list uses FDR <0.05 and absolute log2 fold change >=1.
No custom Python differential-expression fallback is included. The retained
implementation directly imports gene-level estimated counts; it should not be
described as a tximport/tximeta analysis because the code does not do that.

[RNA responses](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/rnaseq/03_rnaseq_build_delta_matrix.py) are separately
calculated as `log2((tumour TPM + 1)/(normal TPM + 1))`. The matrix helper operates
on configured pairs and records filtered and unfiltered responses. The
[PCA producer](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/rnaseq/regenerate_figure13_rna_pca_publication.py) consumes
the archived primary-cohort response matrix, selects the 2,000 most variable
genes and standardizes features across patients. One point is one patient's
response, not one tissue library. These response values are not DESeq2 inputs.

[RNA enrichment](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/rnaseq/generate_native_rna_pathway_table_11.py) submits
upregulated and downregulated gene symbols separately to the named GO, Reactome
and KEGG Enrichr libraries. The supplied implementation uses the library
background; reorganization has not rerun it with a measured-gene background.
Formatting a top-five table does not define the full tested enrichment family.

## From Expression to Immune and Stromal Scores

Start with linear gene-level TPM plus an Ensembl-to-symbol mapping, not with WES
variants. [The quanTIseq script](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/run_salmon_quantiseq_11_pdac.R)
sums mapped gene symbols, records excluded IDs and signature coverage, writes
`salmon_gene_symbol_tpm_11_pairs.csv`, and runs the native R `quantiseqr` TIL10
method. Both tissues use the same tumour-mode signature settings and mRNA scaling.

[The multi-method R script](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/run_salmon_tme_methods_11_pdac.R) reads
that symbol TPM matrix and its matched-sample metadata. EPIC receives linear TPM;
ESTIMATE, MCP-counter and xCell receive `log2(TPM + 1)`. It records coverage and
EPIC fit diagnostics, and writes scores in each method's native scale. EPIC and
quanTIseq provide estimated fractions; ESTIMATE, MCP-counter and xCell scores are
not percentages. CIBERSORT is not part of this current supplied rerun.

Paired tumour-normal comparisons use two-sided Wilcoxon signed-rank tests with
normal approximation and BH correction separately within each method. EPIC also
has a sensitivity table restricted to pairs with converged fits in both tissues.
The [paired-plot summary](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/generate_figureC_immune_heatmap_11_pdac.py)
uses the configured converged EPIC subset and expresses paired changes in
within-feature standard-deviation units for display; it does not combine native
scores into cell percentages.

## Hallmark GMT and Tumour Phenotype Groups

The exact implementation is
[generate_figureD_independent_programmes_11_pdac.py](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tme/generate_figureD_independent_programmes_11_pdac.py).
The GMT is a biological reference, not an expression-derived file. Its
`load_hallmark_sets()` function reads the cached Enrichr
`MSigDB_Hallmark_2020` GMT or retrieves the collection through GSEApy, selects the
named sets, and writes the selected GMT and mapping coverage into the private
analysis tree. An unrelated example GMT is not included as a substitute.

The 19 selected programmes are EMT, hypoxia, TGF-beta signalling, angiogenesis,
KRAS signalling up and down, E2F targets, G2-M checkpoint, MYC targets V1, p53,
glycolysis, oxidative phosphorylation, fatty-acid metabolism, cholesterol
homeostasis, unfolded-protein response, interferon-alpha response,
interferon-gamma response, inflammatory response and complement. The exact
case-sensitive names are in `PROGRAMMES` in the script.

Tumour-only symbol expression is transformed as `log2(TPM + 1)`. GSEApy ssGSEA
uses rank sample/correlation normalization, weight 0.25 and 15-500 mapped genes.
Unnormalized enrichment scores are retained and standardized across tumours.

Each of the 12 selected ESTIMATE/MCP-counter features is standardized as
`(value - mean)/sample SD`. Nine equally weighted features define the immune
composite: ESTIMATE ImmuneScore and MCP-counter T cells, CD8 T cells, cytotoxic
lymphocytes, B lineage, NK cells, monocytic lineage, myeloid dendritic cells and
neutrophils. Three equally weighted features define the stromal composite:
ESTIMATE StromalScore and MCP-counter endothelial cells and fibroblasts.
Standardized Hallmark EMT is the third axis. Missing required values stop
assignment; they are not imputed. Zero-variance standardized rows become zero.

Scores at or above their cohort median are high; scores below it are low.
Stromal/EMT-high, immune-low requires high stromal and EMT and low immune scores.
Immune-high, stromal-low requires high immune and low stromal scores irrespective
of EMT. Remaining tumours are mixed/intermediate. The script exports the feature
definitions, raw scores, z-scores, composites, medians and individual assignments.
These are descriptive tumour-only categories, not validated diagnostic subtypes
or tumour-normal response classes.

## Proteomics and MOFA

[The protein-matrix helper](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/proteomics/04_proteomics_build_protein_matrix.py)
parses protein groups and sample labels from the supplied quantified report.
It is not the missing original proteomics differential-abundance/search pipeline.
The [normalization helper](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/proteomics/07_qc_normalize_inputs.py) documents
the supplied integration/QC transformation, not a replacement for extraction-method
or precursor-level plots.

[The MOFA implementation](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/run_mofa_revised_11.py) aligns three
patient-level views: binary retained WES genes; tumour-normal RNA responses; and
tumour-only protein `log2(abundance + 1)` values. Missing protein abundances become
zero before transformation. The continuous views are standardized, Bernoulli
and Gaussian likelihoods are used as appropriate, and the supplied model has
five factors. The script writes input matrices, QC filters, model, factors,
loadings, clustering and interpretation outputs.

[Candidate scoring](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/regenerate_figure11_composite_score.py) takes
maximum absolute loadings, maps features to symbols, percentile-normalizes within
modality and combines evidence across modalities. [Reconstruction summaries](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/generate_remaining_mofa_publication_figures_11.py)
and [independent-seed fitting](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/mofa/generate_mofa_seed_stability_11.py)
are separate helpers. Sample-exclusion sensitivity in the main model script
reclusters retained rows of existing factor scores; it does not refit MOFA.
Result-specific factor names and exclusion membership stay in private configuration.

## External and Periampullary Comparisons

[TCGA downloading](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tcga/18_pull_tcga_paad_xena.py) and
[matrix preparation](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/tcga/19_prepare_tcga_paad_matrices.py) supply public
RNA, RPPA and mutation inputs. The ancestry/modality script additionally requires
public clinical metadata with genetic-ancestry labels; that file is not supplied
by the Xena downloader. Joint expression PCA, mutation-frequency differences and
pathway-theme counts are descriptive cross-platform comparisons, not matched
ancestry-association tests.

[PDAC/periampullary comparisons](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/periampullary/compare_pdac_ampulla_14.py)
and [pathway plots](https://github.com/tkcaccia/genome-PDAC2026/blob/main/scripts/periampullary/generate_pdac_ampulla_exploratory_figures.py)
keep the two cohorts separate. Their lists are private configuration, never
hardcoded patient identifiers in public code.
