# GSVA/ssGSEA Pathway and Programme Scoring

This folder shows how the RNA expression matrix was converted into one activity score per biological programme and tumour. The inputs and outputs are RNA-derived. WES mutation calls do not enter this calculation.

The scored programmes include PDAC classical/progenitor and basal/squamous expression, cancer-associated fibroblast (CAF) and extracellular matrix (ECM), epithelial-mesenchymal transition (EMT), invasion, hypoxia, angiogenesis, metabolism, cytolytic activity, antigen presentation and immune checkpoints.

## Start Here: What You Need

Use a tab-separated gene-by-sample expression matrix:

| gene | sample_1 | sample_2 | sample_3 |
| --- | ---: | ---: | ---: |
| KRT19 | 10.2 | 8.8 | 14.1 |
| COL1A1 | 6.3 | 11.5 | 4.2 |

Requirements:

1. The first column contains HGNC gene symbols or Ensembl gene IDs. For Ensembl IDs, supply the matching GTF with `--gtf`; the script records its checksum and performs the mapping.
2. Every other column is one RNA-seq sample.
3. Values are normalized counts or a continuous log-like expression measure, for example limma-voom log2 counts per million, DESeq2 variance-stabilizing transformation, or log2(TPM + 1). `--transform auto` applies log2(x + 1) when the input scale resembles non-negative normalized counts.
4. Do not gene-wise z-scale the matrix before scoring. The script needs the expression ranking across genes within each sample.
5. Do not supply unnormalized raw integer counts. Run library-size/composition normalization first.

The `DESeq2_normalized_counts.tsv` output produced by `../rnaseq/run_standard_de_from_star_readspergene.R` can be used with the GTF that matches the RNA-seq annotation. The script collapses multiple Ensembl IDs mapping to one symbol by summing normalized expression before transformation. Always inspect `<out-prefix>.gene_set_coverage.tsv`; expression filtering can remove biologically important marker genes even when identifier mapping succeeds.

For the corrected PDAC2026 audit, the same documented gene-level TPM fallback
used for immune deconvolution was supplied here with `--gene-column
gene_symbol --transform log2p1`. This retained 61,510 mapped gene symbols and
scored 14 tumours. It improved classical/progenitor coverage from 5/15 to 15/15
genes and basal/squamous coverage from 6/14 to 14/14 genes. The corrected
scores did not change any of the 14 final 3/3/8 phenotype labels.

Install the required packages in the analysis environment if they are absent:

```r
install.packages("data.table")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GSVA")
```

## What The GMT File Is

A GMT file is a biological reference, not a result calculated from the expression matrix. Each tab-separated line contains:

```text
PROGRAMME_NAME    description_or_source    GENE1    GENE2    GENE3    ...
```

For this repository, the patient-data-safe curated fallback is:

```text
../../templates/pdac_programme_gene_sets_example.gmt
```

That file is version-controlled and can be opened in any text editor. It was assembled from transparent PDAC, stromal, EMT, metabolic and immune marker lists; no R script creates it from patient data. The scoring script `run_gsva_ssgsea_programme_scores.R` **reads** this GMT and matches its HGNC symbols to the expression-matrix row names.

This distinction matters:

- The expression matrix says how strongly each gene is expressed in each sample.
- The GMT says which genes define each biological programme.
- GSVA/ssGSEA combines those two inputs to produce programme-by-sample scores.

The bundled GMT is an explicit curated fallback, not a claim that the listed markers reproduce every full published subtype classifier. For broader public pathway collections, retrieve Hallmark, Reactome or other Molecular Signatures Database (MSigDB) collections with `msigdbr`, record its database version, and convert the selected `gs_name`/`gene_symbol` pairs to GMT. Do not mix an updated MSigDB collection with the curated fallback without recording which sets were used.

## Which R Script To Run

Run `run_gsva_ssgsea_programme_scores.R`. It is the script that reads both the normalized/log expression matrix and the GMT file and calculates the scores:

```bash
Rscript run_gsva_ssgsea_programme_scores.R \
  --expression path/to/gene_tpm.tsv \
  --gmt ../../templates/pdac_programme_gene_sets_example.gmt \
  --gene-column gene_symbol \
  --transform log2p1 \
  --metadata path/to/sample_metadata.tsv \
  --sample-column sample_id \
  --condition-column condition \
  --tumour-label Tumour \
  --out-prefix results/pdac_programme_scores
```

Metadata are optional. When they are supplied, only columns whose metadata condition is `Tumour` are scored. Without metadata, every sample column is scored.

## Outputs

The script writes:

- `<out-prefix>.ssgsea_scores.tsv`: programme-by-sample single-sample gene-set enrichment analysis (ssGSEA) scores.
- `<out-prefix>.gsva_scores.tsv`: programme-by-sample gene set variation analysis (GSVA) scores, if that mode succeeds.
- `<out-prefix>.programme_classes.tsv`: low/intermediate/high calls made separately for each programme using cohort tertiles.
- `<out-prefix>.programme_summary.tsv`: programme-level score summaries and class counts.
- `<out-prefix>.gene_set_coverage.tsv`: the total and expression-matched gene count for every GMT programme.
- `<out-prefix>.expression_used.tsv`: symbol-mapped, collapsed and transformed expression used by GSVA/ssGSEA.
- `<out-prefix>.method_notes.txt`: input paths and MD5 checksums, annotation, transformation, dimensions, R version and GSVA version.

The score values are relative enrichment/activity measures. They are not RNA counts, percentages, immune-cell fractions or clinical diagnostic classes.

If gene-level TPM must first be reconstructed from STAR gene counts, run
`../rnaseq/make_gene_tpm_from_counts.R` with the matching GTF. Complete Salmon
transcript-level TPM remains preferable when it is available for every sample.

## How The Scores Reach Phenotype Assignment

The programme score matrix must be merged with the method-specific immune/stromal matrices. The supplied assembler accepts both feature-by-sample outputs from the reusable runners and archived sample-by-feature tables, transposes where required, and performs the join:

```bash
python ../phenotype_assignment/assemble_tme_score_table.py \
  --programme-table results/pdac_programme_scores.ssgsea_scores.tsv \
  --score-table estimate=results/immune/estimate_scores.tsv \
  --score-table mcp_counter=results/immune/mcp_counter_scores.tsv \
  --score-table epic=results/immune/epic_scores.tsv \
  --score-table xcell=results/immune/xcell_scores.tsv \
  --score-table quantiseq=results/immune/quantiseq_scores.tsv \
  --metadata path/to/sample_metadata.tsv \
  --output results/tumour_score_table.tsv
```

Then `../phenotype_assignment/assign_tme_phenotype_groups.py` standardizes selected score columns across tumours and assigns the exploratory tumour-microenvironment groups.

## Manuscript Table And Heatmap

After the programme scores and phenotype labels have been integrated into a tumour-level table, use:

```bash
Rscript create_programme_score_table_figure.R \
  --score-table path/to/integrated_tumour_programme_scores.tsv \
  --sample-column sample_id \
  --phenotype-column phenotype_group \
  --msi-column validated_MSI_MMR_status \
  --out-dir results/programme_score_summary
```

This writes a programme summary table, a long-format score table, and PNG/PDF heatmaps. The heatmap can include tumour phenotype and microsatellite-instability/mismatch-repair annotation bars.

## Interpretation

- Low/intermediate/high calls are cohort-relative tertiles, not universal cutoffs.
- GSVA/ssGSEA scores should not be numerically averaged with raw ESTIMATE, MCP-counter, CIBERSORT, EPIC, xCell or quanTIseq outputs before method-aware standardization.
- With a small cohort, these scores support descriptive biology and hypothesis generation rather than clinical classification.

See also `../immune_infiltration/README.md`, `../phenotype_assignment/README.md` and `../rnaseq/phenotype_group_comparison_limma_template.R`.

Official software documentation: [Bioconductor GSVA](https://bioconductor.org/packages/GSVA/) and [msigdbr](https://igordot.github.io/msigdbr/).
