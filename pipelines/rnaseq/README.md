# RNA-seq

Production RNA expression analysis with `nf-core/rnaseq 3.24.0`.

Files:

- `run_rnaseq.sh`: production launcher wrapper
- `rnaseq.config`: Nextflow process/resource configuration
- `launch_example.sh`: command used to launch/resume the production run
- `_common.sh` and `validate_samplesheet.py`: copied shared helpers required by the wrapper
- `run_standard_de_from_star_readspergene.R`: constructs STAR gene-count matrices and runs paired limma-voom plus DESeq2 sensitivity analysis
- `phenotype_group_comparison_limma_template.R`: compares expression between already assigned tumour phenotype groups; it does not calculate immune/stromal/EMT scores
- `de_diagnostics_sensitivity.R`: checks sample QC, pairing, PCA, limma-vs-DESeq2 effect-size concordance, significant-gene overlap and discordant genes

## What Each RNA-seq Script Answers

| Script | Biological question | Main required inputs | Main outputs |
| --- | --- | --- | --- |
| `run_rnaseq.sh` | Process RNA-seq reads and create expression outputs | nf-core/rnaseq samplesheet and references | STAR/Salmon expression and QC outputs |
| `run_standard_de_from_star_readspergene.R` | Which genes differ between matched tumour and normal samples? | STAR output directory; optional phenotype assignment table | Count matrix, metadata, paired limma-voom table, paired DESeq2 table, normalized/log expression |
| `de_diagnostics_sensitivity.R` | Do limma-voom and DESeq2 agree in effect size and significance? | Count matrix, metadata, limma table, DESeq2 table | QC tables, PCA, limma-vs-DESeq2 comparison, discordant genes |
| `phenotype_group_comparison_limma_template.R` | Among tumours only, which genes differ between pre-defined phenotype groups? | Count matrix plus metadata containing `phenotype_group` | Tumour-normal limma table, tumour phenotype-group limma table or skipped-file explanation |

The paired tumour-normal immune-infiltration comparison is not performed in this
RNA-seq folder. It is implemented in
`../immune_infiltration/paired_tumour_normal_immune_comparison.R`, after
deconvolution scores have already been calculated from the expression matrix.

## Important: Scoring Is Upstream of Group Comparison

The tumour phenotype comparison uses labels such as `StromalHigh_EMTHigh_ImmuneLow` and `ImmuneHigh_StromalLow`. These labels come from upstream immune/stromal/EMT scoring and integrated interpretation. The limma script only compares expression between groups after the labels already exist.

Where the labels come from:

- The starting point is the normalized RNA-seq expression matrix generated from the nf-core/rnaseq/STAR gene-count outputs.
- Immune/stromal context was estimated with tools such as ESTIMATE, MCP-counter, CIBERSORT and immunedeconv methods where available.
- CAF, EMT, hypoxia, angiogenesis and pathway programmes were scored from curated gene sets using the expression matrix. A patient-data-safe GSVA/ssGSEA implementation is provided in `../pathway_scoring/run_gsva_ssgsea_programme_scores.R`.
- Those score tables were integrated at tumour level. Tumours with high stromal/CAF/EMT signal and relatively low immune signal were labelled `StromalHigh_EMTHigh_ImmuneLow`; tumours with stronger immune signal and low stromal signal were labelled `ImmuneHigh_StromalLow`; unclear cases were labelled `Intermediate`.

The labels are therefore derived phenotype summaries used for an exploratory tumour-only comparison. They are not raw percentages, and they are not calculated inside the limma script itself.

In other words:

```text
expression matrix -> immune/stromal/EMT scores -> phenotype_group labels -> limma group comparison
```

If a metadata file has no `phenotype_group` column, the phenotype-group comparison cannot be run.

A safe example phenotype-assignment file is available at:

```text
templates/phenotype_assignment_template.tsv
```
