# genome-PDAC2026

Reusable, patient-data-safe pipeline code and methods documentation for the PDAC2026 African pancreatic cancer multi-omics analysis.

This repository is intended to document how the analyses were performed, not to distribute private data or patient-level results. Raw FASTQ/BAM/CRAM/VCF/MAF/CNV/SV files, expression matrices, reports, figures and real result tables are deliberately excluded.

## Website

The GitHub Pages documentation lives in `docs/` and is designed to be published from the `main` branch using GitHub Pages.

Main pages:

- `docs/index.md`: project overview and pipeline index
- `docs/pipelines.md`: exhaustive pipeline and method map

## Reproducible Quick Start

Create an isolated environment, copy the configuration template, and point all outputs to a restricted directory outside the repository:

```bash
conda env create -f env/environment.yml
conda activate pdac2026-downstream
Rscript env/install_R_packages.R
cp config/config.example.yaml config/config.yaml
./run_all.sh config/config.yaml
```

`run_all.sh` records a success, failure or memory-gate status for every independent step and continues after non-critical failures. The example configuration is intentionally non-runnable until its placeholder input and output paths are replaced.

## Code Layout

- `pipelines/rnaseq/`: nf-core/rnaseq launcher/configuration and differential-expression templates
- `pipelines/pathway_scoring/`: GSVA/ssGSEA programme scoring from normalized RNA-seq expression matrices
- `pipelines/immune_infiltration/`: expression-to-immune/stromal scoring plus paired tumour-normal comparison
- `pipelines/phenotype_assignment/`: documented immune/stromal/EMT phenotype group assignment logic
- `pipelines/driver_mutation_review/`: safe example showing how TP53 mutation evidence was extracted from VEP-annotated Sarek somatic calls
- `pipelines/sarek_tumor_normal/`: nf-core/sarek tumour-normal WES calling
- `pipelines/sarek_germline/`: nf-core/sarek germline WES calling
- `pipelines/sarek_sv_cna/`: SV/CNA recovery, interchromosomal translocation summaries and circos-style plotting
- `pipelines/rnafusion/`: nf-core/rnafusion launcher/configuration
- `pipelines/cosmic_signatures/`: SigProfilerAssignment SBS/DBS/ID signature assignment
- `pipelines/cosmic_sv_cna_signatures/`: SigProfilerAssignment CNV/SV signature assignment
- `pipelines/cosmic_annotation/`: local COSMIC database annotation of somatic calls
- `pipelines/cohort_autodraft/`: samplesheet/autodraft helper code
- `pipelines/shared_runtime/`: shared shell and validation helpers
- `config/`: patient-data-free execution configuration template
- `env/`: conda and package requirements
- `run_all.sh`: configuration-driven downstream RNA orchestrator

Machine administration, credentials, remote watchdogs, mount repair and retired drive paths are intentionally excluded because they are not scientific methods and are unsafe to reuse.

## Provenance Status

The repository documents methods; it is not evidence that every historical result file remains available. The August 2026 audit verified the current RNA reruns (paired differential expression, gene-level TPM, GSVA/ssGSEA, ESTIMATE, MCP-counter, EPIC, xCell, quanTIseq and phenotype assignment). Primary Sarek VCF/MAF files, detailed ASCAT/Manta outputs, nf-core/rnafusion outputs and SigProfiler result directories were not present on the connected storage during that audit. Their code is retained for reproducibility, but biological claims from those unrecovered primary outputs must not be described as independently re-audited.

## RNA-seq Downstream Analysis

RNA-seq downstream analysis has two different layers that should not be mixed up. First, expression matrices are used to calculate immune, stromal, EMT, CAF and pathway-related scores. Those scores are then integrated to assign broad tumour phenotype groups. Only after those groups exist does the limma-voom phenotype comparison test which genes differ between groups.

The reproducible RNA data flow is implemented in:

- `pipelines/rnaseq/differential_expression_limma_template.R`
- `pipelines/rnaseq/phenotype_group_comparison_limma_template.R`
- `pipelines/rnaseq/run_standard_de_from_star_readspergene.R`
- `pipelines/immune_infiltration/run_immune_stromal_scores_from_expression.R`
- `pipelines/immune_infiltration/paired_tumour_normal_immune_comparison.R`
- `pipelines/pathway_scoring/run_gsva_ssgsea_programme_scores.R`
- `pipelines/phenotype_assignment/assemble_tme_score_table.py`
- `pipelines/phenotype_assignment/assign_tme_phenotype_groups.py`
- `docs/pipelines.md`

For reporting, the audited paired tumour-normal analysis used DESeq2 as the prespecified primary count model with edgeR robust quasi-likelihood and limma-voom as sensitivity analyses. Historical custom logCPM/t-test fallback outputs must remain labelled exploratory and must not be described as DESeq2, edgeR or limma.

Limma-voom and DESeq2 can disagree in FDR significance even when their fold changes are highly correlated. Report both prespecified outputs, filtering rules and diagnostic concordance rather than selecting the method with the larger significant-gene count.

Important distinction for students: the phenotype-group limma script does not calculate immune, stromal or EMT scores. Starting with RNA expression, `run_immune_stromal_scores_from_expression.R` calculates method-specific deconvolution scores, while `run_gsva_ssgsea_programme_scores.R` combines normalized/log expression with the version-controlled GMT file in `templates/pdac_programme_gene_sets_example.gmt`. `assemble_tme_score_table.py` merges those outputs by sample, and `assign_tme_phenotype_groups.py` creates the labels. Only after that assignment exists in metadata does the limma script compare expression between groups.

The GMT is not generated from the patient's expression matrix. It is a separate biological reference in which each row names a programme and lists its member genes. The expression matrix supplies measured expression, the GMT supplies biological membership, and GSVA/ssGSEA produces programme activity scores from their intersection.

When the expression matrix contains Ensembl IDs, the scoring script requires the matching GTF and records the expression, GMT and GTF checksums. This prevents an example GMT or later code revision from being mistaken for the reference used to produce an archived result.

Tumour-normal immune/stromal comparisons are a separate analysis. The immune script takes deconvolution score matrices from ESTIMATE, MCP-counter, CIBERSORT LM22, EPIC, xCell or quanTIseq, matches tumour and normal samples by patient, calculates tumour-minus-normal deltas, and performs paired Wilcoxon tests with FDR correction across features within each method.

## Safety

This repository should remain code-only. Before committing, run:

```bash
python3 scripts/audit_public_repository.py
git status --short
find . -type f \( -name "*.tsv" -o -name "*.csv" -o -name "*.vcf*" -o -name "*.maf" -o -name "*.bam" -o -name "*.cram" -o -name "*.fastq*" \)
```

No patient-data files should be staged.
