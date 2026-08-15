# Pipeline and Method Details

## Overview

All production analyses were performed on the remote Ubuntu workstation. This repository contains reusable code and documentation only. Patient data and result tables are excluded.

The project combined standard nf-core pipelines with downstream R/Python/Bash scripts:

- `nf-core/rnaseq 3.24.0` for RNA-seq processing and expression quantification
- `nf-core/sarek 3.8.1` for germline and tumour-normal WES analyses
- `nf-core/rnafusion 4.1.0` for RNA fusion detection
- `SigProfilerAssignment 1.1.3` for mutational signature assignment
- Local downstream scripts for COSMIC annotation, KRAS review, MSI/MMR review, immune/stromal analyses, circos-style plots and integrated summaries

## RNA-seq Quantification

Code location: `pipelines/rnaseq/`

Purpose:

- Process bulk RNA-seq samples
- Generate expression matrices and QC outputs
- Provide gene-level counts/expression for downstream phenotype analyses

Main tools/resources:

- `nf-core/rnaseq 3.24.0`
- GRCh38 reference genome
- GENCODE v46 annotation
- STAR/Salmon outputs depending on pipeline configuration

Important scripts:

- `run_rnaseq.sh`: production wrapper
- `rnaseq.config`: resource and process configuration
- `validate_samplesheet.py`: input samplesheet validation
- `differential_expression_limma_template.R`: paired tumour-normal limma-voom template
- `phenotype_group_comparison_limma_template.R`: phenotype-group and interaction comparison template

## Differential Expression

The statistically preferred downstream workflow is `limma-voom` with TMM normalization from `edgeR`. For matched tumour-normal comparisons, the model is:

```r
~ patient_id + condition
```

where `conditionTumour` estimates tumour versus matched normal while controlling for patient pairing.

For phenotype-group analyses, the recommended approach is:

- Tumour-only group comparison for stromal-high/EMT-high/immune-low versus immune-high/stromal-low tumours
- Paired tumour-normal delta comparison when both tumour and normal exist
- Optional interaction model when sample size is sufficient

The historical Python fallback should be treated as exploratory and not described as DESeq2.

## Germline and Somatic WES

Code locations:

- `pipelines/sarek_germline/`
- `pipelines/sarek_tumor_normal/`

Purpose:

- Run reproducible WES germline and tumour-normal calling
- Annotate variants with VEP
- Generate somatic SNV/indel calls for downstream driver review

Main tools/resources:

- `nf-core/sarek 3.8.1`
- GRCh38
- Ensembl VEP 115
- Singularity containers
- Mutect2 and Strelka for somatic SNV/indel calls

Downstream summaries:

- Mutation burden
- KRAS/TP53/CDKN2A/SMAD4 review
- DDR/MMR gene review
- Exploratory germline/somatic integrated oncoprint

## Copy Number and Structural Variants

Code location: `pipelines/sarek_sv_cna/`

Purpose:

- Recover or run WES-compatible SV/CNA analysis
- Summarize ASCAT copy-number segments
- Summarize Manta SV calls
- Identify interchromosomal translocations
- Generate per-patient circos-style plots

Important scripts:

- `run_sarek_sv_cna_remote.sh`
- `summarize_interchromosomal_translocations.py`
- `make_patient_circos.py`
- `run_enhanced_circos_remote.sh`

Outputs generated during the study included gene-level CNA summaries, SV burden summaries, translocation summaries and circos-style plots. The output files themselves are not included here.

## RNA Fusion

Code location: `pipelines/rnafusion/`

Purpose:

- Detect RNA fusion candidates from RNA-seq data
- Support review of KRAS-wild-type or alternative-driver tumours

Main pipeline:

- `nf-core/rnafusion 4.1.0`

Tools noted in the manuscript:

- Arriba
- FusionCatcher
- Salmon

`fusionreport` was not used in the completed workflow because its database step failed despite local database availability. This should be documented as a limitation of the fusion workflow rather than hidden.

## COSMIC Annotation

Code location: `pipelines/cosmic_annotation/`

Purpose:

- Annotate somatic calls against local COSMIC database tables
- Compare caller-level COSMIC annotation across Mutect2 and Strelka

Important scripts:

- `annotate_mutect2_with_cosmic.py`
- `annotate_strelka_with_cosmic.py`
- `compare_cosmic_callers.py`

Notes:

- COSMIC database files require appropriate licensing/access.
- COSMIC-derived patient result tables should not be committed.

## Mutational Signatures

Code location: `pipelines/cosmic_signatures/`

Purpose:

- Run SigProfilerAssignment on PASS somatic variants
- Fit COSMIC SBS/DBS/ID signatures
- Summarize per-sample contributions

Main scripts:

- `run_sigprofiler_cosmic_assignment.py`
- `summarize_sigprofiler_assignment.py`

Interpretation:

- Signature results are exploratory because this is a small WES cohort.
- Broad hypermutation should not be claimed without orthogonal validation.

## CNA and SV Signatures

Code location: `pipelines/cosmic_sv_cna_signatures/`

Purpose:

- Convert ASCAT segments and Manta calls into CNV48/SV32-style matrices
- Run SigProfilerAssignment-style fitting for copy-number and structural-variant signatures

Important scripts:

- `run_sigprofiler_sv_cna_signatures.py`
- `summarize_sigprofiler_sv_cna.py`

Interpretation:

- These are exploratory WES-derived signatures, not definitive whole-genome scar measurements.

## Immune and Stromal Deconvolution

Methods used in the downstream study:

- ESTIMATE
- MCP-counter
- CIBERSORT LM22 with 100 permutations and quantile normalization disabled
- EPIC
- xCell
- quanTIseq

Input:

- Normalized RNA-seq expression matrices, transformed as required by each method

Interpretation:

- Scores are method-specific and usually not percentages unless the method explicitly reports fractions.
- ESTIMATE score is a relative combined stromal/immune admixture score, not a percent purity.
- Negative scores are possible and indicate lower expression of the corresponding stromal/immune signature relative to the method scale.

## KRAS, MSI/MMR and Hypermutation Review

Purpose:

- Avoid assigning KRAS-wild-type status from a single summary table
- Integrate variant calls, filtered calls and hotspot coverage
- Review MSI/MMR status across all patients

Main principles:

- KRAS-mutated supported: PASS or otherwise well-supported coding hotspot evidence
- KRAS-wild-type supported: adequate hotspot coverage and no coding KRAS hotspot variant
- Filtered/subclonal candidate: candidate exists but requires cautious interpretation
- MSI/MMR candidate: MSIsensor-pro signal plus MMR gene or mutation-pattern support

Final manuscript framing:

- Patient 23 is the strongest computational KRAS-wild-type MSI/MMR-deficient candidate.
- Patient 35 is borderline MSI-elevated but not equivalent to patient 23.
- Broad cohort-level hypermutation is not a robust conclusion after strict filtering.

## Integrated Reporting

The final integrated interpretation combined:

- KRAS status
- TP53/CDKN2A/SMAD4 status
- DDR/MMR status
- TMB
- CNV burden
- SV burden
- RNA subtype/phenotype
- Immune/stromal scores
- CAF/EMT/hypoxia/pathway scores
- Exploratory actionable annotations

These outputs are research-grade and hypothesis-generating. They should not be used as clinical reports.
