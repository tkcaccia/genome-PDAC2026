# Pipeline and Method Details

## Overview

All production analyses were performed on the remote Ubuntu workstation. This repository contains reusable code and documentation only. Patient data and result tables are excluded.

The project combined standard nf-core pipelines with downstream R/Python/Bash scripts:

- `nf-core/rnaseq 3.24.0` for RNA-seq processing and expression quantification
- `nf-core/sarek 3.8.1` for germline and tumour-normal WES analyses
- `nf-core/rnafusion 4.1.0` for RNA fusion detection
- `SigProfilerAssignment 1.1.3` for mutational signature assignment
- Local downstream scripts for COSMIC annotation, KRAS review, MSI/MMR review, immune/stromal analyses, circos-style plots and integrated summaries

## Data Flow

| Source pipeline/tool | Main file classes used | Downstream transformation | Final analysis use |
| --- | --- | --- | --- |
| `nf-core/rnaseq 3.24.0` | STAR `ReadsPerGene.out.tab`, Salmon quantification outputs, RNA-seq QC summaries | Merge STAR gene-count files into a count matrix; generate normalized/log-expression matrices | Differential expression, RNA-seq QC, subtype scoring, pathway scoring, immune/stromal deconvolution and RNA-fusion context |
| `nf-core/sarek 3.8.1` germline | Normal-sample germline VCF/annotation outputs and QC summaries | Filter/harmonize germline-relevant calls for driver, DDR and MMR panels | Germline-context review, DDR/MMR interpretation and integrated oncoprint annotations |
| `nf-core/sarek 3.8.1` tumour-normal | Mutect2 and Strelka somatic VCFs, VEP annotations, alignment/QC outputs | Filter calls, convert to MAF-like tables and summarize by patient/sample | Somatic SNV/indel summary, driver review, KRAS validation, TMB estimation and COSMIC annotation |
| `nf-core/sarek 3.8.1` CNA/SV | ASCAT segments, purity/ploidy outputs, Manta candidate SV calls | Summarize ASCAT segments into gene-level CNA and burden features; convert Manta calls into affected-gene/SV-type tables | CNA burden, key-gene CNA status, purity/ploidy, SV burden, interchromosomal translocation summaries and exploratory CNA/SV signatures |
| `nf-core/rnafusion 4.1.0` | Arriba and FusionCatcher candidate fusion calls; Salmon transcript-level outputs | Filter by gene-pair/breakpoint interpretability, caller support, artifact review and DNA-SV concordance where available | Exploratory RNA fusion review and potentially clinically relevant fusion-candidate prioritization |
| `SigProfilerAssignment 1.1.3` | SBS96/DBS/ID catalogues from somatic variants; CNV48 and SV32-style matrices from ASCAT/Manta summaries | Fit COSMIC v3.5 signatures; use exome mode for SBS/DBS/ID analyses | Exploratory SBS, DBS, indel, CNA and SV signature contribution summaries |
| `MSIsensor-pro` | Tumour-normal WES alignments and exome-filtered microsatellite site list | Calculate microsatellite instability scores per tumour and integrate with MMR-gene evidence | MSI/MMR classification, hypermutation review and patient-prioritization tables |
| Custom R/Python/Bash scripts | Pipeline outputs, metadata, sample-pairing tables and curated PDAC gene sets | Harmonize identifiers, check pairing, run downstream statistics, score gene sets, integrate deconvolution outputs and generate figures | Integrated genotype-phenotype summaries, manuscript-style figures, exploratory annotations and safe reusable pipeline export |

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
- `run_standard_de_from_star_readspergene.R`: end-to-end STAR gene-count matrix construction with paired limma-voom and paired DESeq2 sensitivity analysis
- `../pathway_scoring/run_gsva_ssgsea_programme_scores.R`: GSVA/ssGSEA pathway and programme scoring template from normalized expression matrices

## Differential Expression

The statistically preferred downstream workflow is `limma-voom` with TMM normalization from `edgeR`, with `DESeq2` used as a sensitivity analysis when raw integer counts are available. For matched tumour-normal comparisons, the model is:

```r
~ patient_id + condition
```

where `conditionTumour` estimates tumour versus matched normal while controlling for patient pairing.

For phenotype-group analyses, the recommended approach is:

- Tumour-only group comparison for stromal-high/EMT-high/immune-low versus immune-high/stromal-low tumours
- Paired tumour-normal delta comparison when both tumour and normal exist
- Optional interaction model when sample size is sufficient

The historical Python fallback should be treated as exploratory and not described as DESeq2.

## Immune/Stromal/EMT Scoring and Phenotype Assignment

The immune/stromal/EMT phenotype groups were not created by the limma differential-expression scripts. They came from upstream RNA-seq score integration. They are derived tumour-phenotype labels used to organise the cohort for exploratory expression comparisons; they are not separate raw measurements, not percentages, and not formal PDAC transcriptional subtype names.

The workflow is:

```text
nf-core/rnaseq expression outputs
  -> normalized/log expression matrices
  -> immune, stromal, CAF, EMT, hypoxia and pathway scores
  -> integrated tumour phenotype labels
  -> limma-voom phenotype-group comparison
```

The scoring layer used normalized RNA-seq expression matrices from nf-core/rnaseq/STAR-derived gene counts. In practical terms, the inputs were gene-by-sample expression tables generated after RNA-seq processing, not VCF/CNV/SV files. The scoring step then combined several method families:

- ESTIMATE for relative immune, stromal and combined ESTIMATE scores.
- MCP-counter for immune and stromal population abundance-like scores, including fibroblast-related signal.
- CIBERSORT LM22 for model-derived immune-cell fractions.
- EPIC, xCell and quanTIseq through immunedeconv for additional immune/stromal deconvolution.
- Curated gene-set scoring for CAF, EMT, hypoxia, angiogenesis, pathway and PDAC-related programmes using GSVA/ssGSEA-style scoring. The safe reusable script is `pipelines/pathway_scoring/run_gsva_ssgsea_programme_scores.R`, and an example non-patient GMT file is provided at `templates/pdac_programme_gene_sets_example.gmt`.

These outputs were interpreted on their own method-specific scales. They were not treated as the same kind of number and were not all interpreted as percentages.

For paired tumour-normal immune-infiltration analysis, each method-specific score table was analysed separately with `pipelines/immune_infiltration/paired_tumour_normal_immune_comparison.R`. The script matches each tumour sample to its normal sample by `patient_id`, calculates tumour-minus-normal deltas for every immune, stromal or tumour-microenvironment feature, then performs paired Wilcoxon signed-rank tests and paired t-tests. Benjamini-Hochberg FDR correction is applied across features within each method. This is the analysis referred to in the manuscript when describing tumour-normal immune/stromal contrasts.

The score outputs were then reviewed together at patient/tumour level. The assignment was based on concordant patterns across the immune-deconvolution outputs and curated gene-set scores, not on a single hard cutoff from one package. In practical terms:

- Tumours with consistently high fibroblast, stromal, CAF, ECM and EMT-like signals, together with relatively lower immune-cell signal, were labelled stromal-high/EMT-high/immune-low.
- Tumours with stronger immune-infiltration signal and comparatively lower stromal/fibroblast signal were labelled immune-high/stromal-low.
- Tumours that did not clearly fit either extreme were labelled intermediate.

The broad group labels used by downstream scripts are:

- `StromalHigh_EMTHigh_ImmuneLow`
- `ImmuneHigh_StromalLow`
- `Intermediate`

For the PDAC2026 cohort, the final reviewed tumour phenotype split was:

- 3 tumours classified as `StromalHigh_EMTHigh_ImmuneLow`.
- 3 tumours classified as `ImmuneHigh_StromalLow`.
- 8 tumours classified as `Intermediate` or mixed.

The reproducible, patient-data-safe example is provided in `pipelines/phenotype_assignment/assign_tme_phenotype_groups.py`. It takes tumour-level immune, stromal and EMT score columns, z-scales each feature across tumours, averages related features into meta-scores, calculates immune-high/stromal-low and stromal/EMT-high/immune-low contrast scores, and assigns the broad phenotype groups either by cohort-relative quantile thresholds or by ranked top-N extremes.

The phenotype-group limma scripts use this already assigned `phenotype_group` metadata column. They do not calculate ESTIMATE, MCP-counter, CIBERSORT, EPIC, xCell, quanTIseq, CAF or EMT scores themselves. A safe example metadata template is provided in `templates/phenotype_assignment_template.tsv`.

The canonical spellings are `ImmuneHigh_StromalLow` and `StromalHigh_EMTHigh_ImmuneLow`. Older notes or manually edited files may contain spelling variants such as `ImmuneHigh_StromaLow` or `StromaHigh_ENTHigh_ImmuneLow`; the R templates now normalise those variants to the canonical labels before running the comparison.

## Germline and Somatic WES

Code locations:

- `pipelines/sarek_germline/`
- `pipelines/sarek_tumor_normal/`
- `pipelines/driver_mutation_review/`

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

### TP53 Mutation Identification Example

Code location: `pipelines/driver_mutation_review/`

The TP53 mutation review used VEP-annotated somatic calls from the tumour-normal Sarek workflow. The patient-data-safe example script is:

- `identify_tp53_mutations_example.py`

The script demonstrates the logic used to extract TP53 evidence:

- read a tab-separated VEP-annotated somatic variant table;
- select rows where the annotated gene symbol is `TP53`;
- classify protein-impacting consequences such as missense, frameshift, stop-gained and splice-site variants;
- prioritize `PASS` somatic calls as strongest evidence;
- keep lower-confidence rows separately for audit/review rather than silently discarding them;
- write both detailed TP53 rows and a conservative patient-level TP53 mutation-status table.

This example also explains the manuscript distinction between a strict `somatic TP53 mutation` and a broader integrated `TP53 alteration`, where the latter may include copy-number or germline-context flags retained in the patient-level molecular summary.

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

Interpretation limits:

- These are WES-compatible ASCAT and Manta inferences, not WGS-grade copy-number/SV catalogues.
- ASCAT on exome data is useful for allele-specific copy-number, purity/ploidy and gene-level CNA summaries at captured loci.
- Manta on exome data can nominate candidate SVs supported by reads in sufficiently covered regions, but breakpoint sensitivity is lower and less uniform than WGS.
- Broad genome-wide structural-pattern, chromothripsis-like, translocation-burden and genomic-scar conclusions should be treated as exploratory unless confirmed by WGS or orthogonal assays.

## RNA Fusion

Code location: `pipelines/rnafusion/`

Purpose:

- Detect RNA fusion candidates from RNA-seq data
- Support review of KRAS-wild-type or alternative-driver tumours

Main pipeline:

- `nf-core/rnafusion 4.1.0`

Final successful tools:

- Arriba
- FusionCatcher
- Salmon

`fusionreport` was not used in the completed workflow because its database step failed despite local database availability. This should be documented as a limitation of the fusion workflow rather than hidden.

Candidate filtering/prioritization:

- Review gene-pair and breakpoint interpretability.
- Review caller support and remove likely recurrent artifacts or low-confidence calls.
- Prioritize PDAC-relevant or potentially clinically relevant genes including NTRK1, NTRK2, NTRK3, ALK, RET, ROS1, FGFR2, NRG1, ERBB2, BRAF and RAF1.
- Cross-check RNA fusion candidates against DNA structural-variant calls where available.
- Treat unvalidated calls as exploratory candidates, not confirmed clinical fusions.

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
- CNV48/SV32 matrices derived from ASCAT/Manta WES outputs have lower and less uniform genome-wide sensitivity than WGS-derived matrices.
- Use these signatures for prioritization and hypothesis generation, not definitive structural-scar assignment.

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

- Outputs must be interpreted method-by-method rather than pooled as a single immune percentage.
- ESTIMATE scores are relative immune, stromal and combined admixture scores, not percent purity.
- MCP-counter outputs are method-specific abundance scores.
- xCell outputs are enrichment-style scores and should not be treated as direct cell fractions.
- CIBERSORT, EPIC and quanTIseq can report model-derived fractions, but those fractions depend on each method's reference signature, normalization and mixture assumptions.
- Negative scores are possible and indicate lower expression of the corresponding stromal/immune signature relative to the method scale.
- Nominal paired tumour-normal trends should be separated from FDR-significant findings after multiple-testing correction.

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

Final analysis framing:

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
- Exploratory potentially clinically relevant annotations

These outputs are research-grade and hypothesis-generating. They should not be used as clinical reports or treatment recommendations; any candidate finding may warrant orthogonal validation and appropriate clinical review before further use.
