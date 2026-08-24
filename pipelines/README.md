# PDAC2026

Code, configuration, and launch notes for the pancreatic ductal adenocarcinoma short-read analysis completed in 2026.

## Summary

This repository captures the code that was used to prepare, launch, monitor, resume, and recover four production analyses on the remote Ubuntu workstation:

1. `nf-core/rnaseq` for bulk RNA-seq expression quantification
2. `nf-core/sarek` germline calling on WES data
3. `nf-core/sarek` tumor-normal somatic calling on WES data
4. `nf-core/rnafusion` for fusion detection on RNA-seq data
5. `SigProfilerAssignment` for COSMIC mutational signature assignment on the tumor-normal WES callset
6. COSMIC database-driven annotation of the completed somatic callsets
7. `SigProfilerAssignment` for COSMIC `CNV48` and `SV32` assignment on the completed `ASCAT` and `Manta` outputs
8. Per-patient circos-style visualization of SNV burden, CNV segments, SV links, and interchromosomal translocations

The production run finished successfully. Final status:

| Analysis | Pipeline | Status |
| --- | --- | --- |
| RNA expression | `nf-core/rnaseq 3.24.0` | Complete |
| Germline WES | `nf-core/sarek 3.8.1` | Complete |
| Tumor-normal WES | `nf-core/sarek 3.8.1` | Complete |
| RNA fusion | `nf-core/rnafusion 4.1.0` | Complete |
| COSMIC signatures | `SigProfilerAssignment 1.1.3` | Complete |
| COSMIC CNV/SV signatures | `SigProfilerAssignment 1.1.3` | Complete |
| COSMIC annotation (Mutect2) | Local COSMIC v103 tables | Complete |
| COSMIC annotation (Strelka) | Local COSMIC v103 tables | Complete |
| Per-patient circos plots | Local Matplotlib renderer | Complete |

## What This Repo Contains

- `workspace_setup/`: bootstrap and system setup scripts used to prepare the remote workspace
- `cohort_autodraft/`: cohort-specific samplesheet/autodraft generation code
- `shared_runtime/`: shared shell and samplesheet validation helpers used by the pipeline launchers
- `rnaseq/`: RNA-seq wrapper, config, and launch example
- `pathway_scoring/`: normalized/log expression plus curated GMT to GSVA/ssGSEA programme scores
- `immune_infiltration/`: TPM-like expression to method-specific immune/stromal scores, followed by paired tumour-normal comparison
- `phenotype_assignment/`: merge score matrices, calculate meta-scores and assign immune/stromal/EMT tumour phenotypes
- `driver_mutation_review/`: patient-data-safe TP53 mutation-identification example from VEP-annotated somatic calls
- `sarek_germline/`: germline Sarek wrapper, config, and launch example
- `sarek_tumor_normal/`: tumor-normal Sarek wrapper, config, and launch example
- `rnafusion/`: RNA fusion wrapper, config, and launch example
- `cosmic_signatures/`: COSMIC signature assignment code built on the completed somatic WES output
- `cosmic_sv_cna_signatures/`: COSMIC `CNV48`/`SV32` assignment code built on the completed `ASCAT` and `Manta` output
- `sarek_sv_cna/`: SV/CNA Sarek recovery code plus per-patient circos visualization and translocation summary scripts
- `cosmic_annotation/`: COSMIC database-driven post-processing for somatic variant and fusion interpretation
- `monitoring_and_recovery/`: watchdog, recovery, and remote orchestration code

## Reference/Runtime Summary

- Reference build: `GRCh38`
- RNA annotation: `GENCODE v46`
- VEP cache: `Ensembl VEP 115`
- Container runtime used in production: `Singularity`
- Work directory during production: `/media/user/PDAC_SEQ_analysis/work`

## Production Output Locations

These outputs are not stored in GitHub; the repository contains only code and documentation.

- RNA-seq expression: `/home/user/PDAC_SEQ_archive/rnaseq_expression`
- Germline WES: `/media/user/PDAC_SEQ_analysis/results/sarek_germline`
- Tumor-normal WES: `/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal`
- RNA fusion: `/home/user/PDAC_SEQ_native_results/rnafusion_pdac`
- COSMIC signatures: `/media/user/PDAC_SEQ_analysis/results/cosmic_signatures_sigprofiler_assignment_1_1_3`
- COSMIC CNV/SV signatures: `/media/user/PDAC_SEQ_analysis/results/cosmic_sv_cna_signatures_sigprofiler_assignment_1_1_3`
- COSMIC annotation (Mutect2): `/media/user/PDAC_SEQ_analysis/results/cosmic_annotation_mutect2_v103`
- COSMIC annotation (Strelka): `/media/user/PDAC_SEQ_analysis/results/cosmic_annotation_strelka_v103`
- Per-patient circos plots: `/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal_sv_cna/circos_plots_enhanced`

## Important Notes

- Raw patient data are not included here.
- Result files are not included here.
- Passwords, credentialed SSH helpers, and machine-specific secrets were intentionally excluded.
- Several scripts contain the exact production paths used on the remote workstation because the goal of this repo is to preserve the real operational code.
