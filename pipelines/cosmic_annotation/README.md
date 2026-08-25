# COSMIC Annotation

This folder contains post-processing code that can use locally licensed COSMIC
databases to annotate PDAC somatic calls on the analysis workstation.

Current implemented workflow:

- annotate `Mutect2` tumor-normal `PASS` somatic variants using:
  - `Cosmic_MutantCensus`
  - `Cosmic_CancerGeneCensus`
  - `Cosmic_CancerGeneCensusHallmarksOfCancer`
  - `Cosmic_ResistanceMutations`
- annotate `Strelka` tumor-normal `PASS` somatic variants with the same COSMIC resources for caller-to-caller comparison

Outputs produced when the required variant calls and databases are supplied:

- combined variant-level annotation table across the cohort
- recurrent gene summary
- exact COSMIC match summary
- resistance-mutation summary

Planned or easily extensible next steps:

- annotate RNA fusions using `Cosmic_Fusion` and `Cosmic_Breakpoints`
- compare `Mutect2` and `Strelka` COSMIC support rates
- cohort-level summaries restricted to CGC genes, hallmark genes, or resistance-associated variants

Files:

- `annotate_mutect2_with_cosmic.py`: main Python annotation engine
- `run_cosmic_mutect2_annotation.sh`: remote launcher that unpacks the required COSMIC archives and runs the annotation
- `annotate_strelka_with_cosmic.py`: Strelka SNV/indel annotation engine using the same COSMIC tables
- `run_cosmic_strelka_annotation.sh`: remote launcher for the Strelka comparison workflow
- `compare_cosmic_callers.py`: summarizes caller-level and sample-level differences between Mutect2 and Strelka COSMIC annotations
- `run_cosmic_caller_comparison.sh`: remote launcher for the Mutect2-versus-Strelka comparison summary

Audit boundary: the August 2026 connected storage contained the local COSMIC
database collection and recovered summary fields but not the primary annotated
variant tables. The code is preserved, but historical variant-level COSMIC
matches cannot be independently re-audited until source calls are restored or
the annotation is regenerated.
