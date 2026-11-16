# WES Analysis

This folder contains the WES methods and publication-panel producers. Raw inputs,
real sample manifests, mutation tables and figures are not public repository files.

| Step | Source | Purpose and private outputs |
| --- | --- | --- |
| Read QC | `raw/20_fastqc_2_5.sbatch` | FastQC for each paired-end library |
| Reference intervals | `raw/01_prepare_capture_2_5.sh`, `remap_bed_to_refseq.py` | Match vendor Exome 2.5 contigs to the GRCh38.p14 dictionary and create intervals |
| Alignment | `raw/02_align_markdup_2_5.sbatch` | BWA-MEM, coordinate sorting, Picard duplicate marking and BAM indexing |
| Coverage | `raw/05_coverage_qc_2_5.sbatch` | Picard target metrics and corrected 1-based KRAS hotspot depths |
| Somatic calling | `raw/03_deepsomatic_2_5.sbatch` | Matched tumour-normal FFPE_WES calling with DeepSomatic 1.10.0 |
| Annotation | `raw/04_vep_2_5.sbatch` | Offline VEP/cache 101, GRCh38, SIFT/PolyPhen and split-VEP fields |
| Variant parsing | `wes_somatic_01_wes_build_mutation_matrix_FIXED.py` | PASS/HIGH-MODERATE calls and available DP/AD/VAF checks |
| Cohort separation | `reanalyze_wes_cohorts.py` | Separate primary PDAC/periampullary tables, counts and mutation matrices |
| Gene sensitivity rules | `complete_wes_outputs.py` | Gene-family exclusion rules, retained matrices and gene summaries; Table S1 |
| Main panels | `generate_primary_wes_panels.py` | Figure 1A-C and variant-level KRAS Table 1 |
| Hotspot evidence | `extract_deepsomatic_kras_fields.py` | FORMAT-level evidence from annotated calls |
| MSI catalogue | `raw/26_msi_scan.sbatch` | Full-reference microsatellite scan with pinned MSIsensor-pro 1.3.0 |
| MSI scoring | `raw/27_msi_paired.sbatch` | Paired `msi -c 20 -b 4` |
| MSI summary | `summarize_msi.py` | Site counts, percentages and study-defined categories for Table S3 |

Figure 1A reads the quality/impact-filtered burden table before gene-family
exclusions; Figures 1B-C use gene-excluded tables. These are different retained
sets and must be identified in the figure legend. Variant counts are not TMB
without an appropriate callable-territory denominator.

No missing AD/VAF fallback is added in this export. The parser is retained as
source, not presented as a newly validated caller. DeepSomatic PASS and secondary
filters do not independently establish biological absence of a KRAS hotspot.

See [input configuration](../../docs/inputs.md), [publication index](../../docs/figures-and-tables.md)
and [Nextflow](../../workflow/nextflow/README.md). The raw jobs expect the private
run layout described there; the Python export reads private configured roots.
