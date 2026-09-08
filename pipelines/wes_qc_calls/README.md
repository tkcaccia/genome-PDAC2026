# WES depth and genomic-call audit

This workflow summarizes targeted whole-exome sequencing (WES) quality and
tumour-normal genomic calls without modifying primary data. It was written to
prevent two common errors: reporting genome-wide average depth for an exome and
using a whole-genome-style tumour mutational burden (TMB) denominator when a
caller was run with `targeted=false`.

## Inputs

The script expects two result roots:

1. A DRAGEN enrichment-results root containing one sample directory for every
   tumour and normal. Required files include `*.summary.csv`,
   `*.target_bed_coverage_metrics.csv`, and
   `*.target_bed_somatic_callable_regions.bed`. BAM and hard-filtered VCF
   presence is audited but the files are not copied or rewritten.
2. A Hartwig Oncoanalyser root containing SAGE and PURPLE tumour-normal VCFs,
   PURPLE purity/QC files, and driver catalogues for each patient.

The provider sample names are expected to follow `<patient_number><N|T>WESFastq`.
The script converts them to `<patient_number>N_WES` or
`<patient_number>T_WES` internally.

## What is calculated

- Mean target depth, target bases at >=1x/10x/15x/20x/50x/100x, mapping,
  duplicate, Q30, mismatch, enrichment, uniformity, and callability metrics for
  every WES specimen.
- The intersection of tumour and matched-normal somatic-callable BED regions,
  reported in bases, megabases, and as a percentage of the Twist exome target.
- SAGE and PURPLE VCF record counts before filtering, after `PASS`, inside the
  matched callable target, and after both restrictions.
- A strict read-evidence sensitivity count requiring tumour depth >=20, normal
  depth >=10, tumour alternate depth >=3, tumour variant allele frequency
  (VAF) >=0.05, normal alternate depth <=1, and normal VAF <=0.02.
- SNV, MNV, insertion, deletion, nonsynonymous, hotspot, depth, and variant
  allele frequency summaries for final callable calls.
- A detailed PURPLE `PASS` variant table with gene and consequence annotations.
- DRAGEN germline call counts and PURPLE somatic/germline driver-catalogue
  counts. Germline findings remain research candidates until clinical
  confirmation.

The output named `pass_variants_per_pair_callable_mb` is a transparent exome
call-burden measure. It is **not** labelled clinical TMB because its denominator
is the full matched callable Twist target rather than a separately validated
coding territory. Any original Oncoanalyser TMB value is retained only as an
unused provenance field when `targeted=false`.

## Run

```bash
python3 summarize_wes_qc_and_calls.py \
  --enrichment-root /path/to/DRAGEN_WES_enrichment_results \
  --onco-root /path/to/Oncoanalyser_Results \
  --output-dir /path/to/wes_qc_and_calls
```

Default review thresholds are 50x tumour mean target depth, 30x normal mean
target depth, 80% of target bases at >=20x, and 70% paired callable target.
They are descriptive QC flags, not clinical acceptance criteria. Override them
with `--min-tumour-depth`, `--min-normal-depth`, `--min-pct-20x`, and
`--min-pair-callable-pct` when a prespecified study protocol requires different
thresholds. A separate descriptive cohort-outlier flag is calculated on the
log10 strict-burden scale using a median absolute deviation robust z-score; its
default threshold is 2.5 and can be changed with `--outlier-robust-z`.

## Outputs

- `wes_sample_qc.tsv`: one row per WES specimen.
- `wes_matched_pair_qc.tsv`: tumour-normal depth, callable territory, PURPLE
  purity/QC, and descriptive review flags.
- `somatic_call_summary.tsv`: caller-level call counts and quality evidence.
- `somatic_filter_counts.tsv`: VCF filter frequencies.
- `somatic_pass_pair_callable_variants.tsv`: final PURPLE PASS calls inside the
  matched callable target, including gene, consequence, depth, and VAF.
- `germline_and_driver_call_summary.tsv`: germline and driver-catalogue counts.
  Catalogue rows are summarized separately from rows with driver likelihood
  >=0.5; they are not automatically described as clinically reportable.
- `patient_level_wes_genomics_summary.tsv`: one concise row per patient joining
  depth, paired callability, purity/QC, strict somatic burden, driver-catalogue
  counts, missing-call status, and interpretation limitations.
- `input_completeness.tsv`: expected-input presence by patient.
- `wes_qc_and_call_audit_summary.md`: concise audit narrative.
- `output_manifest.tsv`: generated output inventory.

## Interpretation safeguards

- WES depth is reported only for the capture target, not across the genome.
- Somatic counts are restricted to positions callable in both tumour and normal.
- Caller `PASS` status is preserved; filtered records are not promoted to calls.
- PURPLE purity values and QC warnings remain visible alongside call counts.
- PURPLE driver-catalogue counts are labelled unvalidated in the patient-level
  summary because low-purity or non-somatic PURPLE states can inflate them.
- WES-derived copy-number and structural-variant results should not be described
  as having whole-genome resolution.
- Output tables contain patient-level results and must not be committed to the
  public repository.
