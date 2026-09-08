# Preanalytical Nucleic Acid Quality Control

This folder documents extraction-level DNA and RNA quality control performed
before library preparation. These measurements describe nucleic-acid yield,
spectrophotometric ratios and fragmentation. They are not sequencing coverage,
library complexity, tumour cellularity or computational tumour purity.

`summarize_provider_nucleic_acid_qc.py` parses the tabular format used in the
CPGR nucleic-acid QC report supplied for the PDAC2026 cohort. It extracts:

- DNA concentration, A260/280, A260/230, DNA integrity number (DIN), quality
  score and provider pass/fail status;
- RNA concentration, A260/280, A260/230, RNA integrity number (RIN), DV200 and
  provider pass/fail status;
- aggregate medians, ranges and pass counts;
- an optional comparison of tumour DNA extraction QC between ASCAT boundary and
  non-boundary fits.

```bash
python3 summarize_provider_nucleic_acid_qc.py \
  --pdf /path/to/private_provider_qc_report.pdf \
  --ascat-audit /path/to/private/ascat_purity_patient_level_audit.tsv \
  --output-dir /path/to/private/preanalytical_qc
```

If the PDF is parsed on a workstation that already has `pdfplumber`, the two
private extracted tables can be summarized later on a restricted server that
does not have PDF libraries:

```bash
python3 summarize_provider_nucleic_acid_qc.py \
  --dna-qc-tsv /path/to/private/provider_dna_qc_private.tsv \
  --rna-qc-tsv /path/to/private/provider_rna_qc_private.tsv \
  --ascat-audit /path/to/private/ascat_purity_patient_level_audit.tsv \
  --output-dir /path/to/private/preanalytical_qc
```

The parser validates that it recovered 28 unique DNA and 28 unique RNA records.
The sample-level tables and any group comparison derived from them must remain
outside the public repository. An aggregate table should be shared only after
manual disclosure review.

Provider pass/fail labels should be interpreted according to the thresholds in
the source report. A sample can fail because one spectrophotometric ratio falls
outside the provider range while still meeting the fragmentation criterion and
remaining technically processable. The original report should therefore remain
the authoritative source for laboratory decisions.
