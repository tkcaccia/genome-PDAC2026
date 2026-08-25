# COSMIC Signatures

This folder preserves the workflow configured to assign known COSMIC mutational
signatures to PASS somatic variants from matched tumour-normal whole-exome
sequencing (WES).

Approach used:

- intended caller input: filtered tumour-normal `Mutect2` variant call format
  (VCF) files produced by nf-core/sarek
- additional filtering: retain `PASS` variants only before signature fitting
- tool: `SigProfilerAssignment 1.1.3`
- reference signatures: `COSMIC v3.5`
- genome build: `GRCh38`
- assay mode: `exome=True` because the somatic data are WES, not WGS

Outputs are generated separately for:

- `SBS96`
- `DBS`
- `ID`

Additional summary outputs can be generated when assignment tables exist:

- long-form sample/signature activity table
- per-sample top-signature table
- cohort-level signature totals

Files:

- `run_sigprofiler_cosmic_assignment.sh`: end-to-end remote launcher for environment setup, PASS-only VCF preparation, and assignment
- `run_sigprofiler_cosmic_assignment.py`: Python driver that runs SigProfilerAssignment for all requested contexts
- `summarize_sigprofiler_assignment.py`: converts activity tables into cohort-friendly summary TSVs
- `run_sigprofiler_summary.sh`: remote launcher for the summary step

Set `PDAC2026_RUNTIME_ROOT` and `PDAC2026_RESULTS_ROOT` before running the
shell wrapper. Both must point outside this public repository.

Audit boundary: the August 2026 connected storage did not contain the primary
SigProfilerAssignment output directory. The code and versioned method are
preserved, but historical signature contributions cannot be independently
re-audited until the source VCFs and assignment outputs are restored or
regenerated.
