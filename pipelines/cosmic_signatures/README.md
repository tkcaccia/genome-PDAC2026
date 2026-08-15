# COSMIC Signatures

This folder contains the code used to assign known COSMIC mutational signatures to the completed tumor-normal WES callset.

Approach used:

- caller input: `Mutect2` tumor-normal filtered VCFs from the finished Sarek run
- additional filtering: retain `PASS` variants only before signature fitting
- tool: `SigProfilerAssignment 1.1.3`
- reference signatures: `COSMIC v3.5`
- genome build: `GRCh38`
- assay mode: `exome=True` because the somatic data are WES, not WGS

Outputs are generated separately for:

- `SBS96`
- `DBS`
- `ID`

Additional summary outputs can be generated from the finished assignment tables:

- long-form sample/signature activity table
- per-sample top-signature table
- cohort-level signature totals

Files:

- `run_sigprofiler_cosmic_assignment.sh`: end-to-end remote launcher for environment setup, PASS-only VCF preparation, and assignment
- `run_sigprofiler_cosmic_assignment.py`: Python driver that runs SigProfilerAssignment for all requested contexts
- `summarize_sigprofiler_assignment.py`: converts the finished activity tables into cohort-friendly summary TSVs
- `run_sigprofiler_summary.sh`: remote launcher for the summary step
