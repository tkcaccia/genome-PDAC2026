## COSMIC SV/CNA Signatures

This folder contains the SigProfiler-based downstream analysis for the missing structural-variant and copy-number branches of the PDAC tumor-normal WES cohort.

Inputs expected from the dedicated Sarek SV/CNA rerun:

- `Manta` somatic SV VCFs
- `ASCAT` copy-number segments

The workflow does four things:

1. normalizes per-sample `ASCAT` segments into a merged tabular file compatible with `SigProfilerMatrixGenerator`
2. converts `Manta` somatic SV VCFs into plain `.vcf` files for `SVMatrixGenerator`
3. builds `CNV48` and `SV32` matrices
4. runs `SigProfilerAssignment` against the latest COSMIC reference set for both contexts and summarizes the per-sample activities

Operational notes from the production run:

- the published Sarek somatic Manta files are named `*.manta.somatic_sv.vcf.gz`, so the collector supports both the generic `somaticSV` pattern and the real Sarek publish pattern
- the shell launcher is safe to run either from a checked-out script path or when streamed over SSH stdin by setting the script directory defensively
- the code is intended to preserve methods only; result tables and patient-level outputs stay on the analysis workstation

Expected remote result root:

- `/media/user/PDAC_SEQ_analysis/results/cosmic_sv_cna_signatures_sigprofiler_assignment_1_1_3`

Typical remote inputs:

- `/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal_sv_cna`

Key scripts:

- [run_sigprofiler_sv_cna_signatures.sh](/Users/stefano/Documents/SEQ/PDAC2026/cosmic_sv_cna_signatures/run_sigprofiler_sv_cna_signatures.sh)
- [run_sigprofiler_sv_cna_signatures.py](/Users/stefano/Documents/SEQ/PDAC2026/cosmic_sv_cna_signatures/run_sigprofiler_sv_cna_signatures.py)
- [summarize_sigprofiler_sv_cna.py](/Users/stefano/Documents/SEQ/PDAC2026/cosmic_sv_cna_signatures/summarize_sigprofiler_sv_cna.py)
