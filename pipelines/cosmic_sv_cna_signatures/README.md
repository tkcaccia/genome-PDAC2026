## COSMIC SV/CNA Signatures

This folder contains the SigProfiler-based downstream analysis for the missing structural-variant and copy-number branches of the PDAC tumor-normal WES cohort.

Important interpretation note:

These are WES-derived CNV/SV signature workflows. They are useful for exploratory prioritization, but they are not equivalent to whole-genome sequencing based CNV/SV signature or genomic-scar analyses. Exome capture reduces genome-wide breakpoint sensitivity and copy-number resolution outside targeted regions.

Inputs expected from the dedicated Sarek SV/CNA rerun:

- `Manta` somatic SV VCFs
- `ASCAT` copy-number segments

The workflow does four things:

1. normalizes per-sample `ASCAT` segments into a merged tabular file compatible with `SigProfilerMatrixGenerator`
2. converts `Manta` somatic SV VCFs into plain `.vcf` files for `SVMatrixGenerator`
3. builds `CNV48` and `SV32` matrices
4. runs `SigProfilerAssignment` against the configured COSMIC v3.5 reference set for both contexts and summarizes the per-sample activities

The resulting CNV48 and SV32 activities should be described as exploratory WES-compatible signatures, not definitive WGS-grade structural-scar measurements.

Operational notes retained from the historical run:

- the published Sarek somatic Manta files are named `*.manta.somatic_sv.vcf.gz`, so the collector supports both the generic `somaticSV` pattern and the real Sarek publish pattern
- the shell launcher is safe to run either from a checked-out script path or when streamed over SSH stdin by setting the script directory defensively
- the code is intended to preserve methods only; result tables and patient-level outputs stay on the analysis workstation

Runtime paths are supplied through `PDAC2026_RUNTIME_ROOT`,
`PDAC2026_RESULTS_ROOT` and optional `SAREK_SV_CNA_ROOT` environment
variables. Patient-derived inputs and outputs must remain outside this
repository.

Key scripts:

- [`run_sigprofiler_sv_cna_signatures.sh`](run_sigprofiler_sv_cna_signatures.sh)
- [`run_sigprofiler_sv_cna_signatures.py`](run_sigprofiler_sv_cna_signatures.py)
- [`summarize_sigprofiler_sv_cna.py`](summarize_sigprofiler_sv_cna.py)

Audit boundary: the August 2026 connected storage did not contain the primary
ASCAT/Manta matrices or SigProfiler output directory. The method is preserved,
but historical CNV48/SV32 contributions cannot be independently re-audited
until those files are recovered or regenerated.
