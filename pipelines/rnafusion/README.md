# RNA Fusion

Reproducible wrapper and configuration for RNA fusion analysis with
`nf-core/rnafusion 4.1.0`.

Files:

- `run_rnafusion.sh`: launcher wrapper
- `rnafusion.config`: Nextflow process/resource configuration
- `launch_example.sh`: command used to launch/resume the fusion run
- `_common.sh` and `validate_samplesheet.py`: copied shared helpers required by the wrapper

Archived configuration note:

- The archived launch configuration selected `arriba,fusioncatcher,salmon`.
- Arriba and FusionCatcher are configured for candidate fusion detection.
- Salmon supports transcript-level quantification within the workflow.
- `fusionreport` was omitted from the archived configuration after a historical
  downloader failure.
- Candidate fusions should be interpreted as exploratory unless supported by caller evidence, interpretable gene-pair/breakpoint annotation, manual artifact review, relevance to PDAC or potentially clinically relevant genes, and DNA structural-variant concordance where available. These candidates are not treatment recommendations and may warrant orthogonal validation before further use.

Audit boundary: the August 2026 connected storage did not contain the
nf-core/rnafusion result directory, execution report, trace or module logs.
Consequently, the audit could verify the wrapper and configured modules but not
which modules completed for each sample or any historical fusion candidate.
