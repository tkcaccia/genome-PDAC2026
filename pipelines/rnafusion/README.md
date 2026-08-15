# RNA Fusion

Production fusion analysis with `nf-core/rnafusion 4.1.0`.

Files:

- `run_rnafusion.sh`: production launcher wrapper
- `rnafusion.config`: Nextflow process/resource configuration
- `launch_example.sh`: command used to launch/resume the fusion run
- `_common.sh` and `validate_samplesheet.py`: copied shared helpers required by the wrapper

Important production note:

- The final successful fusion tool set was `arriba,fusioncatcher,salmon`.
- Arriba and FusionCatcher were used for candidate fusion detection.
- Salmon-derived outputs supported transcript-level quantification within the workflow.
- `fusionreport` was removed from the production run because its downloader step failed even when the local database was present.
- Candidate fusions should be interpreted as exploratory unless supported by caller evidence, interpretable gene-pair/breakpoint annotation, manual artifact review, relevance to PDAC/actionable genes, and DNA structural-variant concordance where available.
