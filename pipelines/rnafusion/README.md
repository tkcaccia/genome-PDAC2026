# RNA Fusion

Production fusion analysis with `nf-core/rnafusion 4.1.0`.

Files:

- `run_rnafusion.sh`: production launcher wrapper
- `rnafusion.config`: Nextflow process/resource configuration
- `launch_example.sh`: command used to launch/resume the fusion run
- `_common.sh` and `validate_samplesheet.py`: copied shared helpers required by the wrapper

Important production note:

- The final successful fusion tool set was `arriba,fusioncatcher,salmon`.
- `fusionreport` was removed from the production run because its downloader step failed even when the local database was present.
