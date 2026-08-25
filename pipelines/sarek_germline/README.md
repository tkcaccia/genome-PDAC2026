# Sarek Germline

Reusable germline whole-exome sequencing (WES) workflow wrapper for
`nf-core/sarek 3.8.1`.

Files:

- `run_sarek.sh`: launcher wrapper
- `sarek.config`: Nextflow process/resource configuration
- `launch_example.sh`: command used to launch/resume the germline run
- `_common.sh` and `validate_samplesheet.py`: copied shared helpers required by the wrapper

Audit boundary: the August 2026 connected storage did not contain the primary
germline variant call format files or execution reports. This folder preserves
the configured method and does not assert that historical calls were
independently re-audited.
