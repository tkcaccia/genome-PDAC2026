# Sarek Tumor-Normal

Reusable matched tumour-normal whole-exome sequencing (WES) workflow wrapper
for `nf-core/sarek 3.8.1`.

Files:

- `run_sarek.sh`: launcher wrapper
- `sarek.config`: Nextflow process/resource configuration
- `launch_example.sh`: command used to launch/resume the tumor-normal run
- `_common.sh` and `validate_samplesheet.py`: copied shared helpers required by the wrapper

Audit boundary: the August 2026 connected storage retained recovered
patient-level summaries but not the primary Sarek VCF/MAF files or execution
reports. This folder preserves the configured method and does not make the
historical result files currently auditable.
