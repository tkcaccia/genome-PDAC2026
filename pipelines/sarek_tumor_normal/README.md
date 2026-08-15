# Sarek Tumor-Normal

Production tumor-normal somatic WES analysis with `nf-core/sarek 3.8.1`.

Files:

- `run_sarek.sh`: production launcher wrapper
- `sarek.config`: Nextflow process/resource configuration
- `launch_example.sh`: command used to launch/resume the tumor-normal run
- `_common.sh` and `validate_samplesheet.py`: copied shared helpers required by the wrapper
