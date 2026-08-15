# RNA-seq

Production RNA expression analysis with `nf-core/rnaseq 3.24.0`.

Files:

- `run_rnaseq.sh`: production launcher wrapper
- `rnaseq.config`: Nextflow process/resource configuration
- `launch_example.sh`: command used to launch/resume the production run
- `_common.sh` and `validate_samplesheet.py`: copied shared helpers required by the wrapper
- `run_standard_de_from_star_readspergene.R`: constructs STAR gene-count matrices and runs paired limma-voom plus DESeq2 sensitivity analysis
- `de_diagnostics_sensitivity.R`: checks sample QC, pairing, PCA, limma-vs-DESeq2 effect-size concordance, significant-gene overlap and discordant genes
