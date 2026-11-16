# Nextflow Workflow

The DSL2 graph exposes the supplied site-aware orchestration for WES QC,
alignment/duplicate marking, target and hotspot coverage, paired DeepSomatic,
VEP, full-reference MSI, WES summaries, MOFA and figure export. Prepared RNA and
protein inputs are reused; it does not process raw RNA or mass-spectrometry data.

## Required Private Assets

`stage.py` intentionally expects the complete private handoff package. Set
`--package_root` to its extracted folder, `--asset_root` to the complete original
large-asset tree and `--manifest` to the private sample manifest. None is supplied
in this repository. The full package includes its own execution-layout scripts
and figure-completion dependencies. The adapter runs those packaged scripts;
it does not reconstruct that private package from the grouped public export.

The grouped public Python files use private JSON configuration. Their source
hashes identify the packaged origins; they are a code-only inspection/export
surface, not a second numerically validated run. This distinction prevents a
false claim that the repository alone contains patient inputs or all site assets.

The site environment requires module trees, Picard, the VEP 101 cache, Singularity,
the DeepSomatic image, the pinned MSI binary/source tree and a matching Python
environment. Adapt site installation paths before using another cluster.
The recorded WES/MOFA Python versions are in
[the environment snapshot](../../env/requirements_recorded_wes_mofa.txt); it is
not an independently resolved lockfile or the RNA/TME environment.
The scan job records provenance using `git rev-parse` and therefore expects
the MSI source checkout, not just a copied binary directory.

## Commands

Read-only completed-output audit, with reports and the Nextflow work directory
placed in private storage:

```bash
cd /path/to/private/nextflow_execution
nextflow run /path/to/genome-PDAC2026/workflow/nextflow/main.nf \
  --mode audit --asset_root /path/to/private/completed_run \
  --reports_dir /path/to/private/audit_reports
```

A fresh scientific run must be requested deliberately:

```bash
nextflow run /path/to/genome-PDAC2026/workflow/nextflow/main.nf -profile slurm \
  --mode full --python /path/to/analysis/python \
  --package_root /path/to/private/extracted_handoff \
  --asset_root /path/to/private/large_assets \
  --manifest /path/to/private/sample_manifest.tsv \
  --run_root /path/to/NEW_workspace/results/bed_2.5 \
  --reports_dir /path/to/private/production_reports
```

Preparation refuses an existing run root. Preserve the new run and Nextflow work
directory for `-resume`; never run two controllers against one output tree.
Logs and reports can contain private paths/identifiers and must not be committed.

## Execution Evidence

The supplied scientific results were generated using Bash/Slurm. The supplied
Nextflow records cover a 25-task stub run, local/Slurm completed-output audits
and a staging adapter check. They do not establish a fresh end-to-end scientific
Nextflow run or numerical equivalence of every output. No supplied workflow or
analysis was executed while preparing this public export. Its configuration now
requires explicit private paths and defaults to audit mode.
