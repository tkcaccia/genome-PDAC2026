#!/usr/bin/env bash
set -Eeuo pipefail

# Replace these example paths with paths on the analysis host. Keep the private
# samplesheet, data, work, results, and logs outside the Git repository.
export RNAFUSION_PIPELINE_ROOT=/large-disk/code/nf-core-rnafusion-4.1.0
export RNAFUSION_REFERENCES=/large-disk/references/rnafusion-4.1.0
export RNAFUSION_REFERENCE_MANIFEST=/large-disk/manifests/rnafusion_reference_manifest.tsv
export RNAFUSION_SAMPLESHEET=/private-config/rnafusion_samplesheet.csv
export RNAFUSION_RESULTS_ROOT=/large-disk/results/rnafusion-4.1.0
export RNAFUSION_WORK_ROOT=/large-disk/work/rnafusion-4.1.0
export RNAFUSION_STATE_ROOT=/large-disk/state/rnafusion-4.1.0
export RNAFUSION_LOG_ROOT=/large-disk/logs/rnafusion-4.1.0

# This workstation profile waits until 24 GiB RAM and 150 GiB disk are free,
# runs one sample/task at a time, and invokes all four selected methods.
export RNAFUSION_MIN_AVAILABLE_GIB=24
export RNAFUSION_MIN_DISK_GIB=150
export RNAFUSION_TOOLS=arriba,fusioncatcher,starfusion,salmon

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
"$SCRIPT_DIR/run_sample_queue.sh"
