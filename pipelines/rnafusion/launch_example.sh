#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
: "${RNAFUSION_SAMPLESHEET:?Set RNAFUSION_SAMPLESHEET}"
: "${RNAFUSION_OUTDIR:?Set RNAFUSION_OUTDIR outside this repository}"

"$SCRIPT_DIR/run_rnafusion.sh" \
  --samplesheet "$RNAFUSION_SAMPLESHEET" \
  --tools arriba,fusioncatcher,salmon \
  --outdir "$RNAFUSION_OUTDIR" \
  -resume
