#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
: "${RNASEQ_SAMPLESHEET:?Set RNASEQ_SAMPLESHEET}"
: "${RNASEQ_OUTDIR:?Set RNASEQ_OUTDIR outside this repository}"

"$SCRIPT_DIR/run_rnaseq.sh" \
  --samplesheet "$RNASEQ_SAMPLESHEET" \
  --outdir "$RNASEQ_OUTDIR" \
  -resume
