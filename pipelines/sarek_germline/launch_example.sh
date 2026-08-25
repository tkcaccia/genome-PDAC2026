#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
: "${SAREK_SAMPLESHEET:?Set SAREK_SAMPLESHEET}"
: "${SAREK_INTERVALS:?Set SAREK_INTERVALS}"
: "${SAREK_GERMLINE_OUTDIR:?Set SAREK_GERMLINE_OUTDIR outside this repository}"

"$SCRIPT_DIR/run_sarek.sh" \
  --mode germline \
  --samplesheet "$SAREK_SAMPLESHEET" \
  --intervals "$SAREK_INTERVALS" \
  --outdir "$SAREK_GERMLINE_OUTDIR" \
  -resume
