#!/usr/bin/env bash
set -Eeuo pipefail

MUTECT2="${1:-${COSMIC_MUTECT2_TABLE:-}}"
STRELKA="${2:-${COSMIC_STRELKA_TABLE:-}}"
OUT_DIR="${3:-${COSMIC_COMPARISON_OUTDIR:-}}"
if [[ -z "$MUTECT2" || -z "$STRELKA" || -z "$OUT_DIR" ]]; then
  echo "Usage: $0 MUTECT2_TABLE STRELKA_TABLE OUT_DIR" >&2
  exit 2
fi
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="$SCRIPT_DIR/compare_cosmic_callers.py"

python3 "$PY_SCRIPT" \
  --mutect2 "$MUTECT2" \
  --strelka "$STRELKA" \
  --output-dir "$OUT_DIR"
