#!/usr/bin/env bash
set -Eeuo pipefail

MUTECT2="${1:-/media/user/PDAC_SEQ_analysis/results/cosmic_annotation_mutect2_v103/mutect2_pass_cosmic_annotation.tsv}"
STRELKA="${2:-/media/user/PDAC_SEQ_analysis/results/cosmic_annotation_strelka_v103/strelka_pass_cosmic_annotation.tsv}"
OUT_DIR="${3:-/media/user/PDAC_SEQ_analysis/results/cosmic_annotation_caller_comparison_v103}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="$SCRIPT_DIR/compare_cosmic_callers.py"

python3 "$PY_SCRIPT" \
  --mutect2 "$MUTECT2" \
  --strelka "$STRELKA" \
  --output-dir "$OUT_DIR"
