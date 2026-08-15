#!/usr/bin/env bash
set -Eeuo pipefail

RUN_ROOT="${1:-/media/user/PDAC_SEQ_analysis/results/cosmic_signatures_sigprofiler_assignment_1_1_3}"
OUT_DIR="${2:-$RUN_ROOT/summary}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="$SCRIPT_DIR/summarize_sigprofiler_assignment.py"

python3 "$PY_SCRIPT" \
  --input-root "$RUN_ROOT/output" \
  --output-dir "$OUT_DIR"
