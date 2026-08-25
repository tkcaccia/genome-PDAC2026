#!/usr/bin/env bash
set -Eeuo pipefail

RUN_ROOT="${1:-${SIGPROFILER_RUN_ROOT:-}}"
if [[ -z "$RUN_ROOT" ]]; then
  echo "Usage: $0 RUN_ROOT [OUT_DIR], or set SIGPROFILER_RUN_ROOT" >&2
  exit 2
fi
OUT_DIR="${2:-$RUN_ROOT/summary}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="$SCRIPT_DIR/summarize_sigprofiler_assignment.py"

python3 "$PY_SCRIPT" \
  --input-root "$RUN_ROOT/output" \
  --output-dir "$OUT_DIR"
