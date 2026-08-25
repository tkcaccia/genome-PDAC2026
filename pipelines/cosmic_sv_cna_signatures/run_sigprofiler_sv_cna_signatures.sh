#!/usr/bin/env bash
set -Eeuo pipefail

BASE="${PDAC2026_RUNTIME_ROOT:?Set PDAC2026_RUNTIME_ROOT to a restricted runtime directory}"
RESULTS_ROOT="${PDAC2026_RESULTS_ROOT:?Set PDAC2026_RESULTS_ROOT outside this repository}"
SAREK_ROOT="${SAREK_SV_CNA_ROOT:-$RESULTS_ROOT/sarek_tumor_normal_sv_cna}"
RUN_ROOT="${1:-$RESULTS_ROOT/cosmic_sv_cna_signatures_sigprofiler_assignment_1_1_3}"
VENV_DIR="$BASE/venvs/sigprofilerassignment-1.1.3"
LOG_DIR="$RUN_ROOT/logs"
LOG_FILE="$LOG_DIR/run.log"
if [[ -n "${BASH_SOURCE[0]:-}" ]]; then
  PYTHON_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
elif [[ -n "${SIGPROFILER_SCRIPT_DIR:-}" ]]; then
  PYTHON_SCRIPT_DIR="$SIGPROFILER_SCRIPT_DIR"
else
  PYTHON_SCRIPT_DIR="$(pwd)"
fi
PYTHON_SCRIPT="$PYTHON_SCRIPT_DIR/run_sigprofiler_sv_cna_signatures.py"
SUMMARY_SCRIPT="$PYTHON_SCRIPT_DIR/summarize_sigprofiler_sv_cna.py"
CPU="${SIGPROFILER_CPU:-8}"
SKIP_COMPLETED="${SIGPROFILER_SKIP_COMPLETED:-1}"
MAKE_PLOTS="${SIGPROFILER_MAKE_PLOTS:-0}"
PROJECT="${SIGPROFILER_PROJECT:-PDAC2026}"

mkdir -p "$RUN_ROOT" "$LOG_DIR" "$BASE/venvs"
: > "$LOG_FILE"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG_FILE"
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    log "ERROR: missing command: $1"
    exit 1
  }
}

require_cmd python3

if [[ ! -d "$VENV_DIR" ]]; then
  log "Creating virtual environment at $VENV_DIR"
  python3 -m venv "$VENV_DIR"
fi

source "$VENV_DIR/bin/activate"

log "Installing SigProfilerAssignment and dependencies"
python -m pip install --upgrade pip >>"$LOG_FILE" 2>&1
python -m pip install "SigProfilerAssignment==1.1.3" >>"$LOG_FILE" 2>&1

args=(
  --sarek-root "$SAREK_ROOT"
  --output "$RUN_ROOT/output"
  --project "$PROJECT"
  --genome-build GRCh38
  --cosmic-version 3.5
  --cpu "$CPU"
)

if [[ "$SKIP_COMPLETED" == "1" ]]; then
  args+=(--skip-completed)
fi

if [[ "$MAKE_PLOTS" == "1" ]]; then
  args+=(--make-plots)
fi

log "Launching SigProfiler CNV48 and SV32 assignment"
python "$PYTHON_SCRIPT" "${args[@]}" >>"$LOG_FILE" 2>&1

log "Summarizing CNV48 and SV32 activities"
python "$SUMMARY_SCRIPT" \
  --input-root "$RUN_ROOT/output" \
  --output-dir "$RUN_ROOT/summary" >>"$LOG_FILE" 2>&1

log "SigProfiler CNV48/SV32 assignment finished"
