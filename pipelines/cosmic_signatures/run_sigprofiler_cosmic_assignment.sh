#!/usr/bin/env bash
set -Eeuo pipefail

BASE="/media/user/SEQ"
RESULTS_ROOT="/media/user/PDAC_SEQ_analysis/results"
SOMATIC_ROOT="$RESULTS_ROOT/sarek_tumor_normal/variant_calling/mutect2"
RUN_ROOT="${1:-$RESULTS_ROOT/cosmic_signatures_sigprofiler_assignment_1_1_3}"
PASS_VCF_DIR="$RUN_ROOT/input/pass_mutect2_vcfs_plain"
VENV_DIR="$BASE/venvs/sigprofilerassignment-1.1.3"
LOG_DIR="$RUN_ROOT/logs"
LOG_FILE="$LOG_DIR/run.log"
PYTHON_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$PYTHON_SCRIPT_DIR/run_sigprofiler_cosmic_assignment.py"
CPU="${SIGPROFILER_CPU:-8}"
CONTEXTS="${SIGPROFILER_CONTEXTS:-96 DINUC ID}"
SKIP_COMPLETED="${SIGPROFILER_SKIP_COMPLETED:-1}"
MAKE_PLOTS="${SIGPROFILER_MAKE_PLOTS:-0}"

mkdir -p "$PASS_VCF_DIR" "$LOG_DIR" "$BASE/venvs"
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
require_cmd pip3
require_cmd bcftools

if [[ ! -d "$VENV_DIR" ]]; then
  log "Creating virtual environment at $VENV_DIR"
  python3 -m venv "$VENV_DIR"
fi

source "$VENV_DIR/bin/activate"

log "Installing latest stable SigProfilerAssignment from PyPI"
python -m pip install --upgrade pip >>"$LOG_FILE" 2>&1
python -m pip install "SigProfilerAssignment==1.1.3" >>"$LOG_FILE" 2>&1

log "Installing SigProfilerMatrixGenerator GRCh38 reference"
python - <<'PY' >>"$LOG_FILE" 2>&1
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install("GRCh38")
PY

log "Preparing PASS-only Mutect2 VCFs"
find "$SOMATIC_ROOT" -type f -name '*.mutect2.filtered.vcf.gz' | sort | while read -r vcf; do
  sample="$(basename "$vcf" .mutect2.filtered.vcf.gz)"
  out="$PASS_VCF_DIR/${sample}.pass.vcf"
  if [[ -s "$out" ]]; then
    log "Keeping existing PASS VCF: $out"
    continue
  fi
  log "Creating PASS VCF for $sample"
  bcftools view -f PASS -Ov -o "$out" "$vcf" >>"$LOG_FILE" 2>&1
done

log "Launching SigProfilerAssignment on PASS-only Mutect2 VCFs"
args=(
  --samples "$PASS_VCF_DIR"
  --output "$RUN_ROOT/output"
  --genome-build GRCh38
  --cosmic-version 3.5
  --cpu "$CPU"
  --contexts $CONTEXTS
)

if [[ "$SKIP_COMPLETED" == "1" ]]; then
  args+=(--skip-completed)
fi

if [[ "$MAKE_PLOTS" == "1" ]]; then
  args+=(--make-plots)
fi

python "$PYTHON_SCRIPT" "${args[@]}" >>"$LOG_FILE" 2>&1

log "SigProfilerAssignment finished"
