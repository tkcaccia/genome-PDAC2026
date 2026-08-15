#!/usr/bin/env bash
set -Eeuo pipefail

BASE="${SEQ_BASE:-/media/user/SEQ}"
RESULTS_ROOT="/media/user/PDAC_SEQ_analysis/results"
WORK_ROOT="${SEQ_WORK:-/media/user/PDAC_SEQ_analysis/work}"
INTERVALS_BED="${PDAC_INTERVALS_BED:-$BASE/refs/optional/PDAC_Twist_ILMN_Exome_2.5_Plus_Panel.hg38.majority.bed}"
CHECK_INTERVAL="${PDAC_WATCHDOG_INTERVAL_SEC:-1800}"
STATE_LOG="$BASE/logs/pdac_watchdog.log"
PIDFILE="$BASE/tmp/pdac_watchdog.pid"
STATUSFILE="$BASE/tmp/pdac_watchdog.status"
RNASEQ_ARCHIVE="${PDAC_RNASEQ_ARCHIVE:-/home/user/PDAC_SEQ_archive/rnaseq_expression}"
RNAFUSION_OUT="${PDAC_RNAFUSION_OUT:-/home/user/PDAC_SEQ_native_results/rnafusion_pdac}"
RUN_ONCE=0

if [[ "${1:-}" == "--once" ]]; then
  RUN_ONCE=1
fi

mkdir -p "$BASE/logs" "$BASE/tmp"

if [[ -f "$BASE/configs/resources.env" ]]; then
  # shellcheck disable=SC1090
  source "$BASE/configs/resources.env"
fi

export PATH="/home/user/.local/bin:$PATH"

if [[ "$RUN_ONCE" -eq 0 ]]; then
  if [[ -f "$PIDFILE" ]]; then
    old_pid="$(cat "$PIDFILE" 2>/dev/null || true)"
    if [[ -n "$old_pid" ]] && kill -0 "$old_pid" 2>/dev/null; then
      echo "Watchdog already running with PID $old_pid"
      exit 0
    fi
  fi

  echo "$$" > "$PIDFILE"
  trap 'rm -f "$PIDFILE"' EXIT
fi

log() {
  printf '[%s] %s\n' "$(date '+%F %T %Z')" "$*" | tee -a "$STATE_LOG"
}

marker_exists() {
  local path
  for path in "$@"; do
    [[ -e "$path" ]] && return 0
  done
  return 1
}

write_status() {
  cat > "$STATUSFILE" <<EOF
timestamp=$(date '+%F %T %Z')
rnaseq_complete=$(rnaseq_complete && echo yes || echo no)
germline_complete=$(germline_complete && echo yes || echo no)
tumornormal_complete=$(tumornormal_complete && echo yes || echo no)
rnafusion_complete=$(rnafusion_complete && echo yes || echo no)
rnaseq_satisfied=$(rnaseq_satisfied && echo yes || echo no)
germline_satisfied=$(germline_satisfied && echo yes || echo no)
tumornormal_satisfied=$(tumornormal_satisfied && echo yes || echo no)
rnaseq_active=$(rnaseq_active && echo yes || echo no)
germline_active=$(germline_active && echo yes || echo no)
tumornormal_active=$(tumornormal_active && echo yes || echo no)
rnafusion_active=$(rnafusion_active && echo yes || echo no)
EOF
}

stage_success_log_exists() {
  local log_name="$1"
  local f
  while IFS= read -r f; do
    if grep -Eq 'Pipeline completed successfully|-\[nf-core/.+\] Pipeline completed successfully-' "$f" 2>/dev/null; then
      return 0
    fi
  done < <(find "$BASE/logs" -type f -name "$log_name" 2>/dev/null | sort)
  return 1
}

rnaseq_complete() {
  stage_success_log_exists "rnaseq.log" || marker_exists \
    "$RNASEQ_ARCHIVE/multiqc/star_salmon/multiqc_report.html"
}

germline_complete() {
  stage_success_log_exists "sarek_germline.log" || marker_exists \
    "$RESULTS_ROOT/sarek_germline/multiqc/multiqc_report.html"
}

tumornormal_complete() {
  stage_success_log_exists "sarek_tumor_normal.log" || marker_exists \
    "$RESULTS_ROOT/sarek_tumor_normal/multiqc/multiqc_report.html"
}

rnafusion_complete() {
  stage_success_log_exists "rnafusion.log" || marker_exists \
    "$RNAFUSION_OUT/multiqc/multiqc_report.html"
}

rnaseq_active() {
  pgrep -af 'run_rnaseq\.sh|rnaseq_expression|rnaseq-3\.24\.0/main\.nf' >/dev/null
}

germline_active() {
  pgrep -af 'sarek_germline|run_sarek\.sh.*germline' >/dev/null
}

tumornormal_active() {
  pgrep -af 'sarek_tumor_normal|run_sarek\.sh.*tumor-normal' >/dev/null
}

rnafusion_active() {
  pgrep -af 'rnafusion_pdac|run_rnafusion\.sh|rnafusion-4\.1\.0/main\.nf' >/dev/null
}

rnaseq_satisfied() {
  rnaseq_complete || germline_complete || germline_active || tumornormal_complete || tumornormal_active || rnafusion_complete || rnafusion_active
}

germline_satisfied() {
  germline_complete || tumornormal_complete || tumornormal_active || rnafusion_complete || rnafusion_active
}

tumornormal_satisfied() {
  tumornormal_complete || rnafusion_complete || rnafusion_active
}

launch_in_logdir() {
  local stage="$1"
  shift
  local ts logdir
  ts="$(date +%Y%m%d_%H%M%S)"
  logdir="$BASE/logs/pdac_watchdog_${ts}_${stage}"
  mkdir -p "$logdir"
  cat > "$logdir/runner.sh" <<'RUN'
#!/usr/bin/env bash
set -Eeuo pipefail
export PATH="/home/user/.local/bin:$PATH"
source /media/user/SEQ/configs/resources.env
RUN
  python3 - "$logdir/runner.sh" "$@" <<'PY'
import shlex
import sys
path = sys.argv[1]
args = sys.argv[2:]
with open(path, "a", encoding="utf-8") as fh:
    fh.write(" ".join(shlex.quote(a) for a in args))
    fh.write("\n")
PY
  chmod +x "$logdir/runner.sh"
  nohup bash "$logdir/runner.sh" > "$logdir/${stage}.log" 2>&1 &
  local pid=$!
  log "Launched ${stage} in $logdir (pid=$pid)"
}

launch_rnaseq() {
  launch_in_logdir rnaseq \
    /media/user/SEQ/scripts/run_rnaseq.sh \
    --samplesheet /media/user/SEQ/samplesheets/rnaseq_samplesheet.PDAC_RNA_fastq_autodraft.csv \
    --outdir /media/user/PDAC_SEQ_analysis/results/rnaseq_expression \
    -resume
}

launch_germline() {
  launch_in_logdir sarek_germline \
    /media/user/SEQ/scripts/run_sarek.sh \
    --mode germline \
    --samplesheet /media/user/SEQ/samplesheets/sarek_samplesheet.PDAC_WES_fastq_autodraft.csv \
    --intervals "$INTERVALS_BED" \
    --outdir /media/user/PDAC_SEQ_analysis/results/sarek_germline \
    -resume
}

launch_tumornormal() {
  launch_in_logdir sarek_tumor_normal \
    /media/user/SEQ/scripts/run_sarek.sh \
    --mode tumor-normal \
    --samplesheet /media/user/SEQ/samplesheets/sarek_samplesheet.PDAC_WES_fastq_autodraft.csv \
    --intervals "$INTERVALS_BED" \
    --outdir /media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal \
    -resume
}

launch_rnafusion() {
  launch_in_logdir rnafusion \
    /media/user/SEQ/scripts/run_rnafusion.sh \
    --samplesheet /media/user/SEQ/samplesheets/rnafusion_samplesheet.PDAC_RNA_fastq_autodraft.csv \
    --tools arriba,fusioncatcher,salmon \
    --outdir /media/user/SEQ/results/rnafusion_pdac \
    -resume
}

run_check() {
  write_status

  if ! rnaseq_satisfied; then
    if rnaseq_active; then
      log "rnaseq still running"
    else
      log "rnaseq incomplete and inactive; relaunching"
      launch_rnaseq
    fi
  elif ! germline_satisfied; then
    if germline_active; then
      log "sarek germline still running"
    else
      log "sarek germline incomplete and inactive; relaunching"
      launch_germline
    fi
  elif ! tumornormal_satisfied; then
    if tumornormal_active; then
      log "sarek tumor-normal still running"
    else
      log "sarek tumor-normal incomplete and inactive; relaunching"
      launch_tumornormal
    fi
  elif ! rnafusion_complete; then
    if rnafusion_active; then
      log "rnafusion still running"
    else
      log "rnafusion incomplete and inactive; relaunching"
      launch_rnafusion
    fi
  else
    log "All pancreatic stages are complete; exiting watchdog"
    exit 0
  fi

  write_status
}

if [[ "$RUN_ONCE" -eq 1 ]]; then
  log "PDAC watchdog one-shot check"
  run_check
  exit 0
fi

log "PDAC watchdog started (interval=${CHECK_INTERVAL}s)"

while true; do
  run_check
  sleep "$CHECK_INTERVAL"
done
