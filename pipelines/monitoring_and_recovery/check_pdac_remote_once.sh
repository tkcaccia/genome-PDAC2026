#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SSH_HELPER="${PDAC_WATCH_SSH_HELPER:-$SCRIPT_DIR/ssh_pdac_remote.expect}"
REMOTE_WATCHDOG_LOCAL="${PDAC_WATCH_REMOTE_SCRIPT:-$SCRIPT_DIR/pdac_remote_watchdog.sh}"
REMOTE_WATCHDOG_PATH="/media/user/SEQ/scripts/pdac_watchdog.sh"
LOCAL_LOG_DIR="${PDAC_WATCH_LOG_DIR:-$HOME/Library/Logs/PDACWatch}"
LOGFILE="$LOCAL_LOG_DIR/pdac_local_watchdog.log"

mkdir -p "$LOCAL_LOG_DIR"

log() {
  printf '[%s] %s\n' "$(date '+%F %T %Z')" "$*" | tee -a "$LOGFILE"
}

run_remote() {
  "$SSH_HELPER" "$1"
}

sync_remote_watchdog() {
  local installer
  installer="$(mktemp "${TMPDIR:-/tmp}/pdac_watchdog_sync.XXXXXX")"
  cat >"$installer" <<EOF
#!/usr/bin/env bash
set -Eeuo pipefail
mkdir -p /media/user/SEQ/scripts /media/user/SEQ/logs /media/user/SEQ/tmp
cat > '$REMOTE_WATCHDOG_PATH' <<'REMOTE_WATCHDOG_EOF'
$(cat "$REMOTE_WATCHDOG_LOCAL")
REMOTE_WATCHDOG_EOF
chmod +x '$REMOTE_WATCHDOG_PATH'
EOF
  chmod 700 "$installer"
  "$SSH_HELPER" script "$installer"
  rm -f "$installer"
}

main() {
  local output

  log "Starting 30-minute PDAC remote check"

  if ! sync_remote_watchdog >>"$LOGFILE" 2>&1; then
    log "Remote sync failed; will retry on the next scheduled run"
    exit 0
  fi

  if ! output="$(
    run_remote "pkill -f '/media/user/SEQ/scripts/pdac_watchdog.sh$' >/dev/null 2>&1 || true; rm -f /media/user/SEQ/tmp/pdac_watchdog.pid; bash '$REMOTE_WATCHDOG_PATH' --once; echo ---STATUS---; cat /media/user/SEQ/tmp/pdac_watchdog.status 2>/dev/null || true; echo ---PGREP---; pgrep -af 'run_rnaseq|run_sarek|run_rnafusion|nextflow|fusioncatcher|arriba|java' || true" 2>&1
  )"; then
    log "Remote health check failed; will retry on the next scheduled run"
    printf '%s\n' "$output" >>"$LOGFILE"
    exit 0
  fi

  log "Remote one-shot check completed"
  printf '%s\n' "$output" >>"$LOGFILE"
}

main "$@"
