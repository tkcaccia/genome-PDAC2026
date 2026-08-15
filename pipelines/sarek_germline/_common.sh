#!/usr/bin/env bash
set -Eeuo pipefail

SEQ_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SEQ_BASE="$(cd -- "$SEQ_SCRIPT_DIR/.." && pwd)"
SEQ_TMP="$SEQ_BASE/tmp"
SEQ_WORK="${SEQ_WORK:-$SEQ_BASE/work}"
SEQ_RESULTS="${SEQ_RESULTS:-$SEQ_BASE/results}"
SEQ_CONFIG="${SEQ_CONFIG:-$SEQ_BASE/configs/resources.env}"
SEQ_BACKUPS="$SEQ_BASE/tmp/script_backups"

# Prefer the user-local Nextflow install when scripts are run from a non-login shell.
if [[ -d "$HOME/.local/bin" ]]; then
  export PATH="$HOME/.local/bin:$PATH"
fi

# Surface workspace-pinned helper binaries such as bgzip and tabix.
if [[ -d "$SEQ_BASE/bin" ]]; then
  export PATH="$SEQ_BASE/bin:$PATH"
fi

seq_log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

seq_die() {
  seq_log "ERROR: $*" >&2
  exit 1
}

seq_require_cmd() {
  command -v "$1" >/dev/null 2>&1 || seq_die "Required command not found: $1"
}

seq_require_file() {
  [[ -s "$1" ]] || seq_die "Required file missing: $1"
}

seq_source_env() {
  seq_require_file "$SEQ_CONFIG"
  # shellcheck disable=SC1090
  source "$SEQ_CONFIG"
  export NXF_HOME NXF_WORK NXF_SINGULARITY_CACHEDIR NXF_APPTAINER_CACHEDIR SINGULARITY_CACHEDIR APPTAINER_CACHEDIR
}

seq_backup_if_exists() {
  local target="$1"
  if [[ -e "$target" || -L "$target" ]]; then
    local stamp
    stamp="$(date +%Y%m%d_%H%M%S)"
    local rel="${target#/}"
    local dest="$SEQ_BACKUPS/$stamp/$rel"
    mkdir -p "$(dirname "$dest")"
    cp -a "$target" "$dest"
  fi
}

seq_strip_csv_comments() {
  local input="$1"
  local output="$2"
  awk 'NF && $0 !~ /^[[:space:]]*#/' "$input" > "$output"
}

seq_runtime_profile() {
  if [[ -n "${SEQ_RUNTIME_PROFILE:-}" ]]; then
    printf '%s\n' "$SEQ_RUNTIME_PROFILE"
    return 0
  fi
  if command -v apptainer >/dev/null 2>&1; then
    printf 'apptainer\n'
  elif command -v singularity >/dev/null 2>&1; then
    printf 'singularity\n'
  elif command -v docker >/dev/null 2>&1; then
    printf 'docker\n'
  else
    seq_die "No supported container runtime detected."
  fi
}
