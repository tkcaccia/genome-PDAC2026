#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SSH_HELPER="/Users/stefano/Documents/SEQ/scripts/ssh_pdac_remote.expect"

exec "$SSH_HELPER" script "$SCRIPT_DIR/run_sarek_sv_cna_remote.sh"
