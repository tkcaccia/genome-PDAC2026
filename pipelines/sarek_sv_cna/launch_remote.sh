#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SSH_HELPER="${PDAC2026_SSH_HELPER:?Set PDAC2026_SSH_HELPER to your private SSH wrapper}"

exec "$SSH_HELPER" script "$SCRIPT_DIR/run_sarek_sv_cna_remote.sh"
