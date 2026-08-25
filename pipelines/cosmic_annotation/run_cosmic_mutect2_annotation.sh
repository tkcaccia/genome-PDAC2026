#!/usr/bin/env bash
set -Eeuo pipefail

BASE="${PDAC2026_RUNTIME_ROOT:?Set PDAC2026_RUNTIME_ROOT to a restricted runtime directory}"
RESULTS_ROOT="${PDAC2026_RESULTS_ROOT:?Set PDAC2026_RESULTS_ROOT outside this repository}"
COSMIC_ARCHIVE_DIR="${1:-$BASE/refs/cosmic_archives_v103}"
COSMIC_EXTRACT_DIR="${2:-$BASE/refs/cosmic_v103}"
VCF_DIR="${MUTECT2_ANNOTATION_DIR:-$RESULTS_ROOT/sarek_tumor_normal/annotation/mutect2}"
OUT_DIR="${3:-$RESULTS_ROOT/cosmic_annotation_mutect2_v103}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="$SCRIPT_DIR/annotate_mutect2_with_cosmic.py"
LOG_DIR="$OUT_DIR/logs"
LOG_FILE="$LOG_DIR/run.log"

mkdir -p "$COSMIC_EXTRACT_DIR" "$OUT_DIR" "$LOG_DIR"
: > "$LOG_FILE"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG_FILE"
}

extract_one() {
  local archive="$1"
  local stem="$2"
  if [[ -s "$COSMIC_EXTRACT_DIR/${stem}.tsv.gz" ]]; then
    log "Keeping existing ${stem}.tsv.gz"
    return 0
  fi
  log "Extracting $archive"
  tar -xf "$COSMIC_ARCHIVE_DIR/$archive" -C "$COSMIC_EXTRACT_DIR" >>"$LOG_FILE" 2>&1
}

extract_one Cosmic_MutantCensus_Tsv_v103_GRCh38.tar Cosmic_MutantCensus_v103_GRCh38
extract_one Cosmic_CancerGeneCensus_Tsv_v103_GRCh38.tar Cosmic_CancerGeneCensus_v103_GRCh38
extract_one Cosmic_CancerGeneCensusHallmarksOfCancer_Tsv_v103_GRCh38.tar Cosmic_CancerGeneCensusHallmarksOfCancer_v103_GRCh38
extract_one Cosmic_ResistanceMutations_Tsv_v103_GRCh38.tar Cosmic_ResistanceMutations_v103_GRCh38

log "Running COSMIC annotation over Mutect2 PASS variants"
python3 "$PY_SCRIPT" \
  --vcf-dir "$VCF_DIR" \
  --cosmic-dir "$COSMIC_EXTRACT_DIR" \
  --output-dir "$OUT_DIR" >>"$LOG_FILE" 2>&1

log "COSMIC Mutect2 annotation finished"
