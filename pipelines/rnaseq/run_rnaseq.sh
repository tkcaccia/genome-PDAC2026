#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
seq_source_env

SAMPLESHEET=""
OUTDIR=""
WORKDIR=""
EXTRA_ARGS=()
RESUME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --work-dir) WORKDIR="$2"; shift 2 ;;
    -resume|--resume) RESUME="-resume"; shift ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

[[ -n "$SAMPLESHEET" ]] || seq_die "Missing --samplesheet"
seq_require_file "$SAMPLESHEET"
python3 "$SCRIPT_DIR/validate_samplesheet.py" --pipeline rnaseq --input "$SAMPLESHEET"

TMP_SHEET="$(mktemp "$SEQ_TMP/rnaseq_samplesheet.XXXXXX.csv")"
seq_strip_csv_comments "$SAMPLESHEET" "$TMP_SHEET"

OUTDIR="${OUTDIR:-${SEQ_RNASEQ_OUTDIR:-$SEQ_RESULTS/rnaseq}}"
WORKDIR="${WORKDIR:-${SEQ_RNASEQ_WORKDIR:-$SEQ_WORK}}"
PROFILE="$(seq_runtime_profile)"
mkdir -p "$OUTDIR" "$WORKDIR"
nextflow run "$SEQ_BASE/pipelines/rnaseq-3.24.0/main.nf" \
  -profile "$PROFILE" \
  -c "$SEQ_BASE/configs/rnaseq.config" \
  --input "$TMP_SHEET" \
  --outdir "$OUTDIR" \
  --fasta "$SEQ_RNASEQ_FASTA" \
  --gtf "$SEQ_RNASEQ_GTF" \
  --star_index "$SEQ_RNASEQ_STAR_INDEX" \
  --salmon_index "$SEQ_RNASEQ_SALMON_INDEX" \
  --igenomes_ignore \
  -w "$WORKDIR" \
  $RESUME \
  "${EXTRA_ARGS[@]}"
