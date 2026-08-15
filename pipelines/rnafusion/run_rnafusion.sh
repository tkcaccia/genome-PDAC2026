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
TOOLS="${SEQ_RNAFUSION_DEFAULT_TOOLS:-arriba,fusioncatcher,salmon}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --tools) TOOLS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --work-dir) WORKDIR="$2"; shift 2 ;;
    -resume|--resume) RESUME="-resume"; shift ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

[[ -n "$SAMPLESHEET" ]] || seq_die "Missing --samplesheet"
seq_require_file "$SAMPLESHEET"
python3 "$SCRIPT_DIR/validate_samplesheet.py" --pipeline rnafusion --input "$SAMPLESHEET"

TMP_SHEET="$(mktemp "$SEQ_TMP/rnafusion_samplesheet.XXXXXX.csv")"
seq_strip_csv_comments "$SAMPLESHEET" "$TMP_SHEET"

if [[ ! -d "$SEQ_RNAFUSION_GENOMES_BASE/GRCh38/gencode_v46/star" ]]; then
  seq_die "RNA fusion STAR references are missing. Re-run sync_rnafusion_refs.sh first."
fi
if [[ ",$TOOLS," == *",fusioncatcher,"* ]] && [[ ! -d "$SEQ_RNAFUSION_GENOMES_BASE/GRCh38/gencode_v46/fusioncatcher" ]]; then
  seq_die "FusionCatcher was requested but its references are missing. Re-run sync_rnafusion_refs.sh with FusionCatcher enabled first."
fi
if [[ ",$TOOLS," == *",starfusion,"* ]] && [[ ! -d "$SEQ_RNAFUSION_GENOMES_BASE/GRCh38/gencode_v46/starfusion" ]]; then
  seq_die "STAR-Fusion was requested but its references are missing. Re-run sync_rnafusion_refs.sh with STAR-Fusion enabled first."
fi

OUTDIR="${OUTDIR:-${SEQ_RNAFUSION_OUTDIR:-$SEQ_RESULTS/rnafusion}}"
WORKDIR="${WORKDIR:-${SEQ_RNAFUSION_WORKDIR:-$SEQ_WORK}}"
PROFILE="$(seq_runtime_profile)"
mkdir -p "$OUTDIR" "$WORKDIR"
nextflow run "$SEQ_BASE/pipelines/rnafusion-4.1.0/main.nf" \
  -profile "$PROFILE" \
  -c "$SEQ_BASE/configs/rnafusion.config" \
  --input "$TMP_SHEET" \
  --outdir "$OUTDIR" \
  --genomes_base "$SEQ_RNAFUSION_GENOMES_BASE" \
  --genome GRCh38 \
  --genome_gencode_version 46 \
  --tools "$TOOLS" \
  --no_cosmic \
  -w "$WORKDIR" \
  $RESUME \
  "${EXTRA_ARGS[@]}"
