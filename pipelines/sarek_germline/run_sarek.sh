#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
seq_source_env

MODE=""
SAMPLESHEET=""
INTERVALS=""
OUTDIR=""
WORKDIR=""
EXTRA_ARGS=()
RESUME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --intervals) INTERVALS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --work-dir) WORKDIR="$2"; shift 2 ;;
    -resume|--resume) RESUME="-resume"; shift ;;
    *) EXTRA_ARGS+=("$1"); shift ;;
  esac
done

[[ -n "$MODE" ]] || seq_die "Missing --mode (tumor-normal or germline)"
[[ -n "$SAMPLESHEET" ]] || seq_die "Missing --samplesheet"
[[ -n "$INTERVALS" ]] || seq_die "Missing --intervals with the exome capture BED / interval list"

seq_require_file "$SAMPLESHEET"
seq_require_file "$INTERVALS"
python3 "$SCRIPT_DIR/validate_samplesheet.py" --pipeline sarek --input "$SAMPLESHEET"

TMP_SHEET="$(mktemp "$SEQ_TMP/sarek_samplesheet.XXXXXX.csv")"
seq_strip_csv_comments "$SAMPLESHEET" "$TMP_SHEET"
python3 - "$TMP_SHEET" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
rows = []
with path.open(newline="") as handle:
    reader = csv.DictReader(handle)
    fieldnames = reader.fieldnames
    if not fieldnames:
        raise SystemExit("Sarek samplesheet is missing a header row")
    for row in reader:
        patient = (row.get("patient") or "").strip()
        if patient and patient.isdigit():
            row["patient"] = f"P{patient}"
        rows.append(row)

with path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
PY

if [[ "$INTERVALS" == *.bed ]]; then
  SORTED_INTERVALS="$(mktemp "$SEQ_TMP/sarek_intervals.XXXXXX.bed")"
  awk 'BEGIN { OFS="\t" } NF >= 3 && $0 !~ /^[[:space:]]*#/ && $0 !~ /^track([[:space:]]|$)/ && $0 !~ /^browser([[:space:]]|$)/ { print $1, $2, $3 }' "$INTERVALS" \
    | LC_ALL=C sort -k1,1V -k2,2n -k3,3n > "$SORTED_INTERVALS"
  INTERVALS="$SORTED_INTERVALS"
fi

case "$MODE" in
  tumor-normal)
    TOOLS="mutect2,strelka,vep"
    ;;
  germline)
    TOOLS="haplotypecaller,vep"
    ;;
  *)
    seq_die "Unsupported mode: $MODE"
    ;;
esac

PON_ARGS=()
if [[ -n "${SEQ_SAREK_PON:-}" && -s "${SEQ_SAREK_PON:-}" ]]; then
  PON_ARGS+=(--pon "$SEQ_SAREK_PON")
fi
if [[ -n "${SEQ_SAREK_PON_TBI:-}" && -s "${SEQ_SAREK_PON_TBI:-}" ]]; then
  PON_ARGS+=(--pon_tbi "$SEQ_SAREK_PON_TBI")
fi
if [[ "$MODE" == "tumor-normal" && "${#PON_ARGS[@]}" -eq 0 ]]; then
  seq_log "WARNING: No controlled MuTect2 Panel of Normals is configured; continuing without --pon."
fi

OUTDIR="${OUTDIR:-${SEQ_SAREK_OUTDIR:-$SEQ_RESULTS/sarek}}"
WORKDIR="${WORKDIR:-${SEQ_SAREK_WORKDIR:-$SEQ_WORK}}"
PROFILE="$(seq_runtime_profile)"
LOCAL_IGENOMES_BASE="${SEQ_SAREK_IGENOMES_BASE:-$SEQ_BASE/refs/igenomes_stub}"
LOCAL_SNPEFF_CACHE="${SEQ_SAREK_SNPEFF_CACHE:-$SEQ_BASE/refs/annotation/snpeff_cache}"
mkdir -p "$OUTDIR" "$WORKDIR" "$LOCAL_IGENOMES_BASE" "$LOCAL_SNPEFF_CACHE"
nextflow run "$SEQ_BASE/pipelines/sarek-3.8.1/main.nf" \
  -profile "$PROFILE" \
  -c "$SEQ_BASE/configs/sarek.config" \
  --input "$TMP_SHEET" \
  --outdir "$OUTDIR" \
  --tools "$TOOLS" \
  --fasta "$SEQ_SAREK_FASTA" \
  --fasta_fai "$SEQ_SAREK_FASTA_FAI" \
  --dict "$SEQ_SAREK_DICT" \
  --dbsnp "$SEQ_SAREK_DBSNP" \
  --dbsnp_tbi "$SEQ_SAREK_DBSNP_TBI" \
  --known_indels "$SEQ_SAREK_KNOWN_INDELS" \
  --known_indels_tbi "$SEQ_SAREK_KNOWN_INDELS_TBI" \
  --known_snps "$SEQ_SAREK_KNOWN_SNPS" \
  --known_snps_tbi "$SEQ_SAREK_KNOWN_SNPS_TBI" \
  --germline_resource "$SEQ_SAREK_GERMLINE_RESOURCE" \
  --germline_resource_tbi "$SEQ_SAREK_GERMLINE_RESOURCE_TBI" \
  "${PON_ARGS[@]}" \
  --igenomes_base "$LOCAL_IGENOMES_BASE" \
  --snpeff_cache "$LOCAL_SNPEFF_CACHE" \
  --vep_cache "$SEQ_VEP_CACHE_ROOT" \
  --vep_cache_version "$SEQ_VEP_CACHE_VERSION" \
  --vep_genome GRCh38 \
  --vep_species homo_sapiens \
  --intervals "$INTERVALS" \
  --igenomes_ignore \
  -w "$WORKDIR" \
  $RESUME \
  "${EXTRA_ARGS[@]}"
