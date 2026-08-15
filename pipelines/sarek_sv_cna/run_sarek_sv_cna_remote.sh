#!/usr/bin/env bash
set -Eeuo pipefail

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
  log "ERROR: $*" >&2
  exit 1
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

download_if_missing() {
  local url="$1"
  local dest="$2"
  if [[ -s "$dest" ]]; then
    log "Reusing $(basename "$dest")"
    return 0
  fi
  mkdir -p "$(dirname "$dest")"
  log "Downloading $(basename "$dest")"
  curl -L --fail --retry 5 --retry-delay 5 -o "$dest" "$url"
}

prepare_ascat_bundle_zip() {
  local input_zip="$1"
  local output_zip="$2"
  local tmp_root="$3"
  local target_prefix="$4"
  local add_chr_prefix="${5:-false}"

  local workdir extracted target_dir
  workdir="$(mktemp -d "$tmp_root/ascat_bundle.XXXXXX")"
  trap 'rm -rf "$workdir"' RETURN

  unzip -q -o "$input_zip" -d "$workdir"
  extracted="$(find "$workdir" -mindepth 1 -maxdepth 1 -type d | head -n 1)"
  [[ -n "$extracted" ]] || die "Could not find extracted directory in $(basename "$input_zip")"

  target_dir="$workdir/$target_prefix"
  mv "$extracted" "$target_dir"

  while IFS= read -r f; do
    [[ -s "$f" ]] || continue
    local base chrom newbase first_line
    base="$(basename "$f")"
    chrom="${base##*_chr}"
    chrom="${chrom%.txt}"
    newbase="${target_prefix}_chr${chrom}.txt"
    if [[ "$base" != "$newbase" ]]; then
      mv "$f" "$(dirname "$f")/$newbase"
      f="$(dirname "$f")/$newbase"
    fi
    if [[ "$add_chr_prefix" == "true" ]]; then
      first_line="$(head -n 1 "$f")"
      case "$first_line" in
        chr*) ;;
        *) sed -i 's/^/chr/' "$f" ;;
      esac
    fi
  done < <(find "$target_dir" -type f -name '*_chr*.txt' | sort)

  mkdir -p "$(dirname "$output_zip")"
  (
    rm -f "$output_zip"
    cd "$target_dir"
    zip -qj "$output_zip" ./*.txt
  )
  log "Prepared $(basename "$output_zip")"
}

run_pipeline() {
  need_cmd python3
  need_cmd awk
  need_cmd sort
  need_cmd curl
  need_cmd unzip
  need_cmd zip
  need_cmd nextflow

  source /media/user/SEQ/configs/resources.env

  local seq_base="${SEQ_BASE:-/media/user/SEQ}"
  local seq_tmp="$seq_base/tmp"
  local samplesheet="/media/user/SEQ/samplesheets/sarek_samplesheet.PDAC_WES_fastq_autodraft.csv"
  local input_override="${PDAC_SV_CNA_INPUT:-/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal_sv_cna/csv/recalibrated_restart_full.csv}"
  local start_step="${PDAC_SV_CNA_STEP:-variant_calling}"
  local intervals="/media/user/SEQ/refs/optional/PDAC_Twist_ILMN_Exome_2.5_Plus_Panel.hg38.majority.bed"
  local outdir="/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal_sv_cna"
  local workdir="/media/user/PDAC_SEQ_analysis/work"
  local ref_base="/media/user/New_Volume3/Lion/PDAC/SEQ_refs/ascat_wes_hg38"
  local profile="${SEQ_RUNTIME_PROFILE:-singularity}"
  local local_igenomes_base="${SEQ_SAREK_IGENOMES_BASE:-$seq_base/refs/igenomes_stub}"
  local local_snpeff_cache="${SEQ_SAREK_SNPEFF_CACHE:-$seq_base/refs/annotation/snpeff_cache}"
  local logdir="${LOGDIR:-$seq_base/logs/pdac_sarek_sv_cna_$(date +%Y%m%d_%H%M%S)}"

  mkdir -p "$seq_tmp" "$outdir" "$workdir" "$local_igenomes_base" "$local_snpeff_cache" "$ref_base/downloads" "$ref_base/processed" "$logdir"

  local tmp_sheet input_sheet
  input_sheet=""
  if [[ -n "$input_override" ]]; then
    input_sheet="$input_override"
  else
    tmp_sheet="$(mktemp "$seq_tmp/sarek_sv_cna_samplesheet.XXXXXX.csv")"
    awk 'NF && $0 !~ /^[[:space:]]*#/' "$samplesheet" > "$tmp_sheet"
    python3 - "$tmp_sheet" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
with path.open(newline="") as handle:
    reader = csv.DictReader(handle)
    fieldnames = reader.fieldnames
    if not fieldnames:
        raise SystemExit("Sarek samplesheet is missing a header row")
    rows = []
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
    input_sheet="$tmp_sheet"
  fi

  local sorted_intervals="$intervals"
  if [[ "$intervals" == *.bed ]]; then
    sorted_intervals="$(mktemp "$seq_tmp/sarek_sv_cna_intervals.XXXXXX.bed")"
    awk 'BEGIN { OFS="\t" } NF >= 3 && $0 !~ /^[[:space:]]*#/ && $0 !~ /^track([[:space:]]|$)/ && $0 !~ /^browser([[:space:]]|$)/ { print $1, $2, $3 }' "$intervals" \
      | LC_ALL=C sort -k1,1V -k2,2n -k3,3n > "$sorted_intervals"
  fi

  local loci_zip="$ref_base/downloads/G1000_loci_WES_hg38.zip"
  local alleles_zip="$ref_base/downloads/G1000_alleles_WES_hg38.zip"
  local gc_zip="$ref_base/downloads/GC_G1000_WES_hg38.zip"
  local rt_zip="$ref_base/downloads/RT_G1000_WES_hg38.zip"
  local alleles_processed_zip="$ref_base/processed/G1000_alleles_WES_hg38.zip"
  local loci_processed_zip="$ref_base/processed/G1000_loci_WES_hg38.zip"

  download_if_missing "https://zenodo.org/records/14008443/files/G1000_loci_WES_hg38.zip?download=1" "$loci_zip"
  download_if_missing "https://zenodo.org/records/14008443/files/G1000_alleles_WES_hg38.zip?download=1" "$alleles_zip"
  download_if_missing "https://zenodo.org/records/14008443/files/GC_G1000_WES_hg38.zip?download=1" "$gc_zip"
  download_if_missing "https://zenodo.org/records/14008443/files/RT_G1000_WES_hg38.zip?download=1" "$rt_zip"
  prepare_ascat_bundle_zip "$alleles_zip" "$alleles_processed_zip" "$seq_tmp" "G1000_alleles_WES_hg38" "false"
  prepare_ascat_bundle_zip "$loci_zip" "$loci_processed_zip" "$seq_tmp" "G1000_loci_WES_hg38" "true"

  {
    log "Launching Sarek SV/CNA branch"
    log "Input: $input_sheet"
    log "Step: ${start_step:-full}"
    log "Intervals: $intervals"
    log "Outdir: $outdir"
    log "Workdir: $workdir"
    log "ASCAT refs: $ref_base"
  } | tee -a "$logdir/run.log"

  local -a step_args=()
  if [[ -n "$start_step" ]]; then
    step_args+=( --step "$start_step" )
  fi

  nextflow run "$seq_base/pipelines/sarek-3.8.1/main.nf" \
    -profile "$profile" \
    -c "$seq_base/configs/sarek.config" \
    --input "$input_sheet" \
    --outdir "$outdir" \
    --tools "manta,ascat" \
    --save_output_as_bam \
    "${step_args[@]}" \
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
    --ascat_alleles "$alleles_processed_zip" \
    --ascat_loci "$loci_processed_zip" \
    --ascat_loci_gc "$gc_zip" \
    --ascat_loci_rt "$rt_zip" \
    --ascat_genome hg38 \
    --igenomes_base "$local_igenomes_base" \
    --snpeff_cache "$local_snpeff_cache" \
    --vep_cache "$SEQ_VEP_CACHE_ROOT" \
    --vep_cache_version "$SEQ_VEP_CACHE_VERSION" \
    --genome GRCh38 \
    --vep_genome GRCh38 \
    --vep_species homo_sapiens \
    --intervals "$sorted_intervals" \
    --igenomes_ignore \
    -w "$workdir" \
    -resume 2>&1 | tee -a "$logdir/run.log"
}

if [[ "${1:-}" == "--foreground" ]]; then
  run_pipeline
  exit 0
fi

source /media/user/SEQ/configs/resources.env
remote_logdir="${LOGDIR:-/media/user/SEQ/logs/pdac_sarek_sv_cna_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$remote_logdir"

nohup bash -lc "export LOGDIR='$remote_logdir'; $(declare -f log); $(declare -f die); $(declare -f need_cmd); $(declare -f download_if_missing); $(declare -f prepare_ascat_bundle_zip); $(declare -f run_pipeline); run_pipeline" \
  > "$remote_logdir/nohup.out" 2>&1 &

log "Launched Sarek SV/CNA in background"
log "PID: $!"
log "Log directory: $remote_logdir"
