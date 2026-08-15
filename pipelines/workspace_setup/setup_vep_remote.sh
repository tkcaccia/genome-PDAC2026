#!/usr/bin/env bash
set -Eeuo pipefail

source /media/user/SEQ/scripts/_common.sh
seq_source_env

OUT="/media/user/SEQ/refs/vep"
LOG_FILE="/media/user/SEQ/logs/setup_vep.log"
VEP_SIF="/media/user/SEQ/containers/ensembl-vep_release_115.2.sif"
VEP_IMAGE="docker://ensemblorg/ensembl-vep:release_115.2"
: > "$LOG_FILE"

mkdir -p "$OUT/cache" "$OUT/fasta"

if [[ ! -s "$VEP_SIF" ]]; then
  if command -v singularity >/dev/null 2>&1; then
    singularity pull "$VEP_SIF" "$VEP_IMAGE" >>"$LOG_FILE" 2>&1
  elif command -v apptainer >/dev/null 2>&1; then
    apptainer pull "$VEP_SIF" "$VEP_IMAGE" >>"$LOG_FILE" 2>&1
  else
    seq_die "No Apptainer/Singularity runtime available for VEP container pull."
  fi
fi

if [[ ! -d "$OUT/cache/homo_sapiens/115_GRCh38" ]]; then
  curl -L --fail --retry 5 --retry-delay 5 -C - -o "$OUT/cache/homo_sapiens_vep_115_GRCh38.tar.gz.part" \
    "https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz" >>"$LOG_FILE" 2>&1
  mv "$OUT/cache/homo_sapiens_vep_115_GRCh38.tar.gz.part" "$OUT/cache/homo_sapiens_vep_115_GRCh38.tar.gz"
  tar -xzf "$OUT/cache/homo_sapiens_vep_115_GRCh38.tar.gz" -C "$OUT/cache" >>"$LOG_FILE" 2>&1
fi

if [[ ! -s "$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" ]]; then
  curl -L --fail --retry 5 --retry-delay 5 -C - -o "$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.part" \
    "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" >>"$LOG_FILE" 2>&1
  mv "$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.part" "$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fi

if [[ ! -s "$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]]; then
  gunzip -c "$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" > "$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fi

if [[ ! -s "$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai" ]]; then
  samtools faidx "$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa" >>"$LOG_FILE" 2>&1
fi

cat > "$OUT/toy_input.vcf" <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	881918	.	G	A	.	PASS	.
VCF

singularity exec "$VEP_SIF" vep \
  --offline \
  --cache \
  --dir_cache "$OUT/cache" \
  --species homo_sapiens \
  --assembly GRCh38 \
  --fasta "$OUT/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
  --input_file "$OUT/toy_input.vcf" \
  --output_file "$OUT/toy_output.vep.txt" \
  --tab >>"$LOG_FILE" 2>&1
