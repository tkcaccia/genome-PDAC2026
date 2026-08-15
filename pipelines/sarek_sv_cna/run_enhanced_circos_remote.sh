#!/usr/bin/env bash
set -euo pipefail

# Reproduce the enhanced per-patient circos plots on the remote PDAC workstation.
# Result PNGs are patient-derived outputs and should not be committed to GitHub.

MUTECT2_ROOT="/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal/variant_calling/mutect2"
VEP_ROOT="/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal/annotation/mutect2"
ASCAT_ROOT="/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal_sv_cna/variant_calling/ascat"
MANTA_ROOT="/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal_sv_cna/variant_calling/manta"
STAR_ROOT="/home/user/PDAC_SEQ_native_results/rnafusion_pdac/star"
GTF="/media/user/SEQ/refs/annotation/gencode.v46.primary_assembly.annotation.gtf"
FAI="/media/user/SEQ/refs/gatk_bundle/Homo_sapiens_assembly38.fasta.fai"
CYTOBAND="/media/user/SEQ/refs/ucsc/hg38.cytoBand.txt"
OUTDIR="/media/user/PDAC_SEQ_analysis/results/sarek_tumor_normal_sv_cna/circos_plots_enhanced"
SCRIPT="/media/user/SEQ/scripts/make_patient_circos.py"

mkdir -p "$(dirname "$CYTOBAND")" "$OUTDIR"

if [[ ! -s "$CYTOBAND" ]]; then
  if command -v curl >/dev/null 2>&1; then
    curl -L --retry 3 \
      -o "${CYTOBAND}.gz" \
      "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
  else
    wget -O "${CYTOBAND}.gz" \
      "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
  fi
  gzip -dc "${CYTOBAND}.gz" > "$CYTOBAND"
fi

python3 "$SCRIPT" \
  --mutect2-root "$MUTECT2_ROOT" \
  --vep-root "$VEP_ROOT" \
  --ascat-root "$ASCAT_ROOT" \
  --manta-root "$MANTA_ROOT" \
  --star-root "$STAR_ROOT" \
  --gtf "$GTF" \
  --cytoband "$CYTOBAND" \
  --fai "$FAI" \
  --output-dir "$OUTDIR"

find "$OUTDIR" -maxdepth 1 -name "*.png" | sort
