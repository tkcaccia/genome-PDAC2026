#!/usr/bin/env bash
set -euo pipefail

# Reproduce the enhanced per-patient circos plots on the remote PDAC workstation.
# Result PNGs are patient-derived outputs and should not be committed to GitHub.

MUTECT2_ROOT="${MUTECT2_ROOT:?Set MUTECT2_ROOT}"
VEP_ROOT="${VEP_ROOT:?Set VEP_ROOT}"
ASCAT_ROOT="${ASCAT_ROOT:?Set ASCAT_ROOT}"
MANTA_ROOT="${MANTA_ROOT:?Set MANTA_ROOT}"
STAR_ROOT="${STAR_ROOT:?Set STAR_ROOT}"
GTF="${GTF:?Set GTF to the matching gene annotation}"
FAI="${FAI:?Set FAI to the reference FASTA index}"
CYTOBAND="${CYTOBAND:?Set CYTOBAND to the hg38 cytoband path}"
OUTDIR="${OUTDIR:?Set OUTDIR outside this repository}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$SCRIPT_DIR/make_patient_circos.py"

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
