#!/usr/bin/env bash
# Publication WES workflow; inputs and outputs stay outside this repository.
#SBATCH --job-name=bed25_intervals
#SBATCH --output=logs/bed25_intervals_%j.log
#SBATCH --partition=Main
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=00:30:00

set -euo pipefail

source /etc/profile.d/lmod.sh
module use /software/modules/common
module use /software/modules/bio

WORKDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DICT="${WORKDIR}/inputs/reference/GCF_000001405.40_GRCh38.p14_genomic.dict"
BED_IN="${WORKDIR}/inputs/intervals/hg38_Twist_Bioscience_for_Illumina_Exome_2_5.bed"
BED_OUT="${WORKDIR}/inputs/intervals/hg38_Twist_Bioscience_for_Illumina_Exome_2_5.refseq.bed"
INTERVAL_LIST="${WORKDIR}/inputs/intervals/hg38_Twist_Bioscience_for_Illumina_Exome_2_5.refseq.interval_list"

module load picard/2.26.11

python3 "${WORKDIR}/scripts/remap_bed_to_refseq.py" "${DICT}" "${BED_IN}" "${BED_OUT}"

picard BedToIntervalList \
    I="${BED_OUT}" \
    O="${INTERVAL_LIST}" \
    SD="${DICT}" \
    DROP_MISSING_CONTIGS=false

echo "BED and interval list ready: ${BED_OUT}"