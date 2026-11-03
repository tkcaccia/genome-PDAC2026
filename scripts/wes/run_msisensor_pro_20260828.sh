#!/usr/bin/env bash
# PDAC2026 analysis: VEP-annotated DeepSomatic calls, retained-variant summaries and paired MSI launcher.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
#SBATCH --job-name=msi_pro_14
#SBATCH --partition=Main
#SBATCH --array=1-14
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs_msi_pro_20260828/msi_%A_%a.out

set -euo pipefail

WORKDIR="${PDAC_MSI_WORKDIR:?Set a private MSI working directory}"
cd "${WORKDIR}"
module load htslib/1.10.2

PATIENT=$(tail -n +2 sample_manifest.tsv | cut -f2 | sort -n | uniq | sed -n "${SLURM_ARRAY_TASK_ID}p")
TUMOR_BAM=$(grep -P "\t${PATIENT}\tT\t" sample_manifest.tsv | cut -f1 | sed 's/$/.dedup.bam/')
NORMAL_BAM=$(grep -P "\t${PATIENT}\tN\t" sample_manifest.tsv | cut -f1 | sed 's/$/.dedup.bam/')
SITES=tools/microsatellites.list
OUT="results/msi_pro_20260828/patient_${PATIENT}"

echo "START patient=${PATIENT} tumor=${TUMOR_BAM} normal=${NORMAL_BAM} host=$(hostname)"

for file in "${TUMOR_BAM}" "${NORMAL_BAM}" "${TUMOR_BAM}.bai" "${NORMAL_BAM}.bai" "${SITES}"; do
    if [[ ! -f "${file}" ]]; then
        echo "Missing required file: ${file}" >&2
        exit 2
    fi
done

./tools/msisensor-pro msi \
    -d "${SITES}" \
    -n "${NORMAL_BAM}" \
    -t "${TUMOR_BAM}" \
    -o "${OUT}" \
    -c 20 \
    -b 4

echo "END patient=${PATIENT} status=0"
