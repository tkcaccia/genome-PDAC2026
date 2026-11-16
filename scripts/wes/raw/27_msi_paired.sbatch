#!/usr/bin/env bash
# Publication WES workflow; inputs and outputs stay outside this repository.
#SBATCH --job-name=msi25_paired
#SBATCH --output=logs/msi25_%A_%a.log
#SBATCH --partition=Main
#SBATCH --array=1-14%4
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
set -euo pipefail
source /etc/profile.d/lmod.sh
module use /software/modules/common /software/modules/bio
module load htslib/1.10.2
ROOT="${SLURM_SUBMIT_DIR}"
TOOL="${ROOT}/inputs/msisensor-pro-1.3.0/binary/msisensor-pro-v1.3.0"
SITES="${ROOT}/inputs/intervals/msisensor_pro_1.3.0_GRCh38p14.list"
MANIFEST="${ROOT}/manifests/sample_manifest.tsv"
PATIENT=$(tail -n +2 "${MANIFEST}" | cut -f2 | sort -nu | sed -n "${SLURM_ARRAY_TASK_ID}p")
[[ -n "${PATIENT}" ]]
TUMOR=$(awk -F '\t' -v p="${PATIENT}" '$2==p && $3=="T" {print $1}' "${MANIFEST}")
NORMAL=$(awk -F '\t' -v p="${PATIENT}" '$2==p && $3=="N" {print $1}' "${MANIFEST}")
TUMOR_BAM="${ROOT}/alignments/${TUMOR}.dedup.bam"
NORMAL_BAM="${ROOT}/alignments/${NORMAL}.dedup.bam"
OUTDIR="${ROOT}/wes/msi_2_5"
mkdir -p "${OUTDIR}"
for INPUT in "${TOOL}" "${SITES}" "${TUMOR_BAM}" "${NORMAL_BAM}" "${TUMOR_BAM}.bai" "${NORMAL_BAM}.bai"; do
    [[ -s "${INPUT}" ]] || { echo "Missing input: ${INPUT}" >&2; exit 1; }
done
[[ ! -e "${OUTDIR}/patient_${PATIENT}" ]] || { echo "Existing result; refusing overwrite" >&2; exit 1; }
TMP=$(mktemp -d "${OUTDIR}/.patient_${PATIENT}.XXXXXX")
echo "Patient=${PATIENT}; paired msi; coverage=20; threads=4; software=v1.3.0"
"${TOOL}" msi -d "${SITES}" -n "${NORMAL_BAM}" -t "${TUMOR_BAM}" -o "${TMP}/patient_${PATIENT}" -c 20 -b 4
[[ -s "${TMP}/patient_${PATIENT}" ]]
for OUTPUT in "${TMP}"/patient_"${PATIENT}"*; do mv "${OUTPUT}" "${OUTDIR}/"; done
rmdir "${TMP}"
echo "Patient ${PATIENT} completed successfully"