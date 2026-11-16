#!/usr/bin/env bash
# Publication WES workflow; inputs and outputs stay outside this repository.
#SBATCH --job-name=fastqc25
#SBATCH --output=logs/fastqc25_%A_%a.log
#SBATCH --partition=Main
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=02:00:00
#SBATCH --array=1-28

set -euo pipefail

source /etc/profile.d/lmod.sh
module use /software/modules/common
module use /software/modules/bio
module load fastqc/0.11.9

WORKDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
MANIFEST="${WORKDIR}/manifests/sample_manifest.tsv"
FASTQ_DIR="${WORKDIR}/inputs/fastq"
OUT_DIR="${WORKDIR}/qc/fastqc_2_5"
mkdir -p "${OUT_DIR}"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MANIFEST}")
if [[ -z "${LINE}" ]]; then
    echo "No manifest row for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi
IFS=$'\t' read -r SAMPLE PATIENT TYPE ZIP_NAME <<< "${LINE}"
ZIP_PATH="${FASTQ_DIR}/${ZIP_NAME}"
TMP_DIR=$(mktemp -d "${OUT_DIR}/.${SAMPLE}.XXXXXX")
trap 'rm -rf "${TMP_DIR}"' EXIT

[[ -s "${ZIP_PATH}" ]] || { echo "Missing ZIP archive: ${ZIP_PATH}" >&2; exit 1; }
unzip -tq "${ZIP_PATH}"
unzip -j -q "${ZIP_PATH}" '*.fastq.gz' -d "${TMP_DIR}"
R1=$(find "${TMP_DIR}" -type f -name '*_R1_*.fastq.gz' ! -name '._*' -print -quit)
R2=$(find "${TMP_DIR}" -type f -name '*_R2_*.fastq.gz' ! -name '._*' -print -quit)
[[ -n "${R1}" && -n "${R2}" ]] || { echo "Could not find R1/R2 in ${ZIP_NAME}" >&2; exit 1; }
gzip -t "${R1}"
gzip -t "${R2}"
fastqc -t 4 --outdir "${TMP_DIR}" "${R1}" "${R2}"

for source in "${TMP_DIR}"/*_fastqc.html "${TMP_DIR}"/*_fastqc.zip; do
    [[ -f "${source}" ]] || continue
    base=$(basename "${source}")
    mv "${source}" "${OUT_DIR}/${SAMPLE}_${base}"
done
echo "${SAMPLE}: FastQC complete."