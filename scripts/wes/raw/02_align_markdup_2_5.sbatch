#!/usr/bin/env bash
# Publication WES workflow; inputs and outputs stay outside this repository.
#SBATCH --job-name=align25
#SBATCH --output=logs/align25_%A_%a.log
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=08:00:00
#SBATCH --array=1-28

set -euo pipefail

source /etc/profile.d/lmod.sh
module use /software/modules/common
module use /software/modules/bio

WORKDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${WORKDIR}"
MANIFEST="${WORKDIR}/manifests/sample_manifest.tsv"
FASTQ_DIR="${WORKDIR}/inputs/fastq"
REFERENCE="${WORKDIR}/inputs/reference/GCF_000001405.40_GRCh38.p14_genomic.fna"
ALIGN_DIR="${WORKDIR}/alignments"
PICARD=(java -jar /software/bio/picard/2.26.11/picard.jar)

module load bwa/github samtools/1.22.1 picard/2.26.11
mkdir -p "${ALIGN_DIR}"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MANIFEST}")
if [[ -z "${LINE}" ]]; then
    echo "No manifest row for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi
IFS=$'\t' read -r SAMPLE PATIENT TYPE ZIP_NAME <<< "${LINE}"
ZIP_PATH="${FASTQ_DIR}/${ZIP_NAME}"
FINAL_BAM="${ALIGN_DIR}/${SAMPLE}.dedup.bam"

if [[ -s "${FINAL_BAM}" && -s "${FINAL_BAM}.bai" && -s "${ALIGN_DIR}/${SAMPLE}.dup_metrics.txt" ]]; then
    echo "${SAMPLE}: deduplicated BAM and metrics already exist; skipping."
    exit 0
fi
[[ -s "${ZIP_PATH}" ]] || { echo "Missing archive: ${ZIP_PATH}" >&2; exit 1; }
[[ -s "${REFERENCE}.bwt" && -s "${REFERENCE}.fai" ]] || { echo "Reference indexes missing for ${REFERENCE}" >&2; exit 1; }

TMP_DIR=$(mktemp -d "${ALIGN_DIR}/.${SAMPLE}.XXXXXX")
trap 'rm -rf "${TMP_DIR}"' EXIT
mkdir -p "${TMP_DIR}/fastq"

echo "Validating archive ${ZIP_NAME}"
unzip -tq "${ZIP_PATH}"
unzip -j -q "${ZIP_PATH}" '*.fastq.gz' -d "${TMP_DIR}/fastq"
R1=$(find "${TMP_DIR}/fastq" -type f -name '*_R1_*.fastq.gz' ! -name '._*' -print -quit)
R2=$(find "${TMP_DIR}/fastq" -type f -name '*_R2_*.fastq.gz' ! -name '._*' -print -quit)
[[ -n "${R1}" && -n "${R2}" ]] || { echo "Could not find R1/R2 in ${ZIP_NAME}" >&2; exit 1; }
gzip -t "${R1}"
gzip -t "${R2}"

SORTED_BAM="${TMP_DIR}/${SAMPLE}.sorted.bam"
DEDUP_BAM="${TMP_DIR}/${SAMPLE}.dedup.bam"
METRICS="${TMP_DIR}/${SAMPLE}.dup_metrics.txt"

bwa mem -t 8 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
    "${REFERENCE}" "${R1}" "${R2}" \
    | samtools sort -@ 8 -T "${TMP_DIR}/sort" -o "${SORTED_BAM}" -
samtools index "${SORTED_BAM}"
"${PICARD[@]}" MarkDuplicates I="${SORTED_BAM}" O="${DEDUP_BAM}" M="${METRICS}"
samtools index "${DEDUP_BAM}"

mv "${DEDUP_BAM}" "${FINAL_BAM}"
mv "${DEDUP_BAM}.bai" "${FINAL_BAM}.bai"
mv "${METRICS}" "${ALIGN_DIR}/${SAMPLE}.dup_metrics.txt"
echo "${SAMPLE}: alignment and duplicate marking complete."