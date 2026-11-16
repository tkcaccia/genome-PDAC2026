#!/usr/bin/env bash
# Publication WES workflow; inputs and outputs stay outside this repository.
#SBATCH --job-name=deepsomatic25
#SBATCH --output=logs/deepsomatic25_%A_%a.log
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --array=1-14

set -euo pipefail

source /etc/profile.d/lmod.sh
module use /software/modules/common
module load singularity/4.4.1

WORKDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${WORKDIR}"
MANIFEST="${WORKDIR}/manifests/sample_manifest.tsv"
REFERENCE="${WORKDIR}/inputs/reference/GCF_000001405.40_GRCh38.p14_genomic.fna"
BED="${WORKDIR}/inputs/intervals/hg38_Twist_Bioscience_for_Illumina_Exome_2_5.refseq.bed"
SIF="${DEEPSOMATIC_SIF:-${WORKDIR}/inputs/reference/deepsomatic_1.10.0.sif}"
CALL_DIR="${WORKDIR}/wes/calls"
mkdir -p "${CALL_DIR}"

PATIENTS=($(tail -n +2 "${MANIFEST}" | cut -f2 | sort -n -u))
PATIENT="${PATIENTS[$((SLURM_ARRAY_TASK_ID - 1))]:-}"
[[ -n "${PATIENT}" ]] || { echo "No patient for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
TUMOR_SAMPLE=$(awk -F '\t' -v p="${PATIENT}" '$2==p && $3=="T" {print $1; exit}' "${MANIFEST}")
NORMAL_SAMPLE=$(awk -F '\t' -v p="${PATIENT}" '$2==p && $3=="N" {print $1; exit}' "${MANIFEST}")
TUMOR_BAM="${WORKDIR}/alignments/${TUMOR_SAMPLE}.dedup.bam"
NORMAL_BAM="${WORKDIR}/alignments/${NORMAL_SAMPLE}.dedup.bam"
OUT_VCF="${CALL_DIR}/patient_${PATIENT}.somatic.vcf.gz"

[[ -s "${TUMOR_BAM}" && -s "${NORMAL_BAM}" ]] || { echo "Missing matched BAMs for patient ${PATIENT}" >&2; exit 1; }
[[ -s "${BED}" ]] || { echo "Missing remapped 2.5 BED: ${BED}; run 01_prepare_capture_2_5.sh first" >&2; exit 1; }
[[ -s "${SIF}" ]] || { echo "Missing DeepSomatic image: ${SIF}" >&2; exit 1; }
if [[ -s "${OUT_VCF}" ]]; then
    echo "Output already exists: ${OUT_VCF}; refusing to overwrite."
    exit 1
fi

TMP_VCF="${CALL_DIR}/.patient_${PATIENT}.${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}.vcf.gz"
trap 'rm -f "${TMP_VCF}"' EXIT

singularity exec "${SIF}" run_deepsomatic \
    --model_type=FFPE_WES \
    --ref="${REFERENCE}" \
    --reads_tumor="${TUMOR_BAM}" \
    --reads_normal="${NORMAL_BAM}" \
    --sample_name_tumor="${PATIENT}T" \
    --sample_name_normal="${PATIENT}N" \
    --output_vcf="${TMP_VCF}" \
    --regions="${BED}" \
    --num_shards=8

mv "${TMP_VCF}" "${OUT_VCF}"
echo "Patient ${PATIENT}: DeepSomatic 1.10.0 FFPE_WES complete using the Exome 2.5 BED."