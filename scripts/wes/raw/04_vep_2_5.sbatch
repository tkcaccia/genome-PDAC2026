#!/usr/bin/env bash
# Publication WES workflow; inputs and outputs stay outside this repository.
#SBATCH --job-name=vep25
#SBATCH --output=logs/vep25_%A_%a.log
#SBATCH --partition=Main
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=03:00:00
#SBATCH --array=1-14

set -euo pipefail

source /etc/profile.d/lmod.sh
module use /software/modules/common
module use /software/modules/bio

WORKDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${WORKDIR}"
MANIFEST="${WORKDIR}/manifests/sample_manifest.tsv"
REFERENCE="${WORKDIR}/inputs/reference/GCF_000001405.40_GRCh38.p14_genomic.fna"
CALL_DIR="${WORKDIR}/wes/calls"
ANNOTATED_DIR="${WORKDIR}/wes/annotated"
mkdir -p "${ANNOTATED_DIR}"

module load vep/101.0 bcftools/1.22

PATIENTS=($(tail -n +2 "${MANIFEST}" | cut -f2 | sort -n -u))
PATIENT="${PATIENTS[$((SLURM_ARRAY_TASK_ID - 1))]:-}"
[[ -n "${PATIENT}" ]] || { echo "No patient for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
IN_VCF="${CALL_DIR}/patient_${PATIENT}.somatic.vcf.gz"
OUT_VCF="${ANNOTATED_DIR}/patient_${PATIENT}.vep.vcf.gz"
TMP_VCF="${ANNOTATED_DIR}/.patient_${PATIENT}.${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}.vep.vcf.gz"

[[ -s "${IN_VCF}" ]] || { echo "Missing input VCF: ${IN_VCF}" >&2; exit 1; }
[[ -s "${REFERENCE}" ]] || { echo "Missing reference FASTA: ${REFERENCE}" >&2; exit 1; }
[[ -d /software/data/bio/vep/homo_sapiens/101_GRCh38 ]] || { echo "VEP cache 101_GRCh38 unavailable" >&2; exit 1; }
if [[ -s "${OUT_VCF}" && -s "${OUT_VCF}.tbi" ]]; then
    echo "Annotated output already exists: ${OUT_VCF}; skipping."
    exit 0
fi
trap 'rm -f "${TMP_VCF}" "${TMP_VCF}.tbi"' EXIT

vep \
    --input_file "${IN_VCF}" \
    --output_file STDOUT \
    --vcf \
    --compress_output bgzip \
    --fasta "${REFERENCE}" \
    --offline \
    --cache \
    --dir_cache /software/data/bio/vep \
    --dir_plugins /software/data/bio/vep/Plugins \
    --cache_version 101 \
    --assembly GRCh38 \
    --no_stats \
    --symbol \
    --biotype \
    --canonical \
    --sift b \
    --polyphen b \
    --fork 4 \
    --vcf_info_field CSQ \
    | bcftools +split-vep -c SYMBOL,IMPACT,Consequence,SIFT,PolyPhen -p "CSQ_" -Oz -o "${TMP_VCF}" -

bcftools index -t "${TMP_VCF}"
mv "${TMP_VCF}" "${OUT_VCF}"
mv "${TMP_VCF}.tbi" "${OUT_VCF}.tbi"
echo "Patient ${PATIENT}: VEP 101 annotation complete."