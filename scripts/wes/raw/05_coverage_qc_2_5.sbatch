#!/usr/bin/env bash
# Publication WES workflow; inputs and outputs stay outside this repository.
#SBATCH --job-name=coverage25
#SBATCH --output=logs/coverage25_%A_%a.log
#SBATCH --partition=Main
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --array=1-28

set -euo pipefail

source /etc/profile.d/lmod.sh
module use /software/modules/common
module use /software/modules/bio

WORKDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${WORKDIR}"
MANIFEST="${WORKDIR}/manifests/sample_manifest.tsv"
REFERENCE="${WORKDIR}/inputs/reference/GCF_000001405.40_GRCh38.p14_genomic.fna"
INTERVAL_LIST="${WORKDIR}/inputs/intervals/hg38_Twist_Bioscience_for_Illumina_Exome_2_5.refseq.interval_list"
ALIGN_DIR="${WORKDIR}/alignments"
QC_DIR="${WORKDIR}/qc/coverage_2_5"
PICARD=(java -jar /software/bio/picard/2.26.11/picard.jar)
KRAS_CONTIG="NC_000012.12"

module load picard/2.26.11 samtools/1.22.1
mkdir -p "${QC_DIR}"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MANIFEST}")
if [[ -z "${LINE}" ]]; then
    echo "No manifest row for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi
IFS=$'\t' read -r SAMPLE PATIENT TYPE ZIP_NAME <<< "${LINE}"
BAM="${ALIGN_DIR}/${SAMPLE}.dedup.bam"
METRICS="${QC_DIR}/${SAMPLE}.hsmetrics.txt"
HOTSPOT="${QC_DIR}/${SAMPLE}.kras_hotspot_depth_v2.tsv"

[[ -s "${BAM}" && -s "${BAM}.bai" ]] || { echo "Missing BAM/index: ${BAM}" >&2; exit 1; }
[[ -s "${INTERVAL_LIST}" ]] || { echo "Missing 2.5 interval list: ${INTERVAL_LIST}" >&2; exit 1; }

if [[ ! -s "${METRICS}" ]]; then
    "${PICARD[@]}" CollectHsMetrics \
        I="${BAM}" \
        O="${METRICS}" \
        R="${REFERENCE}" \
        BAIT_INTERVALS="${INTERVAL_LIST}" \
        TARGET_INTERVALS="${INTERVAL_LIST}"
fi

if [[ ! -s "${HOTSPOT}" ]]; then
    TMP_HOTSPOT="${HOTSPOT}.tmp.${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}"
    trap 'rm -f "${TMP_HOTSPOT}"' EXIT
    printf 'sample\tpatient\ttissue\tgene\tlabel\tcontig\tposition_1based\tdepth\n' > "${TMP_HOTSPOT}"
    for SITE in "KRAS_codon13:25245347" "KRAS_codon12_G12:25245350" "KRAS_codon12_G12_second_base:25245351" "KRAS_codon61_Q61L:25227342"; do
        LABEL="${SITE%%:*}"
        POSITION="${SITE##*:}"
        samtools depth -aa -r "${KRAS_CONTIG}:${POSITION}-${POSITION}" "${BAM}" \
            | awk -v s="${SAMPLE}" -v p="${PATIENT}" -v t="${TYPE}" -v l="${LABEL}" \
                'BEGIN{OFS="\t"} $2 == '"${POSITION}"' {print s,p,t,"KRAS",l,$1,$2,$3}' \
            >> "${TMP_HOTSPOT}"
    done
    mv "${TMP_HOTSPOT}" "${HOTSPOT}"
fi

echo "${SAMPLE}: target-coverage and KRAS hotspot depth QC complete."