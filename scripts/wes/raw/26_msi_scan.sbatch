#!/usr/bin/env bash
# Publication WES workflow; inputs and outputs stay outside this repository.
#SBATCH --job-name=msi25_scan
#SBATCH --output=logs/msi25_scan_%j.log
#SBATCH --partition=Main
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=06:00:00
set -euo pipefail
source /etc/profile.d/lmod.sh
module use /software/modules/common /software/modules/bio
module load htslib/1.10.2
ROOT="${SLURM_SUBMIT_DIR}"
TOOL="${ROOT}/inputs/msisensor-pro-1.3.0/binary/msisensor-pro-v1.3.0"
REF="${ROOT}/inputs/reference/GCF_000001405.40_GRCh38.p14_genomic.fna"
SITES="${ROOT}/inputs/intervals/msisensor_pro_1.3.0_GRCh38p14.list"
[[ ! -e "${SITES}" ]] || { echo "Catalogue already exists; refusing overwrite" >&2; exit 1; }
trap 'rm -f "${SITES}.tmp"' EXIT
"${TOOL}" scan -d "${REF}" -o "${SITES}.tmp" -l 8 -c 5 -m 50 -s 6 -r 5 -p 0
[[ -s "${SITES}.tmp" ]]
mv "${SITES}.tmp" "${SITES}"
sha256sum "${TOOL}" "${REF}" "${SITES}" > "${ROOT}/manifests/msi_sha256.txt"
git -C "${ROOT}/inputs/msisensor-pro-1.3.0" rev-parse HEAD > "${ROOT}/manifests/msi_source_commit.txt"
wc -l "${SITES}"