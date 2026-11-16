#!/usr/bin/env bash
# Publication WES workflow; inputs and outputs stay outside this repository.
#SBATCH --job-name=msi25_summary
#SBATCH --output=logs/msi25_summary_%j.log
#SBATCH --partition=Main
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR}"
${PYTHON:-python3} "${ROOT}/scripts/28_summarize_msi.py"