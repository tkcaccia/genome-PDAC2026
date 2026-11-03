#!/usr/bin/env python3
# PDAC2026 analysis: Salmon gene counts/TPM, paired DESeq2, response PCA and enrichment.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Build the tumour-vs-normal RNA-seq delta matrix.

Input:
- results/rnaseq/tpm_matrix.csv

Outputs:
- results/rnaseq/rnaseq_delta_matrix_unfiltered.csv
- results/rnaseq/rnaseq_delta_matrix.csv
"""

from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = project_root()
DEFAULT_TPM = PROJECT_ROOT / "results" / "rnaseq" / "tpm_matrix.csv"
DEFAULT_OUTDIR = PROJECT_ROOT / "results" / "rnaseq"
PATIENTS = configured_cohort("all")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tpm", type=Path, default=DEFAULT_TPM)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if not args.tpm.exists():
        raise FileNotFoundError(f"Missing TPM matrix: {args.tpm}")

    tpm = pd.read_csv(args.tpm, index_col=0)
    missing = []
    delta_columns: list[pd.Series] = []

    for patient in PATIENTS:
        normal = f"{patient}N"
        tumour = f"{patient}T"
        if normal not in tpm.columns or tumour not in tpm.columns:
            missing.append(patient)
            continue

        delta = pd.Series(
            np.log2((tpm[tumour] + 1) / (tpm[normal] + 1)),
            index=tpm.index,
            name=f"{int(patient):02d}T",
        )
        delta_columns.append(delta)

    if missing:
        print(f"Warning: missing RNA pairs for patients {', '.join(missing)}")

    delta_matrix = pd.concat(delta_columns, axis=1)
    delta_matrix.index.name = "gene"

    expressed = (tpm > 1).sum(axis=1) >= round(0.2 * tpm.shape[1])
    delta_filtered = delta_matrix.loc[expressed]

    delta_matrix.to_csv(args.outdir / "rnaseq_delta_matrix_unfiltered.csv")
    delta_filtered.to_csv(args.outdir / "rnaseq_delta_matrix.csv")

    print(f"Saved unfiltered delta: {args.outdir / 'rnaseq_delta_matrix_unfiltered.csv'}")
    print(f"Saved filtered delta:   {args.outdir / 'rnaseq_delta_matrix.csv'}")
    print(f"Delta shape: {delta_matrix.shape}")
    print(f"Filtered shape: {delta_filtered.shape}")


if __name__ == "__main__":
    main()
