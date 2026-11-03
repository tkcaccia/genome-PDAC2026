#!/usr/bin/env python3
# Current collaborator code export. Salmon gene counts/TPM, paired DESeq2, response PCA and enrichment.
# Source: Original Source and Code/scripts/02_rnaseq_build_tpm_matrix.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Build RNA-seq TPM and count matrices from Salmon quant files.

Inputs:
- rnaseq/quant_files/*.quant.genes.sf

Outputs:
- results/rnaseq/tpm_matrix.csv
- results/rnaseq/counts_matrix.csv
"""

from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

import argparse
import re
from pathlib import Path

import pandas as pd


PROJECT_ROOT = project_root()
DEFAULT_INPUT_DIR = PROJECT_ROOT / "rnaseq" / "quant_files"
DEFAULT_OUTDIR = PROJECT_ROOT / "results" / "rnaseq"
SAMPLE_PATTERN = re.compile(r"(?P<sample>\d+[NT])RNA_.*\.quant\.genes\.sf$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def sample_id_from_path(path: Path) -> str:
    match = SAMPLE_PATTERN.search(path.name)
    if not match:
        raise ValueError(f"Could not derive sample id from {path.name}")
    return match.group("sample")


def sorted_samples(sample_names: list[str]) -> list[str]:
    return sorted(sample_names, key=lambda sample: (int(sample[:-1]), sample[-1]))


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    quant_files = sorted(args.input_dir.glob("*.quant.genes.sf"))
    if not quant_files:
        raise FileNotFoundError(f"No quant files found in {args.input_dir}")

    sample_files = {sample_id_from_path(path): path for path in quant_files}
    sample_order = sorted_samples(list(sample_files))

    tpm_columns: list[pd.Series] = []
    count_columns: list[pd.Series] = []

    for sample in sample_order:
        path = sample_files[sample]
        table = pd.read_csv(path, sep="\t")
        genes = table["Name"].astype(str).str.replace(r"\..*$", "", regex=True)
        tpm = pd.Series(table["TPM"].values, index=genes, name=sample)
        counts = pd.Series(table["NumReads"].values, index=genes, name=sample)
        tpm_columns.append(tpm)
        count_columns.append(counts)

    tpm_matrix = pd.concat(tpm_columns, axis=1)
    counts_matrix = pd.concat(count_columns, axis=1)

    tpm_matrix.index.name = "gene"
    counts_matrix.index.name = "gene"

    tpm_matrix.to_csv(args.outdir / "tpm_matrix.csv")
    counts_matrix.to_csv(args.outdir / "counts_matrix.csv")

    print(f"Saved TPM matrix:   {args.outdir / 'tpm_matrix.csv'}")
    print(f"Saved counts matrix: {args.outdir / 'counts_matrix.csv'}")
    print(f"TPM shape: {tpm_matrix.shape}")
    print(f"Counts shape: {counts_matrix.shape}")


if __name__ == "__main__":
    main()
