#!/usr/bin/env python3
# Current collaborator code export. Protein-group matrix assembly and integration/QC normalization; not the missing proteomics DE analysis.
# Source: Original Source and Code/scripts/04_proteomics_build_protein_matrix.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Build a sample-level proteomics protein-group abundance matrix.

Inputs:
- proteomics/20260116_123225_20241023_FFPE TISSUE_PDAC_SNE_Report.tsv

Outputs:
- results/proteomics/protein_abundance_matrix.csv
- results/proteomics/protein_sample_metadata.csv
"""

from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

import argparse
import re
from pathlib import Path

import pandas as pd


PROJECT_ROOT = project_root()
DEFAULT_INPUT = PROJECT_ROOT / "proteomics" / "20260116_123225_20241023_FFPE TISSUE_PDAC_SNE_Report.tsv"
DEFAULT_OUTDIR = PROJECT_ROOT / "results" / "proteomics"
LABEL_PATTERN = re.compile(r"_(?P<patient>\d+)(?P<condition>[NTP])\.d$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def sample_id_from_label(label: str) -> str | None:
    match = LABEL_PATTERN.search(label)
    if not match:
        return None
    return f"{match.group('patient')}{match.group('condition')}"


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if not args.input.exists():
        raise FileNotFoundError(f"Missing proteomics report: {args.input}")

    usecols = ["R.Condition", "R.Label", "PG.ProteinAccessions", "PG.Quantity", "EG.IsDecoy"]
    table = pd.read_csv(args.input, sep="\t", usecols=usecols, low_memory=False)

    table = table.copy()
    table["sample_id"] = table["R.Label"].astype(str).map(sample_id_from_label)
    table = table[table["sample_id"].notna()]
    table = table[table["R.Condition"].isin(["Normal", "Tumour"])]

    if "EG.IsDecoy" in table.columns:
        table = table[~table["EG.IsDecoy"].astype(str).str.lower().isin(["true", "1"])]

    table["PG.Quantity"] = (
        table["PG.Quantity"]
        .astype(str)
        .str.replace(",", ".", regex=False)
        .pipe(pd.to_numeric, errors="coerce")
    )
    table = table.dropna(subset=["PG.Quantity"])

    metadata = (
        table[["sample_id", "R.Condition", "R.Label"]]
        .drop_duplicates()
        .sort_values(["sample_id", "R.Condition", "R.Label"])
    )

    protein_matrix = (
        table.pivot_table(
            index="PG.ProteinAccessions",
            columns="sample_id",
            values="PG.Quantity",
            aggfunc="median",
        )
        .sort_index(axis=1)
    )

    metadata.to_csv(args.outdir / "protein_sample_metadata.csv", index=False)
    protein_matrix.to_csv(args.outdir / "protein_abundance_matrix.csv")

    print(f"Saved metadata:  {args.outdir / 'protein_sample_metadata.csv'}")
    print(f"Saved matrix:    {args.outdir / 'protein_abundance_matrix.csv'}")
    print(f"Matrix shape: {protein_matrix.shape}")
    print(f"Samples: {protein_matrix.shape[1]}")


if __name__ == "__main__":
    main()
