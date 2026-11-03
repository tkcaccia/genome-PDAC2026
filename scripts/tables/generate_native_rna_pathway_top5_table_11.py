#!/usr/bin/env python3
# Current collaborator code export. Formatting of derived private tables; see the table index for source gaps.
# Source: Original Source and Code/scripts/generate_native_rna_pathway_top5_table_11.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Condense the native RNA DE pathway-enrichment table to a top-5 summary table."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path

import pandas as pd

ROOT = project_root()
OUT = ROOT / "results" / "rnaseq_revised_11_native_deseq2_20260917"
FULL = OUT / "native_rna_de_pathway_enrichment_11_pdac.csv"


def main() -> None:
    df = pd.read_csv(FULL)
    df = df.sort_values("Adjusted P-value").drop_duplicates(subset=["Term"]).head(5)
    top5 = pd.DataFrame({
        "Pathway name": df["Term"].str.replace(r"\s*\(GO:\d+\)", "", regex=True),
        "Gene number": df["Genes"].str.split(";").apply(len),
        "Gene name/s": df["Genes"].str.replace(";", ", "),
        "pValue": df["P-value"].map(lambda value: f"{value:.2e}"),
        "FDR": df["Adjusted P-value"].map(lambda value: f"{value:.2e}"),
    })
    top5.to_csv(OUT / "native_rna_de_pathway_top5_11_pdac.csv", index=False)
    print(top5.to_string(index=False))


if __name__ == "__main__":
    main()
