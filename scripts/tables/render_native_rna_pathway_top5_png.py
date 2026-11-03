#!/usr/bin/env python3
# PDAC2026 analysis: Formatting of derived private tables; see the table index for source gaps.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Render the native RNA DE pathway top-5 table as a PNG artefact."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = project_root()
OUT = ROOT / "results" / "rnaseq_revised_11_native_deseq2_20260917"
CSV = OUT / "native_rna_de_pathway_top5_11_pdac.csv"


def wrap(text: str, width: int) -> str:
    return "\n".join(textwrap.wrap(text, width=width)) or text


def main() -> None:
    df = pd.read_csv(CSV, dtype=str)
    df["Gene name/s"] = df["Gene name/s"].apply(lambda genes: wrap(genes, 78))
    df["Pathway name"] = df["Pathway name"].apply(lambda name: wrap(name, 28))
    row_units = [
        max(cell.count("\n") + 1 for cell in row)
        for row in df.itertuples(index=False)
    ]
    header_h = 0.09
    unit_h = 0.85 / sum(row_units)
    row_heights = [unit * unit_h for unit in row_units]

    fig, ax = plt.subplots(figsize=(15, 1.1 + 0.62 * sum(row_units)))
    ax.axis("off")
    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellLoc="left",
        colLoc="left",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    col_widths = [0.17, 0.08, 0.58, 0.085, 0.085]
    for (row, col), cell in table.get_celld().items():
        cell.set_width(col_widths[col])
        cell.set_height(header_h if row == 0 else row_heights[row - 1])
        cell.PAD = 0.01
        cell.set_edgecolor("#999999")
        cell.set_text_props(va="center")
        if row == 0:
            cell.set_text_props(weight="bold", va="center")
            cell.set_facecolor("#c6e0c1")
        else:
            cell.set_facecolor("#ffffff" if row % 2 else "#f4f7f3")
    ax.set_title("Top 5 RNA-seq DE pathway enrichment: 11 PDAC patients (native paired DESeq2)", fontsize=11, pad=10)
    fig.tight_layout()
    fig.savefig(OUT / "native_rna_de_pathway_top5_11_pdac.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
