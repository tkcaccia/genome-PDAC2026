#!/usr/bin/env python3
# PDAC2026 analysis: Formatting of derived private tables; see the table index for source gaps.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Generate revised 11-patient KRAS and MSI supplementary tables as CSV and PNG."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = project_root()
OUT = ROOT / "results" / "figures" / "final_11_figures"
KRAS = ROOT / "results" / "wes_revised_11_20260916" / "outputs" / "primary_pdac" / "wes_somatic_variants.blacklist_filtered.tsv"
MSI = ROOT / "results" / "wes_revised_11_20260916" / "tables" / "msisensor_pro_revised_cohorts.csv"

PROTEIN_CHANGES = {
    (25245350, "C", "A"): "p.G12V",
    (25227342, "T", "A"): "p.Q61L",
    (25245351, "C", "G"): "p.G12R",
    (25245350, "C", "T"): "p.G12D",
}


def render_table(data: pd.DataFrame, columns: list[str], output: Path, title: str, caption: str, widths: list[float]) -> None:
    display = data[columns].astype(str).copy()
    for column in display.columns:
        if column in {"INTERPRETATION", "Protein change"}:
            display[column] = display[column].map(lambda value: textwrap.fill(value, 42))
    line_counts = [max(str(cell).count("\n") + 1 for cell in row) for row in display.itertuples(index=False)]
    fig, ax = plt.subplots(figsize=(16, 1.0 + 0.48 * len(display)))
    ax.axis("off")
    table = ax.table(cellText=display.values, colLabels=columns, cellLoc="left", colLoc="left", loc="upper center")
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1, 1.8)
    for (row, col), cell in table.get_celld().items():
        cell.set_width(widths[col])
        cell.PAD = 0.02
        cell.set_edgecolor("none")
        if row == 0:
            cell.set_text_props(weight="bold")
            cell.visible_edges = "TB"
            cell.set_linewidth(1.4)
            cell.set_edgecolor("black")
        else:
            cell.set_facecolor("#f2f2f2" if row % 2 else "white")
    for col in range(len(columns)):
        table[(1, col)].visible_edges = "T"
        table[(1, col)].set_edgecolor("black")
        table[(1, col)].set_linewidth(1.0)
        table[(len(display), col)].visible_edges = "B"
        table[(len(display), col)].set_edgecolor("black")
        table[(len(display), col)].set_linewidth(1.4)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    variants = pd.read_csv(KRAS, sep="\t")
    kras = variants[variants["SYMBOL"].eq("KRAS")].copy().sort_values("patient")
    kras["PATIENT"] = kras["patient"].astype(str)
    kras["GENOMIC VARIANT"] = kras.apply(lambda row: f"chr12:{int(row['POS'])} {row['REF']}>{row['ALT']}", axis=1)
    kras["PROTEIN CHANGE"] = kras.apply(lambda row: PROTEIN_CHANGES.get((int(row["POS"]), row["REF"], row["ALT"]), "Not annotated"), axis=1)
    kras["DEPTH"] = kras["DP"].astype(int).map(str)
    kras["ALLELE DEPTH"] = kras["AD"].astype(str)
    kras["VAF"] = pd.to_numeric(kras["VAF"], errors="coerce").mul(100).map(lambda value: f"{value:.1f}%")
    kras_table = kras[["PATIENT", "GENOMIC VARIANT", "PROTEIN CHANGE", "DEPTH", "ALLELE DEPTH", "VAF"]]
    kras_table.to_csv(OUT / "Table_1_KRAS_hotspot_variants_11_pdac.csv", index=False)
    render_table(
        kras_table,
        list(kras_table.columns),
        OUT / "Table_1_KRAS_hotspot_variants_11_pdac.png",
        "Table 1. KRAS hotspot variants identified in the revised 11-patient PDAC WES call set",
        "Table 1. KRAS hotspot variants identified in the revised 11-patient PDAC WES call set. Tumour-only somatic calls are shown with tumour depth, alternate/ref allele depth and variant allele fraction.",
        [0.10, 0.25, 0.18, 0.10, 0.16, 0.10],
    )

    msi = pd.read_csv(MSI)
    pdac = msi[msi["cohort"].eq("primary_pdac")].copy()
    msi_rows = []
    for _, row in pdac.sort_values("msi_score_pct", ascending=False).iterrows():
        interpretation = {"MSS": "Microsatellite stable range", "MSI-Intermediate": "Elevated instability; exploratory intermediate case", "MSI-High": "MSI-high candidate; requires orthogonal validation"}.get(row["category"], "")
        msi_rows.append({"PATIENT": str(row["patient"]), "INFORMATIVE SITES": f"{row['informative_sites']:,}", "UNSTABLE SITES": f"{row['unstable_sites']:,}", "MSI SCORE": f"{row['msi_score_pct']:.2f}%", "MSI CATEGORY": row["category"], "INTERPRETATION": interpretation})
    msi_table = pd.DataFrame(msi_rows)
    msi_table.to_csv(OUT / "Supplementary_Table_MSI_scores_11_pdac.csv", index=False)
    render_table(msi_table, list(msi_table.columns), OUT / "Supplementary_Table_MSI_scores_11_pdac.png", "Table S1. Microsatellite instability scores across the revised 11-patient PDAC cohort", "Table S1. Patient-level MSIsensor-pro outputs for the revised 11-patient PDAC cohort. Scores below 5% were treated as MSS, scores from 5% to below 20% as MSI-intermediate, and scores at or above 20% as MSI-high for exploratory interpretation.", [0.13, 0.15, 0.13, 0.11, 0.15, 0.33])
    print("Generated KRAS and MSI CSV/PNG tables")


if __name__ == "__main__":
    main()
