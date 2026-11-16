#!/usr/bin/env python3
# WES processing, variant summaries, coverage and paired MSI.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""
from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = private_path("handoff_root")
RUN = ROOT / "wes" / "analysis"
OUT = RUN / "figures" / "primary_pdac"
TABLES = RUN / "tables"
PATIENTS = configured_cohort("primary_pdac")
PROTEIN_CHANGES = {
    (25245350, "C", "A"): "p.G12V",
    (25227342, "T", "A"): "p.Q61L",
    (25245351, "C", "G"): "p.G12R",
    (25245350, "C", "T"): "p.G12D",
}


def save_table(data, columns, path):
    display = data[columns].astype(str).copy()
    fig, ax = plt.subplots(figsize=(16, 1.0 + 0.48 * len(display)))
    ax.axis("off")
    table = ax.table(cellText=display.values, colLabels=columns, cellLoc="left", colLoc="left", loc="upper center")
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1, 1.8)
    widths = [0.10, 0.25, 0.18, 0.10, 0.16, 0.10]
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
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    cohort = RUN / "outputs" / "primary_pdac"
    burden = pd.read_csv(cohort / "wes_mutation_burden.csv")
    burden["patient"] = burden["patient"].astype(str)
    burden = burden.set_index("patient").reindex(PATIENTS).reset_index()

    fig, ax = plt.subplots(figsize=(9, 4))
    bars = ax.bar(burden.patient, burden.somatic_variant_count, color="#457b9d", edgecolor="white")
    ax.set_xlabel("Patient")
    ax.set_ylabel("Filtered somatic variants")
    ax.axhline(burden.somatic_variant_count.median(), color="grey", linestyle="--")
    for bar, count in zip(bars, burden.somatic_variant_count):
        ax.text(bar.get_x() + bar.get_width() / 2, count + 2, str(int(count)), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_1A_somatic_mutation_burden.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    matrix = pd.read_csv(cohort / "wes_mutation_matrix.blacklist_filtered.csv", index_col=0)
    matrix.columns = [str(column).lstrip("0") for column in matrix.columns]
    samples = [f"{patient}T" for patient in PATIENTS]
    matrix = matrix.reindex(columns=samples, fill_value=0)
    top = matrix.sum(axis=1).sort_values(ascending=False).head(30).index
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(matrix.loc[top].astype(int), cmap=["#EEEEEE", "#4B3F8C"], cbar=False, linewidths=0.5, linecolor="white", ax=ax)
    ax.set_xticklabels(PATIENTS, rotation=0)
    ax.set_xlabel("Patient")
    ax.set_ylabel("Gene")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_1B_recurrent_mutation_oncoprint.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    variants = pd.read_csv(cohort / "wes_somatic_variants.blacklist_filtered.tsv", sep="\t")
    gene_summary = pd.read_csv(cohort / "gene_summary_table.csv")
    selected_genes = gene_summary.sort_values("variant_count", ascending=False).head(20).sort_values("variant_count")
    gene_counts = selected_genes.set_index("SYMBOL")["variant_count"]
    gene_counts.to_csv(TABLES / "Figure_1C_top20_variant_count.csv")
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh(gene_counts.index, gene_counts.values, color=plt.cm.Blues(np.linspace(0.35, 0.95, len(gene_counts))))
    ax.set_xlabel("Filtered variant count")
    ax.set_ylabel("Gene")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_1C_top20_variant_count.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    kras = variants[variants.SYMBOL.eq("KRAS")].copy().sort_values("patient")
    kras_table = pd.DataFrame({
        "PATIENT": kras.patient.astype(str),
        "GENOMIC VARIANT": kras.apply(lambda row: f"chr12:{int(row.POS)} {row.REF}>{row.ALT}", axis=1),
        "PROTEIN CHANGE": kras.apply(lambda row: PROTEIN_CHANGES.get((int(row.POS), row.REF, row.ALT), "Not annotated"), axis=1),
        "DEPTH": kras.DP.astype(int).astype(str),
        "ALLELE DEPTH": kras.AD.astype(str),
        "VAF": pd.to_numeric(kras.VAF, errors="coerce").mul(100).map(lambda value: f"{value:.1f}%"),
    })
    kras_table.to_csv(TABLES / "Table_1_KRAS_hotspot_variants_11_pdac.csv", index=False)
    save_table(kras_table, list(kras_table.columns), OUT / "Table_1_KRAS_hotspot_variants_11_pdac.png")
    print(f"Primary PDAC panels and KRAS table written under {OUT}")
    print(f"KRAS rows: {len(kras_table)}; patients: {', '.join(kras_table.PATIENT)}")


if __name__ == "__main__":
    main()