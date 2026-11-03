#!/usr/bin/env python3
# PDAC2026 analysis: Separate exploratory PDAC/periampullary comparisons.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Exploratory PDAC versus periampullary comparison across WES, RNA, and proteomics."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = project_root()
OUT = ROOT / "results" / "pdac_ampulla_comparison_20260917"
PDAC = configured_cohort("primary_pdac")
AMPULLA = configured_cohort("periampullary")
ALL = PDAC + AMPULLA


def pca_plot(matrix: pd.DataFrame, name: str, x_label: str) -> None:
    values = matrix.loc[:, ALL].T
    transformed = PCA(n_components=2, random_state=42).fit_transform(StandardScaler().fit_transform(values))
    coords = pd.DataFrame(transformed, index=values.index, columns=["PC1", "PC2"])
    coords["group"] = ["PDAC" if patient in PDAC else "Ampulla" for patient in coords.index]
    coords.to_csv(OUT / "tables" / f"{name}_pca_coordinates.csv")
    fig, ax = plt.subplots(figsize=(7, 6))
    for group, color in [("PDAC", "#457b9d"), ("Ampulla", "#e76f51")]:
        subset = coords[coords.group == group]
        ax.scatter(subset.PC1, subset.PC2, s=90, color=color, edgecolor="black", label=group)
        for patient, row in subset.iterrows():
            ax.annotate(patient, (row.PC1, row.PC2), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set(xlabel="PC1", ylabel="PC2", title=f"{x_label}: PDAC versus ampulla")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "figures" / f"{name}_pca.png", dpi=300)
    plt.close(fig)


def group_tests(matrix: pd.DataFrame, name: str, symbol_column: str | None = None) -> None:
    rows = []
    pdac = matrix.loc[:, PDAC]
    ampulla = matrix.loc[:, AMPULLA]
    for feature in matrix.index:
        a = pd.to_numeric(pdac.loc[feature], errors="coerce").dropna()
        b = pd.to_numeric(ampulla.loc[feature], errors="coerce").dropna()
        if len(a) < 3 or len(b) < 2:
            continue
        statistic, pvalue = mannwhitneyu(a, b, alternative="two-sided", method="exact")
        rows.append({"feature": feature, "pdac_mean": a.mean(), "ampulla_mean": b.mean(), "mean_difference_ampulla_minus_pdac": b.mean() - a.mean(), "p_value": pvalue})
    result = pd.DataFrame(rows).sort_values("p_value")
    ranks = result["p_value"].rank(method="first")
    result["FDR"] = (result["p_value"] * len(result) / ranks).cummin().clip(upper=1)
    result.to_csv(OUT / "tables" / f"{name}_group_comparison.csv", index=False)


def wes_analysis() -> None:
    path = ROOT / "results" / "wes/wes_mutation_matrix.blacklist_filtered.csv"
    matrix = pd.read_csv(path, index_col=0)
    matrix.columns = [str(col).lstrip("0") for col in matrix.columns]
    matrix = matrix.loc[:, [patient + "T" for patient in ALL]]
    burden = matrix.sum(axis=0).rename("mutation_count").to_frame()
    burden["group"] = ["PDAC" if patient[:-1] in PDAC else "Ampulla" for patient in burden.index]
    burden.to_csv(OUT / "tables/wes_mutation_burden_by_patient.csv")
    fig, ax = plt.subplots(figsize=(7, 5))
    for group, color in [("PDAC", "#457b9d"), ("Ampulla", "#e76f51")]:
        values = burden[burden.group == group].mutation_count
        ax.scatter([group] * len(values), values, color=color, s=85, edgecolor="black", label=group)
    ax.set(ylabel="Blacklist-filtered mutated genes", title="WES mutation burden: PDAC versus ampulla")
    fig.tight_layout()
    fig.savefig(OUT / "figures/wes_mutation_burden_comparison.png", dpi=300)
    plt.close(fig)
    freq = pd.DataFrame({"PDAC_frequency": matrix.loc[:, [patient + "T" for patient in PDAC]].sum(axis=1) / len(PDAC), "Ampulla_frequency": matrix.loc[:, [patient + "T" for patient in AMPULLA]].sum(axis=1) / len(AMPULLA)})
    freq["difference_ampulla_minus_pdac"] = freq.Ampulla_frequency - freq.PDAC_frequency
    freq.sort_values(["Ampulla_frequency", "PDAC_frequency"], ascending=False).to_csv(OUT / "tables/wes_gene_frequency_comparison.csv")


def load_rna() -> pd.DataFrame:
    path = ROOT / "results/rnaseq/rnaseq_delta_matrix.csv"
    matrix = pd.read_csv(path, index_col=0)
    matrix.columns = [str(col).lstrip("0").rstrip("T") for col in matrix.columns]
    return matrix.loc[:, ALL]


def load_proteomics() -> pd.DataFrame:
    path = ROOT / "results/proteomics/protein_abundance_matrix.csv"
    matrix = pd.read_csv(path, index_col=0)
    paired = {}
    for patient in ALL:
        normal = np.log2(pd.to_numeric(matrix[f"{patient}N"], errors="coerce") + 1)
        tumour = np.log2(pd.to_numeric(matrix[f"{patient}T"], errors="coerce") + 1)
        paired[patient] = tumour - normal
    return pd.DataFrame(paired, index=matrix.index).replace([np.inf, -np.inf], np.nan).fillna(0)


def main() -> None:
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    wes_analysis()
    rna = load_rna()
    pca_plot(rna, "rna_response", "RNA tumour-normal log2 response")
    group_tests(rna, "rna_response")
    proteomics = load_proteomics()
    pca_plot(proteomics, "proteomics_response", "Proteomics tumour-normal log2 response")
    group_tests(proteomics, "proteomics_response")
    pd.DataFrame({"group": ["PDAC", "Ampulla"], "n": [len(PDAC), len(AMPULLA)], "patients": [",".join(PDAC), ",".join(AMPULLA)]}).to_csv(OUT / "tables/cohort_summary.csv", index=False)
    print(f"Generated comparison outputs under {OUT}")


if __name__ == "__main__":
    main()
