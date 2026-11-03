#!/usr/bin/env python3
# Current collaborator code export. Public-cohort preparation and descriptive comparisons; not ancestry-controlled validation.
# Source: Original Source and Code/scripts/generate_figure13_rna_pca_tcga_11.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Valid joint PCA of 11 PDAC tumour expression profiles and TCGA PAAD RNA expression."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = project_root()
LOCAL = ROOT / "results/rnaseq_revise_11_archive/rnaseq_revised_11_20260916/tables/tpm_matrix_11_pdac.csv"
TCGA = ROOT / "external/tcga_paad/processed/tcga_rna_matrix.csv"
MAPPING = ROOT / "reproducible_pipeline_share_20260713_with_outputs/raw_outputs/results/figures/categories/04_cross_cohort_concordance/rnaseq_ensembl_to_symbol.csv"
OUT = ROOT / "results/figures/final_11_figures/Figure_13_rna_pca_tcga_11_pdac.png"
COORDS = ROOT / "results/figures/final_11_figures/Figure_13_rna_pca_tcga_11_pdac_coordinates.csv"


def main() -> None:
    local = pd.read_csv(LOCAL, index_col=0)
    local = local[[f"{patient}T" for patient in configured_cohort("primary_pdac")]]
    mapping = pd.read_csv(MAPPING, index_col=0)["symbol"]
    local.index = local.index.to_series().map(mapping.to_dict()).fillna(local.index.to_series()).values
    local = local.groupby(level=0).mean()
    tcga = pd.read_csv(TCGA, index_col=0)
    tcga.index = tcga.index.astype(str)
    common = local.index.intersection(tcga.index)
    combined = pd.concat([local.loc[common], tcga.loc[common]], axis=1)
    combined = combined.apply(pd.to_numeric, errors="coerce").fillna(0)
    combined = (combined + 1).map(__import__("numpy").log2)
    top_genes = combined.var(axis=1).nlargest(min(2000, len(combined))).index
    values = combined.loc[top_genes].T
    scaled = StandardScaler().fit_transform(values)
    pca = PCA(n_components=2, random_state=42)
    coordinates = pca.fit_transform(scaled)
    local_ids = list(local.columns)
    tcga_ids = list(tcga.columns)
    samples = local_ids + tcga_ids
    coords = pd.DataFrame(coordinates, index=samples, columns=["PC1", "PC2"])
    coords["cohort"] = ["PDAC (n=11)"] * len(local_ids) + ["TCGA PAAD"] * len(tcga_ids)
    coords.to_csv(COORDS)

    fig, ax = plt.subplots(figsize=(7, 6))
    colors = {"PDAC (n=11)": "#1f77b4", "TCGA PAAD": "#ff7f0e"}
    for cohort, group in coords.groupby("cohort"):
        ax.scatter(group["PC1"], group["PC2"], s=60, color=colors[cohort], edgecolor="white", linewidth=0.5, alpha=0.85, label=cohort)
    for sample, row in coords[coords["cohort"] == "PDAC (n=11)"].iterrows():
        ax.annotate(sample.replace("T", ""), (row["PC1"], row["PC2"]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    ax.set_title("Joint RNA-seq PCA: revised PDAC versus TCGA PAAD\nAbsolute tumour expression on a common log2 scale")
    ax.legend(frameon=True, loc="lower left")
    ax.grid(True, alpha=0.3, linestyle=":")
    fig.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"common_genes={len(common)} features={len(top_genes)} samples={len(coords)} PC1={pca.explained_variance_ratio_[0]:.6f} PC2={pca.explained_variance_ratio_[1]:.6f} output={OUT}")


if __name__ == "__main__":
    main()
