#!/usr/bin/env python3
# PDAC2026 analysis: Salmon gene counts/TPM, paired DESeq2, response PCA and enrichment.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Generate a publication-standard PCA of revised 11-patient RNA response profiles."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = project_root()
RNA = ROOT / "results/rnaseq_revise_11_archive/rnaseq_revised_11_20260916/tables/tumour_normal_log2_delta_matrix_11_pdac.csv"
OUT = ROOT / "results/figures/final_11_figures/Figure_13_rna_pca.png"
COORDS = ROOT / "results/rnaseq_revised_11_native_deseq2_20260917/rna_response_pca_coordinates_11_pdac.csv"


def main() -> None:
    matrix = pd.read_csv(RNA, index_col=0)
    matrix = matrix.loc[:, configured_cohort("primary_pdac")]
    top_genes = matrix.var(axis=1).nlargest(min(2000, len(matrix))).index
    values = matrix.loc[top_genes].T
    scaled = StandardScaler().fit_transform(values)
    model = PCA(n_components=2, random_state=42)
    coordinates = model.fit_transform(scaled)
    coords = pd.DataFrame(coordinates, index=values.index, columns=["PC1", "PC2"])
    coords.to_csv(COORDS)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(coords["PC1"], coords["PC2"], color="#4B3F8C", s=90, edgecolor="black", linewidth=0.7, alpha=0.9)
    offsets = private_value("pca_label_offsets", default={})
    for patient, row in coords.iterrows():
        dx, dy = offsets.get(patient, (7, 7))
        ax.annotate(patient, (row["PC1"], row["PC2"]), xytext=(dx, dy), textcoords="offset points", fontsize=9, color="black")
    ax.set_xlabel(f"PC1 ({model.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({model.explained_variance_ratio_[1] * 100:.1f}%)")
    ax.set_title("RNA-seq tumour-normal response PCA\nRevised 11-patient PDAC cohort")
    ax.grid(True, alpha=0.3, linestyle=":")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"features={len(top_genes)} PC1={model.explained_variance_ratio_[0]:.6f} PC2={model.explained_variance_ratio_[1]:.6f} output={OUT}")


if __name__ == "__main__":
    main()
