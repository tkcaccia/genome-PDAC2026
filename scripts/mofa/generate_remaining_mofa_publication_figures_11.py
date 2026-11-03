#!/usr/bin/env python3
# Current collaborator code export. Three-view MOFA, exported factors/loadings, candidate ranking and sensitivity summaries.
# Source: Original Source and Code/scripts/generate_remaining_mofa_publication_figures_11.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Generate the remaining publication-oriented MOFA figures for the revised 11-patient cohort."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path
from itertools import combinations

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm

ROOT = project_root()
RUN = ROOT / "results/mofa_revised_11_20260916"
OUT = ROOT / "results/figures/final_11_figures"
MODEL = RUN / "mofa/mofa_model.hdf5"
PATIENT_CLUSTERS = RUN / "mofa/clusters/mofa_patient_clusters.csv"

PROGRAMMES = private_value("mofa_factor_annotations")


def training_data():
    wes = pd.read_csv(RUN / "integration/wes_binary.csv", index_col=0)
    rna = pd.read_csv(RUN / "integration/rnaseq_delta.csv", index_col=0)
    prot = pd.read_csv(RUN / "integration/proteomics_log2.csv", index_col=0)
    wes = wes.loc[wes.sum(axis=1) / wes.shape[1] >= 0.10]
    rna = rna.loc[rna.var(axis=1) > 0]
    prot = prot.loc[prot.notna().mean(axis=1) >= 0.20]
    wes = wes.loc[wes.var(axis=1).sort_values(ascending=False).head(600).index]
    zscore = lambda data: data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1).replace(0, 1), axis=0).fillna(0)
    rna = zscore(rna)
    prot = zscore(prot)
    rna = rna.loc[rna.var(axis=1).sort_values(ascending=False).head(2000).index]
    prot = prot.loc[prot.var(axis=1).sort_values(ascending=False).head(1200).index]
    wes = wes.sub(wes.mean(axis=1), axis=0)
    return {"WES": wes, "RNA": rna, "Proteomics": prot}


def load_model():
    with h5py.File(MODEL, "r") as model:
        z = np.asarray(model["expectations/Z/all"][()]).T
        factors = [f"factor_{index + 1}" for index in range(z.shape[1])]
        weights = {}
        for view, key in [("WES", "wes"), ("RNA", "rnaseq"), ("Proteomics", "proteomics")]:
            features = [x.decode() if isinstance(x, bytes) else str(x) for x in model[f"features/{key}"][()]]
            weights[view] = pd.DataFrame(np.asarray(model[f"expectations/W/{key}"][()]).T, index=features, columns=factors)
    return pd.DataFrame(z, columns=factors), weights


def native_r2(data, factors, weights):
    rows = []
    for view, matrix in data.items():
        y = matrix.to_numpy().T
        ss = np.square(y).sum()
        for factor in factors.columns:
            reconstruction = np.outer(factors[factor].values, weights[view][factor].values)
            rows.append({"view": view, "factor": factor, "r2_percent": (1 - np.square(y - reconstruction).sum() / ss) * 100})
    return pd.DataFrame(rows)


def variance_heatmap(data, factors, weights):
    r2 = native_r2(data, factors, weights)
    matrix = r2.pivot(index="factor", columns="view", values="r2_percent").reindex(index=factors.columns, columns=["WES", "RNA", "Proteomics"])
    r2.to_csv(RUN / "tables/mofa_native_variance_explained_11_pdac.csv", index=False)
    fig, ax = plt.subplots(figsize=(8, 6))
    limit = max(abs(matrix.min().min()), abs(matrix.max().max()))
    sns.heatmap(matrix, annot=True, fmt=".2f", cmap="vlag", center=0, vmin=-limit, vmax=limit, linewidths=.6, linecolor="white", cbar_kws={"label": "Native MOFA reconstruction R² (%)"}, ax=ax)
    ax.set_title("Native MOFA variance explained by factor and view\nRevised 11-patient PDAC cohort")
    ax.set_xlabel("Data view")
    ax.set_ylabel("Factor")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_18_mofa_native_variance_explained.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def loading_panels(weights):
    rows = []
    for view, table in weights.items():
        for factor in table.columns:
            selected = table[factor].abs().nlargest(10).index
            for feature in selected:
                rows.append({"view": view, "factor": factor, "feature": feature, "loading": table.loc[feature, factor]})
    top = pd.DataFrame(rows)
    top.to_csv(RUN / "tables/mofa_top_loadings_by_view_factor_11_pdac.csv", index=False)
    fig, axes = plt.subplots(5, 3, figsize=(17, 24), sharex=False)
    for row, factor in enumerate(weights["WES"].columns):
        for col, view in enumerate(["WES", "RNA", "Proteomics"]):
            ax = axes[row, col]
            subset = top[(top.factor == factor) & (top.view == view)].sort_values("loading")
            colours = ["#d1495b" if value < 0 else "#457b9d" for value in subset.loading]
            ax.barh(subset.feature.astype(str), subset.loading, color=colours)
            ax.axvline(0, color="black", linewidth=.7)
            ax.set_title(f"{factor.replace('_', ' ').title()} — {view}")
            ax.tick_params(axis="y", labelsize=8)
            ax.grid(axis="x", alpha=.25, linestyle=":")
    fig.suptitle("Top MOFA feature loadings by factor and data view\nRevised 11-patient PDAC cohort", fontsize=16, y=.995)
    fig.tight_layout(rect=(0, 0, 1, .985))
    fig.savefig(OUT / "Figure_19_mofa_top_loadings_by_view.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def factor_score_heatmap(factors):
    clusters = pd.read_csv(PATIENT_CLUSTERS, index_col=0)["cluster"].reindex(factors.index)
    ordered = factors.assign(cluster=clusters).sort_values("cluster").drop(columns="cluster")
    ordered.to_csv(RUN / "tables/mofa_factor_scores_annotated_11_pdac.csv")
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.heatmap(ordered.T, cmap="vlag", center=0, linewidths=.5, linecolor="white", cbar_kws={"label": "Factor score"}, ax=ax)
    ax.set_title("MOFA factor scores across revised 11-patient PDAC cohort")
    ax.set_xlabel("Patient")
    ax.set_ylabel("Factor")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_20_mofa_factor_scores_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def modality_support():
    ranking = pd.read_csv(RUN / "tables/mofa_top200_composite_score_11_pdac.csv")
    top = ranking.head(200).copy()
    flags = ["wes_supported", "rnaseq_supported", "proteomics_supported"]
    labels = {"wes_supported": "WES", "rnaseq_supported": "RNA", "proteomics_supported": "Proteomics"}
    top["support_set"] = top[flags].apply(lambda row: "+".join(labels[col] for col in flags if bool(row[col])), axis=1)
    counts = top["support_set"].value_counts().rename_axis("modality_support").reset_index(name="count")
    counts.to_csv(RUN / "tables/mofa_modality_support_top200_11_pdac.csv", index=False)
    ordered = counts.sort_values("count")
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(ordered.modality_support, ordered["count"], color="#4C78A8")
    for bar, value in zip(ax.patches, ordered["count"]):
        ax.text(value + .8, bar.get_y() + bar.get_height() / 2, str(value), va="center", fontsize=9)
    ax.set_xlabel("Number of candidates")
    ax.set_ylabel("Supporting modalities")
    ax.set_title("Modality support combinations among top 200 MOFA candidates")
    ax.grid(axis="x", alpha=.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_21_mofa_modality_support_top200.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = training_data()
    factors, weights = load_model()
    variance_heatmap(data, factors, weights)
    loading_panels(weights)
    factor_score_heatmap(factors)
    modality_support()
    print("Generated Figures 18-21 and supporting tables")


if __name__ == "__main__":
    main()
