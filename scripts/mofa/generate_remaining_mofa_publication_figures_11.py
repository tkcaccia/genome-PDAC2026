#!/usr/bin/env python3
# MOFA fitting, molecular integration and publication panels.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""
from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path
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

ROOT = private_path("handoff_root")
RUN = ROOT / "mofa" / "pdac_11"
OUT = RUN / "figures"
MODEL = RUN / "mofa" / "mofa_model.hdf5"
PATIENT_CLUSTERS = RUN / "mofa" / "clusters" / "mofa_patient_clusters.csv"

PROGRAMMES = private_value("mofa_factor_annotations")


def training_data(factors, weights):
    data = {
        "WES": pd.read_csv(RUN / "qc/wes_binary_qc.csv", index_col=0),
        "RNA": pd.read_csv(RUN / "qc/rnaseq_delta_zscore.csv", index_col=0),
        "Proteomics": pd.read_csv(RUN / "qc/proteomics_log2_zscore.csv", index_col=0),
    }
    aligned = {}
    for view, matrix in data.items():
        if view == "WES":
            matrix = matrix.sub(matrix.mean(axis=1), axis=0)
        aligned[view] = matrix.loc[weights[view].index, factors.index]
    return aligned


def load_model():
    with h5py.File(MODEL, "r") as model:
        z = np.asarray(model["expectations/Z/all"][()]).T
        samples = [x.decode() if isinstance(x, bytes) else str(x) for x in model["samples/all"][()]]
        factors = [f"factor_{index + 1}" for index in range(z.shape[1])]
        weights = {}
        for view, key in [("WES", "wes"), ("RNA", "rnaseq"), ("Proteomics", "proteomics")]:
            features = [x.decode() if isinstance(x, bytes) else str(x) for x in model[f"features/{key}"][()]]
            weights[view] = pd.DataFrame(np.asarray(model[f"expectations/W/{key}"][()]).T, index=features, columns=factors)
    return pd.DataFrame(z, index=samples, columns=factors), weights


def native_r2(data, factors, weights):
    rows = []
    z = factors.to_numpy()
    for view, matrix in data.items():
        y = matrix.to_numpy().T
        ss = np.square(y).sum()
        w = weights[view].to_numpy()
        for index, factor in enumerate(factors.columns):
            reconstruction = np.outer(z[:, index], w[:, index])
            rows.append({"view": view, "factor": factor, "r2_percent": (1 - np.square(y - reconstruction).sum() / ss) * 100})
        reconstruction = z @ w.T
        rows.append({"view": view, "factor": "All five factors", "r2_percent": (1 - np.square(y - reconstruction).sum() / ss) * 100})
    return pd.DataFrame(rows)


def variance_heatmap(data, factors, weights):
    r2 = native_r2(data, factors, weights)
    factor_order = list(factors.columns) + ["All five factors"]
    matrix = r2.pivot(index="factor", columns="view", values="r2_percent").reindex(index=factor_order, columns=["WES", "RNA", "Proteomics"])
    matrix.index = [factor.replace("factor_", "Factor ") if factor.startswith("factor_") else factor for factor in matrix.index]
    r2.to_csv(RUN / "tables/mofa_native_variance_explained_11_pdac.csv", index=False)
    fig, ax = plt.subplots(figsize=(8, 6.8))
    limit = max(abs(matrix.min().min()), abs(matrix.max().max()))
    sns.heatmap(matrix, annot=True, fmt=".2f", cmap="vlag", center=0, vmin=-limit, vmax=limit, linewidths=.6, linecolor="white", cbar_kws={"label": "Native MOFA+ reconstruction R² (%)"}, ax=ax)
    ax.set_title("Centered training-set reconstruction by MOFA+ factor\nRevised 11-patient PDAC cohort")
    ax.set_xlabel("Data view")
    ax.set_ylabel("Factor")
    fig.tight_layout()
    fig.savefig(OUT / "Figure_18_mofa_native_variance_explained.png", dpi=300, bbox_inches="tight")
    plt.close(fig)






def modality_support():
    ranking = pd.read_csv(RUN / "mofa/biology/mofa_cross_modal_target_ranking.csv")
    top = ranking.head(200).copy()
    score_columns = {"wes_score": "WES", "rnaseq_score": "RNA", "proteomics_score": "Proteomics"}
    flags = list(score_columns)
    labels = score_columns
    for column in flags:
        top[column] = top[column].notna()
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
    factors, weights = load_model()
    data = training_data(factors, weights)
    variance_heatmap(data, factors, weights)
    pass
    pass
    modality_support()
    print("Generated Figures 18-21 and supporting tables")


if __name__ == "__main__":
    main()
