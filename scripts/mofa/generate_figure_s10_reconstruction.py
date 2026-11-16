#!/usr/bin/env python3
# MOFA fitting, molecular integration and publication panels.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""
from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = private_path("handoff_root")
RUN = ROOT / "mofa" / "pdac_11"
OUT = ROOT / "svg_pdf" / "s10_rerun"
MODEL = RUN / "mofa" / "mofa_model.hdf5"


def decode(values):
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def load_model():
    with h5py.File(MODEL, "r") as model:
        z = np.asarray(model["expectations/Z/all"][()]).T
        samples = decode(model["samples/all"][()])
        factor_names = [f"factor_{index + 1}" for index in range(z.shape[1])]
        factors = pd.DataFrame(z, index=samples, columns=factor_names)
        weights = {}
        for view, key in [("WES", "wes"), ("RNA", "rnaseq"), ("Proteomics", "proteomics")]:
            features = decode(model[f"features/{key}"][()])
            values = np.asarray(model[f"expectations/W/{key}"][()]).T
            weights[view] = pd.DataFrame(values, index=features, columns=factor_names)
    return factors, weights


def load_centered_training_data(factors, weights):
    matrices = {
        "WES": pd.read_csv(RUN / "qc" / "wes_binary_qc.csv", index_col=0),
        "RNA": pd.read_csv(RUN / "qc" / "rnaseq_delta_zscore.csv", index_col=0),
        "Proteomics": pd.read_csv(RUN / "qc" / "proteomics_log2_zscore.csv", index_col=0),
    }
    aligned = {}
    for view, matrix in matrices.items():
        if view == "WES":
            matrix = matrix.sub(matrix.mean(axis=1), axis=0)
        aligned[view] = matrix.loc[weights[view].index, factors.index]
    return aligned


def calculate_r2(data, factors, weights):
    rows = []
    z = factors.to_numpy()
    for view, matrix in data.items():
        y = matrix.to_numpy().T
        denominator = np.square(y).sum()
        w = weights[view].to_numpy()
        for index, factor in enumerate(factors.columns):
            reconstruction = np.outer(z[:, index], w[:, index])
            r2 = 1 - np.square(y - reconstruction).sum() / denominator
            rows.append({"view": view, "factor": factor, "r2_percent": 100 * r2})
        reconstruction = z @ w.T
        r2 = 1 - np.square(y - reconstruction).sum() / denominator
        rows.append({"view": view, "factor": "All five factors", "r2_percent": 100 * r2})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    factors, weights = load_model()
    data = load_centered_training_data(factors, weights)
    results = calculate_r2(data, factors, weights)
    results.to_csv(OUT / "mofa_native_variance_explained_11_pdac.csv", index=False)

    order = list(factors.columns) + ["All five factors"]
    matrix = results.pivot(index="factor", columns="view", values="r2_percent")
    matrix = matrix.reindex(index=order, columns=["WES", "RNA", "Proteomics"])
    matrix.index = [name.replace("factor_", "Factor ") if name.startswith("factor_") else name for name in matrix.index]
    limit = max(abs(matrix.min().min()), abs(matrix.max().max()))

    fig, ax = plt.subplots(figsize=(8, 6.8))
    sns.heatmap(
        matrix,
        annot=True,
        fmt=".2f",
        cmap="vlag",
        center=0,
        vmin=-limit,
        vmax=limit,
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"label": "Native MOFA+ reconstruction R² (%)"},
        ax=ax,
    )
    ax.set_title("Centered training-set reconstruction by MOFA+ factor\nRevised 11-patient PDAC cohort")
    ax.set_xlabel("Data view")
    ax.set_ylabel("Factor")
    fig.tight_layout()
    for extension in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"Figure_S10_native_reconstruction.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(results[results.factor.eq("All five factors")].to_string(index=False))


if __name__ == "__main__":
    main()