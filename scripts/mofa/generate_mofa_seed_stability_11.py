#!/usr/bin/env python3
"""Assess future seed fits against the saved primary model on identical inputs.

This correction does not update the manuscript's existing stability results.
A scientific rerun is required before replacing those results.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import linear_sum_assignment
from mofapy2.run.entry_point import entry_point

from pdac_runtime import private_path
from scripts.mofa.model_inputs import VIEWS, load_saved_inputs, write_input_identity

RUN = private_path("handoff_root") / "mofa" / "pdac_11"
# Never overwrite the supplied fit-to-fit stability tables or manuscript figure.
OUT = RUN / "mofa" / "seed_stability_model_aligned"
SEEDS = (42, 43, 44)


def fit(seed, matrices):
    model = entry_point()
    data = [[matrices[view].T.to_numpy()] for view in VIEWS]
    model.set_data_matrix(
        data=data, likelihoods=["bernoulli", "gaussian", "gaussian"],
        views_names=list(VIEWS), groups_names=["all"],
        samples_names=[list(matrices["wes"].columns)],
        features_names=[matrices[view].index.tolist() for view in VIEWS],
    )
    model.set_data_options(scale_views=False, scale_groups=False,
                           center_groups=True, use_float32=False)
    model.set_model_options(factors=5, spikeslab_factors=False,
                            spikeslab_weights=True, ard_factors=False, ard_weights=True)
    model.set_train_options(
        iter=1500, startELBO=1, freqELBO=5, startSparsity=50, tolerance=1e-5,
        convergence_mode="fast", startDrop=1, freqDrop=5,
        verbose=False, quiet=True, seed=seed,
    )
    model.build()
    model.run()
    return np.asarray(model.model.nodes["Z"].getExpectations()["E"])


def match_factors(reference, scores):
    if scores.shape != reference.shape or not np.isfinite(scores).all():
        raise ValueError("Fit factor-score shape/values differ from the saved reference")
    k = reference.shape[1]
    correlations = np.corrcoef(reference.T, scores.T)[:k, k:]
    if not np.isfinite(correlations).all():
        raise ValueError("Constant/invalid factors prevent correlation-based matching")
    rows, columns = linear_sum_assignment(-np.abs(correlations))
    order = np.argsort(rows)
    rows, columns = rows[order], columns[order]
    signs = np.where(correlations[rows, columns] < 0, -1, 1)
    return scores[:, columns] * signs, columns


def main():
    matrices, saved_factors = load_saved_inputs(RUN)
    reference = saved_factors.to_numpy()
    if reference.shape[1] != 5:
        raise ValueError("Expected the five-factor primary model")
    OUT.mkdir(parents=True, exist_ok=False)
    write_input_identity(OUT / "input_identity.json", matrices, saved_factors)
    rows = []
    for seed in SEEDS:
        aligned, columns = match_factors(reference, fit(seed, matrices))
        for index, factor in enumerate(saved_factors.columns):
            rows.append({
                "seed": seed, "factor": factor,
                "matched_run_factor": f"factor_{columns[index] + 1}",
                "correlation_vs_saved_primary": np.corrcoef(
                    reference[:, index], aligned[:, index])[0, 1],
            })
        pd.DataFrame(aligned, index=saved_factors.index,
                     columns=saved_factors.columns).to_csv(OUT / f"aligned_scores_seed_{seed}.csv")
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "stability_vs_saved_primary.csv", index=False)
    pivot = result.pivot(index="factor", columns="seed", values="correlation_vs_saved_primary")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    sns.heatmap(pivot, annot=True, fmt=".3f", cmap="vlag", center=0,
                vmin=-1, vmax=1, ax=ax,
                cbar_kws={"label": "Correlation with saved primary model"})
    ax.set_title("MOFA stability on identical saved-model features")
    ax.set_xlabel("Seed")
    ax.set_ylabel("Factor")
    fig.tight_layout()
    fig.savefig(OUT / "stability_vs_saved_primary.png", dpi=300)
    plt.close(fig)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
