#!/usr/bin/env python3
# PDAC2026 analysis: Three-view MOFA, exported factors/loadings, candidate ranking and sensitivity summaries.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Assess revised 11-patient MOFA stability across independent random seeds."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr
from mofapy2.run.entry_point import entry_point

ROOT = project_root()
RUN = ROOT / "results/mofa_revised_11_20260916"
OUT = RUN / "mofa/seed_stability"
FIG = ROOT / "results/figures/final_11_figures/Figure_24_mofa_seed_stability.png"
SEEDS = [42, 43, 44]
VIEWS = ["wes", "rnaseq", "proteomics"]


def load_inputs():
    matrices = {}
    for view, filename in [("wes", "wes_binary.csv"), ("rnaseq", "rnaseq_delta.csv"), ("proteomics", "proteomics_log2.csv")]:
        matrix = pd.read_csv(RUN / "integration" / filename, index_col=0)
        if view == "wes": matrix = matrix.loc[matrix.sum(axis=1) / matrix.shape[1] >= .10]
        elif view == "rnaseq": matrix = matrix.loc[matrix.var(axis=1) > 0]
        else: matrix = matrix.loc[matrix.notna().mean(axis=1) >= .20]
        matrix = matrix.loc[matrix.var(axis=1).sort_values(ascending=False).head({"wes": 600, "rnaseq": 2000, "proteomics": 1200}[view]).index]
        if view != "wes": matrix = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, 1), axis=0).fillna(0)
        else: matrix = matrix.astype(float)
        matrices[view] = matrix
    return matrices


def fit(seed, matrices):
    model = entry_point()
    data = [[matrices[view].T.to_numpy()] for view in VIEWS]
    model.set_data_matrix(data=data, likelihoods=["bernoulli", "gaussian", "gaussian"], views_names=VIEWS, groups_names=["all"], samples_names=[list(matrices["wes"].columns)], features_names=[matrices[v].index.tolist() for v in VIEWS])
    model.set_data_options(scale_views=False, scale_groups=False, center_groups=True, use_float32=False)
    model.set_model_options(factors=5, spikeslab_factors=False, spikeslab_weights=True, ard_factors=False, ard_weights=True)
    model.set_train_options(iter=1500, startELBO=1, freqELBO=5, startSparsity=50, tolerance=1e-5, convergence_mode="fast", startDrop=1, freqDrop=5, verbose=False, quiet=True, seed=seed)
    model.build(); model.run()
    z = np.asarray(model.model.nodes["Z"].getExpectations()["E"])
    w = {view: np.asarray(model.model.nodes["W"].getNodes()[index].getExpectations()["E"]).T for index, view in enumerate(VIEWS)}
    return z, w


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    matrices = load_inputs(); fits = {seed: fit(seed, matrices) for seed in SEEDS}
    reference = fits[42][0]
    rows = []
    aligned_scores = {}
    for seed, (scores, weights) in fits.items():
        corr = np.corrcoef(reference.T, scores.T)[:5, 5:]
        row_ind, col_ind = linear_sum_assignment(-np.abs(corr))
        aligned = scores[:, col_ind] * np.sign(corr[row_ind, col_ind])
        aligned_scores[seed] = aligned
        for ref_factor, run_factor, value in zip(row_ind + 1, col_ind + 1, corr[row_ind, col_ind]):
            rows.append({"seed": seed, "reference_factor": f"factor_{ref_factor}", "run_factor": f"factor_{run_factor}", "sample_score_correlation": value})
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "mofa_seed_stability_summary_11_pdac.csv", index=False)
    heat = pd.DataFrame(index=[f"seed_{s}_factor_{i+1}" for s in SEEDS for i in range(5)])
    for seed, scores in aligned_scores.items():
        for i in range(5):
            heat.loc[f"seed_{seed}_factor_{i+1}", "seed_42_factor_1"] = np.nan
    pair_rows=[]
    for seed in SEEDS[1:]:
        for i in range(5):
            pair_rows.append({"seed": seed, "factor": f"factor_{i+1}", "correlation_vs_seed42": pearsonr(aligned_scores[42][:,i], aligned_scores[seed][:,i]).statistic})
    pair = pd.DataFrame(pair_rows)
    pair.to_csv(OUT / "mofa_factor_score_stability_vs_seed42_11_pdac.csv", index=False)
    pivot = pair.pivot(index="factor", columns="seed", values="correlation_vs_seed42")
    fig, ax = plt.subplots(figsize=(7, 4.5)); sns.heatmap(pivot, annot=True, fmt=".3f", cmap="vlag", center=0, vmin=-1, vmax=1, ax=ax, cbar_kws={"label":"Correlation with seed 42"}); ax.set_title("MOFA factor-score stability across random seeds\nRevised 11-patient PDAC cohort"); ax.set_xlabel("Seed"); ax.set_ylabel("Factor"); fig.tight_layout(); fig.savefig(FIG, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(pair.to_string(index=False))

if __name__ == "__main__": main()
