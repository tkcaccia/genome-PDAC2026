#!/usr/bin/env python3
# MOFA fitting, molecular integration and publication panels.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""
from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

ROOT = private_path("handoff_root")
RUN = ROOT / "mofa/pdac_11"
FACTORS = RUN / "mofa/exports/mofa_factors.csv"
BASE_CLUSTERS = RUN / "mofa/clusters/mofa_patient_clusters.csv"
SUMMARY = RUN / "mofa/robustness/mofa_cluster_sensitivity_summary.csv"
OUT = RUN / "figures/Figure_10_cluster_robustness.png"

SCENARIOS = private_value("mofa_exclusion_scenarios")


def align_labels(base: pd.Series, labels: pd.Series) -> pd.Series:
    """Relabel a KMeans solution so its cluster IDs match the base clustering by majority overlap."""
    common = labels.index
    mapping = {}
    for cluster_id in labels.unique():
        members = labels[labels == cluster_id].index
        votes = base.loc[base.index.intersection(members)]
        mapping[cluster_id] = votes.mode().iat[0] if not votes.empty else cluster_id
    return labels.map(mapping)


def main() -> None:
    factors = pd.read_csv(FACTORS, index_col=0)
    base = pd.read_csv(BASE_CLUSTERS, index_col=0)["cluster"]
    summary = pd.read_csv(SUMMARY)

    assignments = pd.DataFrame(index=factors.index, columns=list(SCENARIOS))
    unchanged = True
    for scenario, drop_ids in SCENARIOS.items():
        subset = factors.drop(index=[x for x in drop_ids if x in factors.index], errors="ignore")
        labels = pd.Series(KMeans(n_clusters=2, random_state=42, n_init=20).fit_predict(subset), index=subset.index)
        aligned = align_labels(base, labels)
        assignments.loc[aligned.index, scenario] = aligned.values
        ari = adjusted_rand_score(base.loc[labels.index], labels)
        recorded = summary.loc[summary["scenario"].eq(scenario), "ari_vs_full"].iat[0]
        if abs(ari - recorded) > 1e-10:
            raise ValueError(f"Recomputed ARI disagrees with summary for {scenario}")
        unchanged = unchanged and aligned.eq(base.loc[aligned.index]).all()

    plot_data = assignments.reindex(base.sort_values().index).astype(float)
    plot_data.to_csv(RUN / "tables/Figure_10_cluster_assignments.csv", index_label="patient")

    fig, (ax_heat, ax_sil) = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={"height_ratios": [3, 1]}, sharex=False)
    cmap = ListedColormap(["#66c2a5", "#fc8d62"])
    im = ax_heat.imshow(plot_data.values, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax_heat.set_xticks(range(len(SCENARIOS)))
    ax_heat.set_xticklabels(list(SCENARIOS), rotation=20, ha="right")
    ax_heat.set_yticks(range(len(plot_data)))
    ax_heat.set_yticklabels(plot_data.index)
    ax_heat.set_ylabel("Patient")
    for row in range(plot_data.shape[0]):
        for col in range(plot_data.shape[1]):
            value = plot_data.values[row, col]
            if pd.notna(value):
                ax_heat.text(col, row, int(value), ha="center", va="center", fontsize=8)
    subtitle = "All retained patients keep their cluster label; ARI = 1.0 in every scenario" if unchanged else "Cluster labels aligned to the full-cohort solution"
    ax_heat.set_title(f"Per-patient MOFA cluster assignment across exclusion scenarios\n{subtitle}")
    cbar = fig.colorbar(im, ax=ax_heat, ticks=[0, 1], fraction=0.046, pad=0.04)
    cbar.ax.set_yticklabels(["Cluster 0", "Cluster 1"])

    ax_sil.plot(summary["scenario"], summary["silhouette"], marker="o", color="#4B3F8C")
    for x, y, n in zip(summary["scenario"], summary["silhouette"], summary["n_patients"]):
        ax_sil.annotate(f"n={n}", (x, y), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=8)
    ax_sil.set_ylabel("Silhouette score")
    ax_sil.set_xticks(range(len(summary)))
    ax_sil.set_xticklabels(summary["scenario"], rotation=20, ha="right")
    ax_sil.set_ylim(0.35, 0.45)
    ax_sil.grid(axis="y", alpha=0.3, linestyle="--")

    fig.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUT}")


if __name__ == "__main__":
    main()
