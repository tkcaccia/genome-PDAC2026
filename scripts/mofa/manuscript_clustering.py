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
from sklearn.metrics import silhouette_score


ROOT = private_path("handoff_root")
RUN = ROOT / "mofa/pdac_11"


def main():
    factors = pd.read_csv(RUN / "mofa/exports/mofa_factors.csv", index_col=0)
    clusters = pd.read_csv(RUN / "mofa/clusters/mofa_patient_clusters.csv", index_col=0)["cluster"].reindex(factors.index)
    if clusters.isna().any() or clusters.nunique() != 2 or len(factors) != 11:
        raise ValueError("Expected aligned 11-patient factor scores and two clusters")
    silhouette = silhouette_score(factors, clusters)
    isolated = clusters.value_counts().idxmin()
    display_clusters = clusters.map(lambda label: 0 if label == isolated else 1)
    fig, ax = plt.subplots(figsize=(8, 7))
    for cluster, color in [(0, "#66c2a5"), (1, "#fc8d62")]:
        mask = display_clusters.eq(cluster)
        ax.scatter(factors.loc[mask, "factor_1"], factors.loc[mask, "factor_2"],
                   color=color, s=120, edgecolor="black", linewidth=.7,
                   alpha=.9, label=str(cluster), zorder=3)
    offsets = private_value("pca_label_offsets", {})
    for sample, row in factors.iterrows():
        ax.annotate(sample, (row["factor_1"], row["factor_2"]),
                    xytext=offsets.get(sample, (7, 7)), textcoords="offset points",
                    fontsize=10, color="black", ha="left", va="center")
    ax.set_xlabel("Factor 1")
    ax.set_ylabel("Factor 2")
    ax.set_title(f"MOFA factor clustering (k=2, silhouette={silhouette:.3f})", fontsize=14, fontweight="bold")
    ax.legend(title="Cluster", loc="lower left", frameon=True)
    ax.grid(True, linewidth=.6, alpha=.5)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    output = RUN / "figures/Figure_9_mofa_factor_clustering.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    table = factors.copy()
    table["model_cluster"] = clusters
    table["display_cluster"] = display_clusters
    table.to_csv(RUN / "tables/Figure_9_mofa_factor_clustering.csv", index_label="patient")
    print(f"Saved {output}; silhouette={silhouette:.6f}")
    print("Displayed Cluster 0:", ", ".join(display_clusters[display_clusters.eq(0)].index))


if __name__ == "__main__":
    main()