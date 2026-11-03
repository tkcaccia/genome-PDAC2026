#!/usr/bin/env python3
# PDAC2026 analysis: Three-view MOFA, exported factors/loadings, candidate ranking and sensitivity summaries.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Rank cross-modal candidates by their composite evidence score."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = project_root()
EXPORTS = ROOT / "results/mofa_revised_11_20260916/mofa/exports"
RNA_MAP = ROOT / "reproducible_pipeline_share_20260713_with_outputs/raw_outputs/results/figures/categories/04_cross_cohort_concordance/rnaseq_ensembl_to_symbol.csv"
PROTEIN_MAP = ROOT / "reproducible_pipeline_share_20260713_with_outputs/raw_outputs/results/comparison_african_tcga/proteomics_index_to_symbol.csv"
OUT_FIG = ROOT / "results/figures/final_11_figures/Figure_11_top_20_candidates.png"
OUT_TABLE = ROOT / "results/mofa_revised_11_20260916/tables/mofa_top20_composite_score_11_pdac.csv"

MODALITY_COLORS = {"wes": "#2ecc71", "rnaseq": "#3498db", "proteomics": "#e74c3c"}


def max_abs_loading(path: Path) -> pd.Series:
    loadings = pd.read_csv(path, index_col=0)
    return loadings.abs().max(axis=1)


def to_symbols(scores: pd.Series, mapping_path: Path) -> pd.Series:
    """Map an ID-indexed score series to gene symbols so it can be merged across modalities."""
    mapping = pd.read_csv(mapping_path, index_col=0)["symbol"]
    mapped_index = scores.index.to_series().map(mapping.to_dict()).fillna(scores.index.to_series())
    scores = scores.copy()
    scores.index = mapped_index.values
    return scores.groupby(level=0).max()


def main() -> None:
    scores = pd.DataFrame({
        "wes_score": max_abs_loading(EXPORTS / "mofa_loadings_wes.csv"),
        "rnaseq_score": to_symbols(max_abs_loading(EXPORTS / "mofa_loadings_rnaseq.csv"), RNA_MAP),
        "proteomics_score": to_symbols(max_abs_loading(EXPORTS / "mofa_loadings_proteomics.csv"), PROTEIN_MAP),
    })
    # Percentile-normalize each modality independently so features are comparable across views.
    pct = scores.rank(pct=True)
    n_modalities = scores.notna().sum(axis=1)
    composite = pct.fillna(0).mean(axis=1) * 100
    dominant = scores.idxmax(axis=1).str.replace("_score", "", regex=False)

    ranking = pd.DataFrame({
        "feature": scores.index,
        "composite_score": composite,
        "n_modalities": n_modalities,
        "dominant_modality": dominant,
    }).sort_values(["n_modalities", "composite_score"], ascending=[False, False]).reset_index(drop=True)
    top20 = ranking.head(20).iloc[::-1]
    top20.to_csv(OUT_TABLE, index=False)

    fig, ax = plt.subplots(figsize=(10, 8))
    colors = [MODALITY_COLORS.get(m, "#888888") for m in top20["dominant_modality"]]
    ax.barh(top20["feature"], top20["composite_score"], color=colors)
    ax.set_xlabel("Composite Evidence Score")
    ax.set_title("Top 20 Ranked Biomarker Candidates: 11 PDAC patients\n(color = dominant modality; not a biological rank by modality)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in MODALITY_COLORS.values()]
    ax.legend(handles, [m.upper() for m in MODALITY_COLORS], title="Modality", loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(top20.to_string(index=False))


if __name__ == "__main__":
    main()
