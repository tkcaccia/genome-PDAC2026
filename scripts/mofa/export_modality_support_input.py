#!/usr/bin/env python3
# Numeric candidate support input for Figure S12; no additional panel is generated.
from pdac_runtime import private_path, project_root
from pathlib import Path


import pandas as pd


ROOT = private_path("handoff_root")


WORKSPACE_ROOT = project_root()


RUN = ROOT / "mofa" / "pdac_11"


EXPORTS = RUN / "mofa" / "exports"


RNA_MAP = WORKSPACE_ROOT / "reproducible_pipeline_share_20260713_with_outputs/raw_outputs/results/figures/categories/04_cross_cohort_concordance/rnaseq_ensembl_to_symbol.csv"


PROTEIN_MAP = WORKSPACE_ROOT / "reproducible_pipeline_share_20260713_with_outputs/raw_outputs/results/comparison_african_tcga/proteomics_index_to_symbol.csv"


OUT_TABLE = RUN / "tables/mofa_top200_composite_score_11_pdac.csv"


def scores(path):
    return pd.read_csv(path, index_col=0).abs().max(axis=1)


def symbols(series, path):
    mapping = pd.read_csv(path, index_col=0)["symbol"]
    mapped = series.index.to_series().map(mapping.to_dict()).fillna(series.index.to_series())
    series = series.copy()
    series.index = mapped.values
    return series.groupby(level=0).max()


def main():
    modality_scores = pd.DataFrame({
        "wes_score": scores(EXPORTS / "mofa_loadings_wes.csv"),
        "rnaseq_score": symbols(scores(EXPORTS / "mofa_loadings_rnaseq.csv"), RNA_MAP),
        "proteomics_score": symbols(scores(EXPORTS / "mofa_loadings_proteomics.csv"), PROTEIN_MAP),
    })
    percentiles = modality_scores.rank(pct=True)
    ranking = pd.DataFrame({
        "feature": modality_scores.index,
        "composite_score": percentiles.fillna(0).mean(axis=1) * 100,
        "n_modalities": modality_scores.notna().sum(axis=1),
        "dominant_modality": modality_scores.idxmax(axis=1).str.replace("_score", "", regex=False),
        "wes_supported": modality_scores["wes_score"].notna(),
        "proteomics_supported": modality_scores["proteomics_score"].notna(),
        "rnaseq_supported": modality_scores["rnaseq_score"].notna(),
    }).sort_values(["n_modalities", "composite_score"], ascending=[False, False]).reset_index(drop=True)
    ranking["rank"] = ranking.index + 1
    top = ranking.head(200)
    top.to_csv(OUT_TABLE, index=False)



if __name__ == "__main__":
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    main()
