#!/usr/bin/env python3
# PDAC2026 analysis: Three-view MOFA, exported factors/loadings, candidate ranking and sensitivity summaries.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Generate KRAS-stratified and integrated feature figures for revised 11-PDAC cohort."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

ROOT = project_root()
RUN = ROOT / "results/mofa_revised_11_20260916"
OUT = ROOT / "results/figures/final_11_figures"
TAB = RUN / "tables"
PATIENTS = [f"{patient}T" for patient in configured_cohort("primary_pdac")]
KRAS_MUTANT = [f"{patient}T" for patient in configured_cohort("kras_hotspot_retained")]
KRAS_WT = [f"{patient}T" for patient in configured_cohort("kras_hotspot_not_retained")]


def load_matrix(name):
    return pd.read_csv(RUN / "integration" / name, index_col=0).loc[:, PATIENTS]


def bh(pvalues):
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values)
    adjusted = np.empty_like(values)
    adjusted[order] = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    return np.clip(adjusted, 0, 1)


def integrated_heatmap():
    wes = load_matrix("wes_binary.csv")
    rna = load_matrix("rnaseq_delta.csv")
    prot = load_matrix("proteomics_log2.csv")
    feature_sets = {
        "KRAS": ["KRAS"],
        "TP53": ["TP53"],
        "SMAD4": ["SMAD4"],
        "CDKN2A": ["CDKN2A"],
        "GNAS": ["GNAS"],
        "RNF43": ["RNF43"],
    }
    candidates = []
    for name, features in feature_sets.items():
        for feature in features:
            if feature in wes.index: candidates.append(("WES", feature, wes.loc[feature]))
            if feature in rna.index: candidates.append(("RNA", feature, rna.loc[feature]))
            if feature in prot.index: candidates.append(("Proteomics", feature, prot.loc[feature]))
    top_rna = rna.var(axis=1).nlargest(12).index
    for feature in top_rna: candidates.append(("RNA", str(feature), rna.loc[feature]))
    top_prot = prot.var(axis=1).nlargest(12).index
    for feature in top_prot: candidates.append(("Proteomics", str(feature), prot.loc[feature]))
    matrix = pd.DataFrame({f"{view}:{feature}": values for view, feature, values in candidates}).T
    matrix = matrix.loc[~matrix.index.duplicated()]
    matrix.to_csv(TAB / "kras_integrated_feature_heatmap_11_pdac.csv")
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, 1), axis=0)
    fig, ax = plt.subplots(figsize=(13, 10))
    sns.heatmap(z, cmap="vlag", center=0, linewidths=.3, linecolor="white", ax=ax, cbar_kws={"label": "Feature z-score across patients"})
    ax.set_title("Integrated mutation, RNA-response and proteomics features\nRevised 11-patient PDAC cohort")
    ax.set_xlabel("Patient")
    ax.set_ylabel("Feature / modality")
    fig.tight_layout(); fig.savefig(OUT / "Figure_22_integrated_feature_heatmap_11_pdac.png", dpi=300, bbox_inches="tight"); plt.close(fig)


def kras_stratified():
    rna = load_matrix("rnaseq_delta.csv")
    prot = load_matrix("proteomics_log2.csv")
    rows = []
    for view, matrix in [("RNA", rna), ("Proteomics", prot)]:
        for feature in matrix.index:
            mutant = matrix.loc[feature, KRAS_MUTANT].dropna()
            wild = matrix.loc[feature, KRAS_WT].dropna()
            if len(mutant) >= 3 and len(wild) >= 3:
                _, pvalue = mannwhitneyu(mutant, wild, alternative="two-sided", method="exact")
                rows.append({"view": view, "feature": feature, "mutant_mean": mutant.mean(), "wildtype_mean": wild.mean(), "difference_mutant_minus_wildtype": mutant.mean() - wild.mean(), "p_value": pvalue})
    result = pd.DataFrame(rows); result["FDR"] = result.groupby("view")["p_value"].transform(lambda x: bh(x.values)); result.sort_values(["view", "FDR", "p_value"]).to_csv(TAB / "kras_mutant_vs_wildtype_feature_comparison_11_pdac.csv", index=False)
    top = result.sort_values("p_value").head(30).copy(); top["label"] = top["feature"].astype(str).str.slice(0, 32); top["signed_minus_log10_p"] = -np.log10(top.p_value.clip(lower=1e-300)) * np.sign(top.difference_mutant_minus_wildtype); top = top.sort_values("signed_minus_log10_p")
    fig, ax = plt.subplots(figsize=(10, 9)); ax.barh(top.label, top.signed_minus_log10_p, color=top.view.map({"RNA":"#e76f51", "Proteomics":"#2a9d8f"})); ax.axvline(0, color="black", linewidth=.8); ax.set_xlabel("Signed -log10 nominal p-value"); ax.set_ylabel("Feature"); ax.set_title("KRAS-mutant versus KRAS-wild-type RNA/proteomics features\nExploratory: n=5 versus n=6"); ax.grid(axis="x", alpha=.25, linestyle=":"); fig.tight_layout(); fig.savefig(OUT / "Figure_23_kras_mutant_vs_wildtype_features_11_pdac.png", dpi=300, bbox_inches="tight"); plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True); TAB.mkdir(parents=True, exist_ok=True)
    integrated_heatmap(); kras_stratified(); print("Generated integrated heatmap and KRAS-stratified outputs")

if __name__ == "__main__": main()
