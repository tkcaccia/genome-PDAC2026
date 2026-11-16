#!/usr/bin/env python3
# MOFA fitting, molecular integration and publication panels.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""
from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

ROOT = private_path("handoff_root")
RUN = ROOT / "mofa" / "pdac_11"
OUT = RUN / "figures"
TAB = RUN / "tables"
PATIENTS = [f"{p}T" for p in configured_cohort("primary_pdac")]
KRAS_HOTSPOT_POSITIONS = {"25245347", "25245350", "25245351", "25227342"}


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
    calls = pd.read_csv(ROOT / "wes" / "analysis" / "outputs" / "primary_pdac" / "wes_somatic_variants.tsv", sep="\t")
    kras = calls[calls["SYMBOL"].eq("KRAS") & calls["POS"].astype(str).isin(KRAS_HOTSPOT_POSITIONS)]
    mutant_patients = set(kras["patient"].astype(str))
    kras_mutant = [sample for sample in PATIENTS if sample[:-1] in mutant_patients]
    kras_wildtype = [sample for sample in PATIENTS if sample[:-1] not in mutant_patients]
    pd.DataFrame({
        "sample": PATIENTS,
        "patient": [sample[:-1] for sample in PATIENTS],
        "group": ["KRAS-hotspot mutant" if sample in kras_mutant else "No detected KRAS hotspot" for sample in PATIENTS],
    }).to_csv(TAB / "kras_hotspot_groups_11_pdac.csv", index=False)
    rows = []
    for view, matrix in [("RNA", rna), ("Proteomics", prot)]:
        for feature in matrix.index:
            mutant = matrix.loc[feature, kras_mutant].dropna()
            wild = matrix.loc[feature, kras_wildtype].dropna()
            if len(mutant) >= 3 and len(wild) >= 3:
                _, pvalue = mannwhitneyu(mutant, wild, alternative="two-sided", method="exact")
                rows.append({"view": view, "feature": feature, "mutant_mean": mutant.mean(), "wildtype_mean": wild.mean(), "difference_mutant_minus_wildtype": mutant.mean() - wild.mean(), "p_value": pvalue})
    result = pd.DataFrame(rows); result["FDR"] = result.groupby("view")["p_value"].transform(lambda x: bh(x.values)); result.sort_values(["view", "FDR", "p_value"]).to_csv(TAB / "kras_mutant_vs_wildtype_feature_comparison_11_pdac.csv", index=False)
    top = result.sort_values("p_value").head(30).copy(); top["label"] = top["feature"].astype(str).str.slice(0, 32); top["signed_minus_log10_p"] = -np.log10(top.p_value.clip(lower=1e-300)) * np.sign(top.difference_mutant_minus_wildtype); top = top.sort_values("signed_minus_log10_p")
    fig, ax = plt.subplots(figsize=(10, 9)); ax.barh(top.label, top.signed_minus_log10_p, color=top.view.map({"RNA":"#e76f51", "Proteomics":"#2a9d8f"})); ax.axvline(0, color="black", linewidth=.8); ax.set_xlabel("Signed -log10 nominal p-value"); ax.set_ylabel("Feature"); ax.set_title(f"KRAS-hotspot mutant versus no detected hotspot\nExploratory: n={len(kras_mutant)} versus n={len(kras_wildtype)}"); ax.grid(axis="x", alpha=.25, linestyle=":"); fig.tight_layout(); fig.savefig(OUT / "Figure_23_kras_mutant_vs_wildtype_features_11_pdac.png", dpi=300, bbox_inches="tight"); plt.close(fig)
    print("KRAS-hotspot patients:", [sample[:-1] for sample in kras_mutant], "non-mutant:", [sample[:-1] for sample in kras_wildtype])


def main():
    OUT.mkdir(parents=True, exist_ok=True); TAB.mkdir(parents=True, exist_ok=True)
    integrated_heatmap(); kras_stratified(); print("Generated integrated heatmap and KRAS-stratified outputs")

if __name__ == "__main__": main()
