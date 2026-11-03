#!/usr/bin/env python3
# Current collaborator code export. Final supplied vector presentation scripts, consuming private derived tables.
# Source: PDAC_Manuscript_Figs_Sne/generate_study_vectors.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Regenerate manuscript plots as native SVG/PDF from archived study tables."""

from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


ROOT = figure_root()
DATA = private_path("study_tables")
PATIENTS = configured_cohort("primary_pdac")
SAMPLES = [f"{patient}T" for patient in PATIENTS]
plt.rcParams["svg.fonttype"] = "none"


def table(relative: str, **kwargs) -> pd.DataFrame:
    return pd.read_csv(DATA / relative, **kwargs)


def export(relative: str, fig: plt.Figure) -> None:
    base = ROOT / relative
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".svg"), format="svg", bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
    plt.close(fig)


def wes_figures() -> None:
    burden = table("WES/wes_mutation_burden.csv")
    burden["patient"] = burden["patient"].astype(str)
    burden = burden.set_index("patient").reindex(PATIENTS)
    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.bar(burden.index, burden["somatic_variant_count"], color="#264653", edgecolor="#1d3557")
    ax.set(xlabel="Patient", ylabel="Filtered somatic variants", title="WES somatic mutation burden: revised 11-patient PDAC cohort")
    for bar, value in zip(bars, burden["somatic_variant_count"]):
        ax.annotate(f"{int(value):,}", (bar.get_x() + bar.get_width() / 2, value), xytext=(0, 4), textcoords="offset points", ha="center", fontsize=8)
    export("Figure 1/Figure 1A_mutation burden per patient", fig)

    matrix = table("WES/wes_mutation_matrix.blacklist_filtered.csv", index_col=0).fillna(0)
    matrix.columns = [str(c).lstrip("0").rstrip("T") for c in matrix.columns]
    matrix = matrix.reindex(columns=PATIENTS, fill_value=0)
    top = matrix.sum(axis=1).sort_values(ascending=False).head(30).index
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(matrix.loc[top].astype(int), cmap=["#eeeeee", "#4b3f8c"], cbar=False, linewidths=.4, linecolor="white", ax=ax)
    ax.set(xlabel="Patient", ylabel="Gene", title="Recurrent mutation oncoprint: revised 11-patient PDAC cohort")
    export("Figure 1/Figure 1B oncoprint", fig)

    genes = table("WES/gene_summary_table.csv").sort_values("variant_count", ascending=False).head(20).sort_values("variant_count")
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.barh(genes["SYMBOL"], genes["variant_count"], color=plt.cm.magma(np.linspace(.12, .9, len(genes))))
    ax.set(xlabel="Total variants", ylabel="Gene", title="Top 20 genes by variant count: revised 11-patient PDAC cohort")
    export("Figure 1/Figure 1C gene summary", fig)


def rnaseq_figures() -> None:
    de = table("RNAseq/native_deseq2_results_11_pdac.csv")
    de["padj"] = pd.to_numeric(de["padj"], errors="coerce").clip(lower=1e-300)
    de["log2FoldChange"] = pd.to_numeric(de["log2FoldChange"], errors="coerce")
    valid = de.dropna(subset=["padj", "log2FoldChange"])
    sig = valid[valid["padj"] < .05]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(valid["log2FoldChange"], -np.log10(valid["padj"]), s=8, color="#b8b8b8", alpha=.55, linewidths=0)
    ax.scatter(sig["log2FoldChange"], -np.log10(sig["padj"]), s=9, c=np.where(sig["log2FoldChange"] > 0, "#d1495b", "#0077b6"), alpha=.75, linewidths=0)
    ax.axhline(-np.log10(.05), color="#666666", linestyle="--", linewidth=.8)
    ax.axvline(0, color="#666666", linewidth=.7)
    ax.set(xlabel="DESeq2 log2 fold change (tumour / normal)", ylabel="-log10 adjusted p-value", title="Native paired DESeq2: 11 PDAC patients")
    export("Figure 2/Figure 2A volcano plot", fig)

    delta = table("RNAseq/tumour_normal_log2_delta_matrix_11_pdac.csv", index_col=0)
    delta.columns = delta.columns.astype(str)
    delta = delta[PATIENTS]
    top = delta.var(axis=1).nlargest(min(2000, len(delta))).index
    values = StandardScaler().fit_transform(delta.loc[top].T)
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(values)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(coords[:, 0], coords[:, 1], s=70, color="#4b3f8c", edgecolor="white", linewidth=.6)
    for i, patient in enumerate(PATIENTS):
        ax.annotate(patient, coords[i], xytext=(5, 5), textcoords="offset points", fontsize=8)
    ax.set(xlabel=f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)", ylabel=f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)", title="RNA-seq tumour-normal response PCA: revised 11-patient PDAC cohort")
    ax.grid(True, alpha=.25, linestyle=":")
    export("Figure 2/Figure 2B PCA", fig)

    fractions = table("RNAseq/quantiseq_cell_fractions_11_pairs.csv")
    cell_types = [col for col in fractions.columns if col != "Sample"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    colors = sns.color_palette("tab20", len(cell_types))
    for ax, condition, suffix in zip(axes, ["N", "T"], ["Normal", "Tumour"]):
        rows = fractions[fractions["Sample"].astype(str).str.endswith(condition)].copy()
        rows["patient"] = rows["Sample"].astype(str).str[:-1]
        rows = rows.set_index("patient").reindex(PATIENTS)
        bottom = np.zeros(len(PATIENTS))
        for cell, color in zip(cell_types, colors):
            values = rows[cell].to_numpy(dtype=float) * 100
            ax.bar(PATIENTS, values, bottom=bottom, label=cell, color=color, width=.82)
            bottom += values
        ax.set_title(suffix)
        ax.set_xlabel("Patient")
        ax.tick_params(axis="x", rotation=45)
    axes[0].set_ylabel("Estimated cell fraction (%)")
    axes[1].legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.suptitle("quanTIseq immune-cell composition: 11 paired PDAC samples")
    fig.tight_layout()
    export("Figure 2/Figure 2C immune cell composition", fig)

    hallmark = table("RNAseq/figureD_hallmark_ssgsea_row_zscores_11_pdac.csv", index_col=0)
    tme = table("RNAseq/figureD_tme_feature_row_zscores_11_pdac.csv", index_col=0)
    fig, axes = plt.subplots(1, 2, figsize=(16, 10), gridspec_kw={"width_ratios": [1.45, 1]})
    sns.heatmap(hallmark, cmap="vlag", center=0, ax=axes[0], cbar_kws={"label": "Row z-score"})
    sns.heatmap(tme, cmap="vlag", center=0, ax=axes[1], cbar_kws={"label": "Row z-score"})
    axes[0].set_title("Hallmark programmes")
    axes[1].set_title("Immune and stromal features")
    for ax in axes:
        ax.set_xlabel("Tumour sample")
        ax.tick_params(axis="x", rotation=45)
    fig.suptitle("Independent tumour programmes and microenvironment profiles")
    fig.tight_layout()
    export("Figure 3/Figure 3", fig)

    immune = table("RNAseq/immune_stromal_paired_feature_summary_11_pdac.csv")
    heat = np.vstack([np.fromstring(value.strip("[]"), sep=" ") for value in immune["standardized_patient_changes"]])
    fig = plt.figure(figsize=(17, max(17, .28 * len(immune))))
    grid = fig.add_gridspec(1, 2, width_ratios=[3.4, 1.5], wspace=.04)
    ax_heat = fig.add_subplot(grid[0, 0])
    ax_summary = fig.add_subplot(grid[0, 1], sharey=ax_heat)
    cmap = sns.color_palette("vlag", as_cmap=True)
    cmap.set_bad("#d9d9d9")
    image = ax_heat.imshow(np.ma.masked_invalid(heat), aspect="auto", interpolation="nearest", cmap=cmap, norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5))
    ax_heat.set_xticks(np.arange(len(PATIENTS)), PATIENTS)
    ax_heat.set_yticks(np.arange(len(immune)), immune["feature_label"], fontsize=7)
    ax_heat.set_xlabel("Patient")
    ax_heat.set_title("Tumour minus normal change, standardized within each feature")
    fig.colorbar(image, ax=ax_heat, fraction=.025, pad=.015).set_label("Within-feature SD units")
    values = immune["median_paired_change_sd"].to_numpy(dtype=float)
    nominal = immune["p_value_nominal"].to_numpy(dtype=float) < .05
    y = np.arange(len(immune))
    ax_summary.axvline(0, color="black", linewidth=.7)
    ax_summary.scatter(values[~nominal], y[~nominal], color="#666666", s=17, label="Nominal p >= 0.05")
    ax_summary.scatter(values[nominal], y[nominal], color="#d1495b", s=24, label="Nominal p < 0.05")
    ax_summary.set_xlabel("Median paired change (SD units)")
    ax_summary.set_title("Paired-change summary")
    ax_summary.tick_params(axis="y", left=False, labelleft=False)
    ax_summary.legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle("Paired immune and stromal estimates across 11 PDAC patients")
    export("Supplementary/Figure S1 immune stromal", fig)


def proteomics_volcano(relative_table: str, relative_output: str, title: str) -> None:
    de = table(relative_table)
    fc_col = "mean_log2FC" if "mean_log2FC" in de else "logFC"
    p_col = "p_value" if "p_value" in de else "P.Value"
    fdr_col = "FDR" if "FDR" in de else "adj.P.Val"
    fc = pd.to_numeric(de[fc_col], errors="coerce")
    fdr = pd.to_numeric(de[fdr_col], errors="coerce").clip(lower=1e-300)
    ok = fc.notna() & fdr.notna()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(fc[ok], -np.log10(fdr[ok]), s=18, c=np.where(fdr[ok] < .05, "#d1495b", "#b8b8b8"), alpha=.75, linewidths=0)
    ax.axhline(-np.log10(.05), color="#666666", linestyle="--", linewidth=.8)
    ax.axvline(0, color="#666666", linewidth=.7)
    ax.set(xlabel="Mean log2 fold change", ylabel="-log10 adjusted p-value", title=title)
    export(relative_output, fig)


def integration_figures() -> None:
    de = table("RNAseq/native_deseq2_results_11_pdac.csv")
    de_symbols = set(de.loc[(pd.to_numeric(de["padj"], errors="coerce") < .05) & de["gene_symbol"].notna(), "gene_symbol"].astype(str))
    matrix = table("WES/wes_mutation_matrix.blacklist_filtered.csv", index_col=0).fillna(0)
    mutated = set(matrix.index[matrix.sum(axis=1) > 0].astype(str))
    overlap = de_symbols & mutated
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.add_patch(plt.Circle((2.4, 2.5), 1.55, color="#5da5da", alpha=.55))
    ax.add_patch(plt.Circle((4.2, 2.5), 1.55, color="#f28e2b", alpha=.55))
    ax.text(1.7, 2.5, f"DE only\n{len(de_symbols - overlap)}", ha="center", va="center")
    ax.text(3.3, 2.5, f"Overlap\n{len(overlap)}", ha="center", va="center")
    ax.text(4.9, 2.5, f"Mutated only\n{len(mutated - overlap)}", ha="center", va="center")
    ax.set_xlim(.2, 6.2); ax.set_ylim(.4, 4.7); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("Differential-expression and mutation overlap")
    export("Figure 5/Figure 5A expression overlap", fig)

    factors = table("MOFA/mofa_factors.csv", index_col=0)
    factors.index = factors.index.astype(str)
    factors = factors[[f"factor_{i}" for i in range(1, 6)]]
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(factors.T, cmap="vlag", center=0, linewidths=.4, linecolor="white", cbar_kws={"label": "Factor score"}, ax=ax)
    ax.set(xlabel="Patient", ylabel="MOFA factor", title="MOFA factor scores across the revised PDAC cohort")
    export("Figure 5/Figure 5B heatmap", fig)

    clusters = table("MOFA/mofa_patient_clusters.csv", index_col=0)
    clusters = clusters.iloc[:, 0].reindex(factors.index)
    fig, ax = plt.subplots(figsize=(7, 6))
    sns.scatterplot(x=factors["factor_1"], y=factors["factor_2"], hue=clusters.astype(str), palette="Set2", s=95, edgecolor="black", ax=ax)
    for sample, row in factors.iterrows():
        ax.annotate(sample.rstrip("T"), (row["factor_1"], row["factor_2"]), xytext=(5, 4), textcoords="offset points", fontsize=8)
    ax.set(xlabel="Factor 1", ylabel="Factor 2", title="Exploratory MOFA patient clustering")
    ax.grid(True, alpha=.25, linestyle=":")
    export("Figure 5/Figure 5C cluster", fig)

    candidates = table("MOFA/mofa_top20_composite_score_11_pdac.csv").sort_values("composite_score")
    colors = candidates["dominant_modality"].map({"wes": "#2a9d8f", "rnaseq": "#4c78a8", "proteomics": "#e76f51"}).fillna("#999999")
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.barh(candidates["feature"], candidates["composite_score"], color=colors)
    ax.set(xlabel="Composite evidence score", ylabel="Feature", title="Top 20 cross-modal candidate features")
    export("Figure 5/Figure 5D candidates", fig)


def comparison_figures() -> None:
    coords = table("RNAseq/Figure_13_rna_pca_tcga_11_pdac_coordinates.csv", index_col=0)
    fig, ax = plt.subplots(figsize=(7, 6))
    for cohort, group in coords.groupby("cohort"):
        ax.scatter(group.PC1, group.PC2, s=48, alpha=.8, label=cohort)
    ax.set(xlabel="PC1", ylabel="PC2", title="Joint RNA-seq PCA: revised PDAC and TCGA PAAD")
    ax.legend(frameon=False)
    export("Figure 6/Figure 6A joint rna pca", fig)

    mutation = table("MOFA/revised_pdac_tcga_mutation_frequency.csv")
    mutation = mutation.sort_values("freq_diff").tail(30)
    fig, ax = plt.subplots(figsize=(8, 9))
    ax.barh(mutation["gene"], mutation["freq_diff"] * 100, color=np.where(mutation["freq_diff"] >= 0, "#e76f51", "#457b9d"))
    ax.axvline(0, color="black", linewidth=.7)
    ax.set(xlabel="PDAC minus TCGA PAAD mutation frequency (percentage points)", ylabel="Gene", title="Revised PDAC versus TCGA PAAD mutation frequency")
    export("Figure 6/Figure 6B pdac vs tcga mut freq", fig)

    theme_df = table("MOFA/revised_11_vs_tcga_pathway_theme_concordance.csv").sort_values(["african_terms", "tcga_terms"], ascending=False)
    themes = theme_df["theme"].tolist()
    palette = dict(zip(themes, sns.color_palette("tab10", n_colors=len(themes))))
    present = (theme_df["african_terms"] > 0) | (theme_df["tcga_terms"] > 0)
    both = (theme_df["african_terms"] > 0) & (theme_df["tcga_terms"] > 0)
    jaccard = both.sum() / present.sum() if present.sum() else 0.0
    rho, pvalue = spearmanr(theme_df["african_terms"], theme_df["tcga_terms"])
    fig, ax = plt.subplots(figsize=(6, 6))
    for _, row in theme_df.iterrows():
        ax.scatter(row["african_terms"], row["tcga_terms"], s=100, color=palette[row["theme"]], edgecolor="k", label=f"{row['theme']} ({int(row['african_terms'])}/{int(row['tcga_terms'])})")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, .5), fontsize=8)
    low = min(theme_df["african_terms"].min(), theme_df["tcga_terms"].min()) - 5
    high = max(theme_df["african_terms"].max(), theme_df["tcga_terms"].max()) + 5
    ax.plot([low, high], [low, high], "--", color="gray", linewidth=.8)
    ax.set_xlabel("African pathway term count (aggregated)")
    ax.set_ylabel("TCGA pathway term count (aggregated)")
    ax.set_title(f"Pathway theme concordance: revised 11-patient PDAC vs TCGA PAAD\nRecurrently mutated genes (>= 2 samples); Jaccard = {jaccard:.2f}; Spearman rho = {rho:.3f} (p = {pvalue:.3g})")
    fig.tight_layout()
    export("Figure 6/Figure 6C pathway theme concordance", fig)

    kras_patients = configured_cohort("all")
    kras_colors = {"G12V": "#4C78A8", "Q61L": "#F58518", "G12R": "#E45756", "G12D": "#54A24B", "Wild type": "#EEEEEE"}
    variants = table("Ampulla/kras_variants_all_cohorts.tsv", sep="\t")
    kras_status = {patient: "Wild type" for patient in kras_patients}
    for _, row in variants.iterrows():
        patient = str(row["patient"])
        if patient not in kras_status:
            continue
        position = int(row["POS"])
        protein = {"A": "G12V", "T": "G12D", "G": "G12R"}.get(str(row["ALT"])) if position in {25245350, 25245351} else None
        if position == 25227342:
            protein = "Q61L"
        kras_status[patient] = protein or "Mutated"
    fig, ax = plt.subplots(figsize=(12, 3.8))
    for index, patient in enumerate(kras_patients):
        category = kras_status[patient]
        ax.add_patch(Rectangle((index, 0), 1, 1, facecolor=kras_colors.get(category, "#984EA3"), edgecolor="white", linewidth=1.5))
        if category != "Wild type":
            ax.text(index + .5, .5, category, ha="center", va="center", fontsize=9, color="white", fontweight="bold")
    ax.axvline(len(PATIENTS), color="black", linewidth=2)
    ax.set_xlim(0, len(kras_patients)); ax.set_ylim(0, 1)
    ax.set_xticks([i + .5 for i in range(len(kras_patients))], [f"{p}T" for p in kras_patients])
    ax.set_yticks([.5], ["KRAS"])
    ax.tick_params(axis="both", length=0)
    ax.set_title("KRAS mutation status: PDAC versus ampulla\nExploratory patient-level comparison", pad=28)
    ax.text(len(PATIENTS) / 2, -.2, "PDAC (n=11)", transform=ax.get_xaxis_transform(), ha="center", va="top", fontweight="bold")
    ax.text(len(PATIENTS) + 1.5, -.2, "Ampulla (n=3)", transform=ax.get_xaxis_transform(), ha="center", va="top", fontweight="bold")
    ax.legend(handles=[Patch(facecolor=color, edgecolor="white", label=label) for label, color in kras_colors.items()], loc="upper center", bbox_to_anchor=(.5, -.34), ncol=5, frameon=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    export("Figure 7/Figure 7A kras mutation", fig)

    shared = table("Ampulla/pdac_ampulla_shared_pathways.csv")
    shared["minus_log10_fdr"] = pd.to_numeric(shared["minus_log10_fdr"], errors="coerce")
    shared = shared.sort_values("minus_log10_fdr").tail(20)
    labels = shared["Term"].astype(str).map(lambda value: textwrap.fill(value[:72], 42))
    colors = shared["modality"].map({"WES": "#264653", "RNA": "#e76f51", "Proteomics": "#2a9d8f"})
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.barh(labels, shared["minus_log10_fdr"], color=colors)
    ax.set(xlabel="-log10 adjusted p-value", ylabel="Shared pathway term", title="Shared pathway signals: PDAC versus ampulla")
    export("Figure 7/Figure 7B shared pathway", fig)


def mofa_supplements() -> None:
    variance = table("MOFA/mofa_native_variance_explained_11_pdac.csv")
    pivot = variance.pivot(index="factor", columns="view", values="r2_percent")
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="vlag", center=0, linewidths=.5, linecolor="white", cbar_kws={"label": "Native MOFA reconstruction R2 (%)"}, ax=ax)
    ax.set(title="Native MOFA variance explained by factor and view", xlabel="Data view", ylabel="Factor")
    export("Supplementary/Figure S9 mofa native variance", fig)

    factors = table("MOFA/mofa_factors.csv", index_col=0)
    factors.index = factors.index.astype(str)
    base = table("MOFA/mofa_patient_clusters.csv", index_col=0).iloc[:, 0].reindex(factors.index)
    stability = table("MOFA/mofa_cluster_sensitivity_summary.csv")
    scenarios = private_value("mofa_exclusion_scenarios")
    assignments = pd.DataFrame(index=factors.index, columns=list(scenarios))
    for scenario, dropped in scenarios.items():
        subset = factors.drop(index=dropped, errors="ignore")
        labels = pd.Series(KMeans(n_clusters=2, random_state=42, n_init=20).fit_predict(subset), index=subset.index)
        mapping = {cluster: base.loc[base.index.intersection(labels[labels == cluster].index)].mode().iat[0] for cluster in labels.unique()}
        aligned = labels.map(mapping)
        assignments.loc[aligned.index, scenario] = aligned.values
    plot_data = assignments.reindex(base.sort_values().index).astype(float)
    fig, (ax_heat, ax_sil) = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={"height_ratios": [3, 1]})
    im = ax_heat.imshow(plot_data.values, aspect="auto", cmap=ListedColormap(["#66c2a5", "#fc8d62"]), vmin=0, vmax=1)
    ax_heat.set_xticks(range(len(scenarios)), list(scenarios), rotation=20, ha="right")
    ax_heat.set_yticks(range(len(plot_data)), plot_data.index)
    ax_heat.set_ylabel("Patient")
    for row in range(plot_data.shape[0]):
        for col in range(plot_data.shape[1]):
            value = plot_data.values[row, col]
            if np.isfinite(value):
                ax_heat.text(col, row, int(value), ha="center", va="center", fontsize=8)
    ax_heat.set_title("Per-patient cluster assignment across exclusion scenarios")
    fig.colorbar(im, ax=ax_heat, ticks=[0, 1], fraction=.046, pad=.04).ax.set_yticklabels(["Cluster 0", "Cluster 1"])
    ax_sil.plot(stability["scenario"], stability["silhouette"], marker="o", color="#4b3f8c")
    for _, row in stability.iterrows():
        ax_sil.annotate(f"n={int(row['n_patients'])}", (row["scenario"], row["silhouette"]), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=8)
    ax_sil.set_ylabel("Silhouette score")
    ax_sil.tick_params(axis="x", rotation=20)
    ax_sil.grid(axis="y", alpha=.25, linestyle=":")
    export("Supplementary/Figure S10 mofa cluster robustness", fig)

    support = table("MOFA/mofa_modality_support_top200_11_pdac.csv").sort_values("count")
    fig, ax = plt.subplots(figsize=(9, 6))
    bars = ax.barh(support["modality_support"], support["count"], color="#4c78a8")
    for bar, value in zip(bars, support["count"]):
        ax.annotate(str(int(value)), (value, bar.get_y() + bar.get_height()/2), xytext=(4, 0), textcoords="offset points", va="center", fontsize=8)
    ax.set(xlabel="Number of candidates", ylabel="Supporting modalities", title="Modality support among top 200 MOFA candidates")
    export("Supplementary/Figure S11 mofa modality", fig)

    heat = table("MOFA/kras_integrated_feature_heatmap_11_pdac.csv", index_col=0)
    fig, ax = plt.subplots(figsize=(11, 9))
    sns.heatmap(heat, cmap="vlag", center=0, linewidths=.25, linecolor="white", ax=ax)
    ax.set(xlabel="Patient", ylabel="Integrated feature", title="Integrated mutation, RNA and protein features")
    export("Supplementary/Figure S12 integrated heatmap", fig)


def ancestry_supplements() -> None:
    summary = table("MOFA/revised_tcga_ancestry_modality_summary.csv")
    summary = summary[summary["quantity"].isin(["RNA", "RPPA", "MC3"])].copy()
    x = np.arange(len(summary)); width = .36
    fig, ax = plt.subplots(figsize=(9, 5.5))
    full = pd.to_numeric(summary["full_tcga"], errors="coerce")
    subset = pd.to_numeric(summary["african_ancestry_subset"], errors="coerce")
    ax.bar(x - width/2, full, width, label="TCGA full", color="#457b9d")
    ax.bar(x + width/2, subset, width, label="AFR/AFR_ADMIX subset", color="#e76f51")
    ax.set_xticks(x, summary["quantity"])
    ax.set(ylabel="Number of patients", title="TCGA PAAD modality availability by genetic ancestry")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=.25, linestyle=":")
    export("Supplementary/Figure S13 tcga ancestry counts", fig)

    overlap = table("MOFA/revised_tcga_ancestry_all_three_patients.csv")
    values = [11, int(subset.iloc[0]), int(subset.iloc[1]), int(subset.iloc[2]), len(overlap)]
    labels = ["Revised PDAC", "TCGA AFR RNA", "TCGA AFR RPPA", "TCGA AFR MC3", "TCGA AFR all three"]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.barh(labels, values, color=["#2a9d8f", "#e76f51", "#e76f51", "#e76f51", "#8c564b"])
    for y, value in enumerate(values):
        ax.annotate(str(value), (value, y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=8)
    ax.set(xlabel="Number of patients", title="Effective sample sizes for ancestry-restricted comparisons")
    ax.grid(axis="x", alpha=.25, linestyle=":")
    export("Supplementary/Figure S14 ancestry overlap", fig)


def normalization_supplement() -> None:
    original = table("Proteomics/proteomics_log2.csv", index_col=0)
    normalized = table("Proteomics/proteomics_log2_zscore.csv", index_col=0)
    sample_order = [sample for sample in SAMPLES if sample in original.columns and sample in normalized.columns]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    for ax, matrix, title in zip(axes, [original, normalized], ["Original responses", "Normalized responses"]):
        values = [pd.to_numeric(matrix[sample], errors="coerce").dropna().to_numpy() for sample in sample_order]
        ax.boxplot(values, tick_labels=sample_order, showfliers=False, patch_artist=True, boxprops={"facecolor": "#9ecae1", "edgecolor": "#356a8a"}, medianprops={"color": "#b23a48", "linewidth": 1.2})
        ax.set_title(title)
        ax.set_xlabel("Patient sample")
        ax.tick_params(axis="x", rotation=45)
        ax.grid(axis="y", alpha=.2, linestyle=":")
    axes[0].set_ylabel("Proteomics feature value")
    fig.suptitle("Proteomics responses before and after row-wise z-score normalization")
    fig.tight_layout()
    export("Supplementary/Figure S2 normalisation", fig)


def write_manifest() -> None:
    native = [
        ("Figure 1/Figure 1A_mutation burden per patient", "native vector from study table", "WES/wes_mutation_burden.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 1/Figure 1B oncoprint", "native vector from study table", "WES/wes_mutation_matrix.blacklist_filtered.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 1/Figure 1C gene summary", "native vector from study table", "WES/gene_summary_table.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 2/Figure 2A volcano plot", "native vector from study table", "RNAseq/native_deseq2_results_11_pdac.csv", "source_code/run_native_deseq2_revised_11.R"),
        ("Figure 2/Figure 2B PCA", "native vector from study table", "RNAseq/tumour_normal_log2_delta_matrix_11_pdac.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 2/Figure 2C immune cell composition", "native vector from study table", "RNAseq/quantiseq_cell_fractions_11_pairs.csv", "source_code/run_salmon_quantiseq_11_pdac.R"),
        ("Figure 3/Figure 3", "native vector from study table", "RNAseq/figureD_hallmark_ssgsea_row_zscores_11_pdac.csv; RNAseq/figureD_tme_feature_row_zscores_11_pdac.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 5/Figure 5A expression overlap", "native vector from study tables", "RNAseq/native_deseq2_results_11_pdac.csv; WES/wes_mutation_matrix.blacklist_filtered.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 5/Figure 5B heatmap", "native vector from study table", "MOFA/mofa_factors.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 5/Figure 5C cluster", "native vector from study tables", "MOFA/mofa_factors.csv; MOFA/mofa_patient_clusters.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 5/Figure 5D candidates", "native vector from study table", "MOFA/mofa_top20_composite_score_11_pdac.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 6/Figure 6A joint rna pca", "native vector from study coordinates", "RNAseq/Figure_13_rna_pca_tcga_11_pdac_coordinates.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 6/Figure 6B pdac vs tcga mut freq", "native vector from study table", "MOFA/revised_pdac_tcga_mutation_frequency.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 6/Figure 6C pathway theme concordance", "native vector from study table", "MOFA/revised_11_vs_tcga_pathway_theme_concordance.csv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 7/Figure 7A kras mutation", "native vector from study table", "Ampulla/kras_variants_all_cohorts.tsv", "scripts/figures/generate_study_vectors.py"),
        ("Figure 7/Figure 7B shared pathway", "native vector from study table", "Ampulla/pdac_ampulla_shared_pathways.csv", "scripts/figures/generate_study_vectors.py"),
        ("Supplementary/Figure S1 immune stromal", "native vector from study table", "RNAseq/immune_stromal_paired_feature_summary_11_pdac.csv", "scripts/figures/generate_study_vectors.py"),
        ("Supplementary/Figure S9 mofa native variance", "native vector from study table", "MOFA/mofa_native_variance_explained_11_pdac.csv", "scripts/figures/generate_study_vectors.py"),
        ("Supplementary/Figure S10 mofa cluster robustness", "native vector from study tables", "MOFA/mofa_factors.csv; MOFA/mofa_patient_clusters.csv; MOFA/mofa_cluster_sensitivity_summary.csv", "scripts/figures/generate_study_vectors.py"),
        ("Supplementary/Figure S11 mofa modality", "native vector from study table", "MOFA/mofa_modality_support_top200_11_pdac.csv", "scripts/figures/generate_study_vectors.py"),
        ("Supplementary/Figure S12 integrated heatmap", "native vector from study table", "MOFA/kras_integrated_feature_heatmap_11_pdac.csv", "scripts/figures/generate_study_vectors.py"),
        ("Supplementary/Figure S13 tcga ancestry counts", "native vector from study table", "MOFA/revised_tcga_ancestry_modality_summary.csv", "scripts/figures/generate_study_vectors.py"),
        ("Supplementary/Figure S14 ancestry overlap", "native vector from study table", "MOFA/revised_tcga_ancestry_modality_summary.csv; MOFA/revised_tcga_ancestry_all_three_patients.csv", "scripts/figures/generate_study_vectors.py"),
        ("Supplementary/Figure S2 normalisation", "native vector from raw and normalized study matrices", "Proteomics/proteomics_log2.csv; Proteomics/proteomics_log2_zscore.csv", "scripts/figures/generate_study_vectors.py"),
    ]
    records = []
    for stem, status, inputs, script in native:
        records.append({"figure": stem, "svg": f"{stem}.svg", "pdf": f"{stem}.pdf", "status": status, "raw_inputs": inputs, "generating_code": script})
    for stem, note in [
        ("Supplementary/Figure S3 CV", "Original numeric chart source was not located; SVG remains a raster-backed conversion."),
        ("Supplementary/Figure S4 peptide distribution", "Original numeric chart source was not located; SVG remains a raster-backed conversion."),
        ("Supplementary/Figure S5 protein distribution", "Original numeric chart source was not located; SVG remains a raster-backed conversion."),
        ("Supplementary/Figure S6 hydrophobicity index", "Original numeric chart source was not located; SVG remains a raster-backed conversion."),
        ("Supplementary/Figure S7 isoelectric point", "Original numeric chart source was not located; SVG remains a raster-backed conversion."),
        ("Figure 4/Figure 4 volcano plot", "The available 175-feature exploratory table has uniform FDR 0.160966 and does not reproduce the source PDF's 1,179 protein groups; retain raster-backed original."),
        ("Supplementary/Figure S8 volcano plot 19", "The 1,179-feature source analysis table was not located; retain raster-backed original."),
    ]:
        pdf = f"raw_inputs/{stem}.pdf"
        report = "raw_inputs/study_data/Proteomics/20260116_123225_20241023_FFPE TISSUE_PDAC_SNE_Report.tsv"
        records.append({"figure": stem, "svg": f"{stem}.svg", "pdf": f"{stem}.pdf", "status": note, "raw_inputs": f"{pdf}; {report}", "generating_code": "scripts/figures/convert_figures.py (conversion only; original analysis code unavailable)"})
    pd.DataFrame(records).to_csv(ROOT / "figure_source_manifest.csv", index=False)


def main() -> None:
    wes_figures()
    rnaseq_figures()
    integration_figures()
    comparison_figures()
    mofa_supplements()
    ancestry_supplements()
    normalization_supplement()
    write_manifest()
    print("Generated native SVG/PDF figures from archived study analysis tables")


if __name__ == "__main__":
    main()