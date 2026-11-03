#!/usr/bin/env python3
# Current collaborator code export. Final supplied vector presentation scripts, consuming private derived tables.
# Source: PDAC_Manuscript_Figs_Sne/Immune_rerun_NN/generate_vectors.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Regenerate immune/TME figures as SVG and PDF from archived study tables."""

from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.gridspec import GridSpec


ROOT = figure_root() / "immune"
TABLES = private_path("immune_tables")
PATIENTS = configured_cohort("primary_pdac")
SAMPLES = [f"{patient}T" for patient in PATIENTS]
METHOD_ORDER = ["ESTIMATE", "MCP-counter", "EPIC", "xCell", "quanTIseq"]
METHOD_FILES = {
    "ESTIMATE": ("estimate_native_scores_11_pairs.csv", "estimate_paired_comparison_11_pdac.csv", None),
    "MCP-counter": ("mcp_counter_native_scores_11_pairs.csv", "mcp_counter_paired_comparison_11_pdac.csv", None),
    "EPIC": ("epic_native_scores_11_pairs.csv", "epic_converged_pairs_sensitivity.csv", configured_cohort("epic_converged")),
    "xCell": ("xcell_native_scores_11_pairs.csv", "xcell_paired_comparison_11_pdac.csv", None),
    "quanTIseq": ("quantiseq_cell_fractions_11_pairs.csv", "quantiseq_paired_comparison_11_pdac.csv", None),
}
plt.rcParams["svg.fonttype"] = "none"


def load(name: str, **kwargs) -> pd.DataFrame:
    return pd.read_csv(TABLES / name, **kwargs)


def save(fig: plt.Figure, relative_png: str) -> None:
    path = ROOT / relative_png
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".svg"), format="svg", bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
    plt.close(fig)


def draw_heatmap(ax, frame: pd.DataFrame, vmin: float = -2, vmax: float = 2) -> None:
    values = frame.to_numpy(dtype=float)
    cmap = sns.color_palette("vlag", as_cmap=True)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            value = values[row, column]
            color = "#d9d9d9" if not np.isfinite(value) else cmap(norm(value))
            ax.add_patch(plt.Rectangle((column, row), 1, 1, facecolor=color, edgecolor="white", linewidth=.25))
    ax.set_xlim(0, values.shape[1]); ax.set_ylim(values.shape[0], 0)
    ax.set_xticks(np.arange(values.shape[1]) + .5, frame.columns.astype(str))
    ax.set_yticks(np.arange(values.shape[0]) + .5, frame.index.astype(str))
    ax.tick_params(length=0)


def draw_colorbar(ax, vmin: float = -2, vmax: float = 2, label: str = "Within-feature SD units") -> None:
    cmap = sns.color_palette("vlag", as_cmap=True)
    steps = 64
    for index in range(steps):
        value = vmin + (index + .5) * (vmax - vmin) / steps
        ax.add_patch(plt.Rectangle((0, index), 1, 1, facecolor=cmap(Normalize(vmin, vmax)(value)), edgecolor="none"))
    ax.set_xlim(0, 1); ax.set_ylim(0, steps)
    ax.set_xticks([])
    ticks = np.linspace(0, steps, 5)
    ax.set_yticks(ticks, [f"{value:g}" for value in np.linspace(vmin, vmax, 5)])
    ax.set_ylabel(label, rotation=90, labelpad=8)
    ax.tick_params(axis="y", labelsize=7, length=2, pad=2)


def immune_stromal() -> None:
    data = load("immune_stromal_paired_feature_summary_11_pdac.csv")
    heat = np.vstack([np.fromstring(str(value).strip("[]"), sep=" ") for value in data["standardized_patient_changes"]])
    fig = plt.figure(figsize=(17, max(17, .28 * len(data))))
    grid = GridSpec(1, 2, figure=fig, width_ratios=[3.4, 1.5], wspace=.04)
    ax_heat = fig.add_subplot(grid[0, 0])
    ax_summary = fig.add_subplot(grid[0, 1], sharey=ax_heat)
    cmap = sns.color_palette("vlag", as_cmap=True)
    cmap.set_bad("#d9d9d9")
    heat_frame = pd.DataFrame(heat, columns=PATIENTS, index=data["feature_label"])
    draw_heatmap(ax_heat, heat_frame, vmin=-2.5, vmax=2.5)
    ax_heat.tick_params(axis="y", labelsize=7)
    ax_heat.set_xlabel("Patient")
    ax_heat.set_title("Patient-level tumour - normal change\n(row-standardized within each method-feature)")
    draw_colorbar(ax_heat.inset_axes([1.01, .08, .025, .84]), vmin=-2.5, vmax=2.5)
    values = data["median_paired_change_sd"].to_numpy(dtype=float)
    nominal = data["p_value_nominal"].to_numpy(dtype=float) < .05
    y = np.arange(len(data))
    ax_summary.axvline(0, color="black", linewidth=.7)
    ax_summary.scatter(values[~nominal], y[~nominal], color="#666666", s=17, label="Nominal p >= 0.05")
    ax_summary.scatter(values[nominal], y[nominal], color="#d1495b", s=24, label="Nominal p < 0.05")
    ax_summary.set_xlabel("Median paired change\n(within-feature SD units)")
    ax_summary.set_title("Paired-change summary\n(nominal Wilcoxon p-values)")
    ax_summary.tick_params(axis="y", left=False, labelleft=False)
    ax_summary.legend(loc="lower right", fontsize=8, frameon=False)
    offsets = np.cumsum([sum(data.method.eq(method)) for method in METHOD_ORDER])
    for boundary in offsets[:-1]:
        ax_heat.axhline(boundary, color="black", linewidth=1)
        ax_summary.axhline(boundary - .5, color="black", linewidth=1)
    epi_rows = data.method.eq("EPIC").to_numpy()
    labels = [label.get_text() + (" dagger" if epi_rows[i] else "") for i, label in enumerate(ax_heat.get_yticklabels())]
    ax_heat.set_yticklabels(labels, fontsize=7)
    fig.suptitle("Paired immune and stromal comparisons: revised 11-patient PDAC cohort", y=.995, fontsize=14, fontweight="bold")
    fig.text(.01, .002, "EPIC uses eight converged matched pairs; other methods use 11 pairs. Red points denote nominal p<0.05; no feature survived method-wise BH FDR<0.05. Standardized changes aid display across native method scales.", fontsize=8, va="bottom")
    fig.subplots_adjust(left=.31, right=.98, top=.96, bottom=.04)
    save(fig, "Figure_C_immune_stromal_paired_11_pdac.png")


def independent_programmes() -> None:
    programmes = load("figureD_hallmark_ssgsea_row_zscores_11_pdac.csv", index_col=0)
    tme = load("figureD_tme_feature_row_zscores_11_pdac.csv", index_col=0)
    composites = load("figureD_tme_composite_scores_11_pdac.csv", index_col=0)
    cutoffs = load("figureD_tme_median_cutoffs_11_pdac.csv").set_index("score")["median_cutoff"]
    fig = plt.figure(figsize=(15, 22), constrained_layout=False)
    grid = fig.add_gridspec(3, 2, width_ratios=[40, 1], height_ratios=[19, 12, 10], hspace=.34, wspace=.08)
    ax_programmes = fig.add_subplot(grid[0, 0]); cax_programmes = fig.add_subplot(grid[0, 1])
    ax_tme = fig.add_subplot(grid[1, 0]); cax_tme = fig.add_subplot(grid[1, 1])
    ax_composite = fig.add_subplot(grid[2, 0]); cax_composite = fig.add_subplot(grid[2, 1])
    cmap = sns.color_palette("vlag", as_cmap=True)
    draw_heatmap(ax_programmes, programmes)
    draw_colorbar(cax_programmes)
    ax_programmes.set_title("A  Hallmark programme scores (ssGSEA)", loc="left", fontsize=12, fontweight="bold", pad=10)
    ax_programmes.set_ylabel(""); ax_programmes.set_xlabel(""); ax_programmes.tick_params(axis="y", labelsize=9, length=0); ax_programmes.tick_params(axis="x", bottom=False, labelbottom=False)
    draw_heatmap(ax_tme, tme)
    draw_colorbar(cax_tme)
    ax_tme.set_title("B  Immune and stromal estimates", loc="left", fontsize=12, fontweight="bold", pad=10)
    ax_tme.set_ylabel(""); ax_tme.set_xlabel("PDAC tumour sample"); ax_tme.tick_params(axis="y", labelsize=9, length=0); ax_tme.tick_params(axis="x", labelrotation=45, labelsize=9); ax_tme.axhline(9, color="black", linewidth=1.2)
    cax_composite.set_axis_off()
    colors = {"Stromal/EMT-high, immune-low": "#d1495b", "Immune-high, stromal-low": "#00798c", "Mixed/intermediate": "#edae49"}
    groups = load("figureD_tme_exploratory_groups_11_pdac.csv").set_index("sample")["exploratory_group"]
    for group, color in colors.items():
        samples = [sample for sample in SAMPLES if groups.get(sample) == group]
        if samples:
            ax_composite.scatter(composites.loc["Immune composite", samples], composites.loc["Stromal composite", samples], color=color, s=65, label=group, zorder=3)
    ax_composite.axvline(cutoffs["Immune composite"], color="#444444", linestyle="--", linewidth=1)
    ax_composite.axhline(cutoffs["Stromal composite"], color="#444444", linestyle="--", linewidth=1)
    for sample in SAMPLES:
        ax_composite.annotate(sample[:-1], (composites.loc["Immune composite", sample], composites.loc["Stromal composite", sample]), xytext=(4, 3), textcoords="offset points", fontsize=8)
    ax_composite.set_title("C  Exploratory cohort-relative TME groups", loc="left", fontsize=12, fontweight="bold", pad=10)
    ax_composite.set_xlabel("Immune composite (mean feature z-score)"); ax_composite.set_ylabel("Stromal composite (mean feature z-score)")
    ax_composite.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8, frameon=False); ax_composite.grid(color="#dddddd", linewidth=.5, zorder=0)
    fig.suptitle("Tumour programme and microenvironment profiles in 11 PDAC cases", fontsize=15, fontweight="bold", y=.99)
    fig.text(.01, .005, "Independent exploratory replacement. Hallmark and immune/stromal heatmap cells are row z-scores across observed tumours, not inferential p-values. Cohort-relative groups are descriptive, not validated subtypes.", fontsize=8, va="bottom", wrap=True)
    fig.subplots_adjust(left=.30, right=.82, top=.965, bottom=.055)
    save(fig, "Figure_D_independent_programmes_TME_11_pdac.png")

    fig, ax = plt.subplots(figsize=(12, 9))
    draw_heatmap(ax, programmes)
    draw_colorbar(ax.inset_axes([1.01, .1, .025, .8]))
    ax.set_title("A  Hallmark programme scores (ssGSEA)")
    ax.set_xlabel("PDAC tumour sample"); ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=45, labelsize=9); ax.tick_params(axis="y", labelsize=9)
    save(fig, "Figure_DA_hallmark_programmes.png")

    fig, ax = plt.subplots(figsize=(12, 7))
    draw_heatmap(ax, tme)
    draw_colorbar(ax.inset_axes([1.01, .1, .025, .8]))
    ax.set_title("B  Immune and stromal estimates")
    ax.set_xlabel("PDAC tumour sample"); ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=45, labelsize=9); ax.tick_params(axis="y", labelsize=9)
    ax.axhline(9, color="black", linewidth=1.2)
    save(fig, "Figure_DB_immune_stromal_estimates.png")

    fig, ax = plt.subplots(figsize=(9, 7))
    for group, color in colors.items():
        samples = [sample for sample in SAMPLES if groups.get(sample) == group]
        if samples:
            ax.scatter(composites.loc["Immune composite", samples], composites.loc["Stromal composite", samples], color=color, s=70, label=group, zorder=3)
    ax.axvline(cutoffs["Immune composite"], color="#444444", linestyle="--", linewidth=1)
    ax.axhline(cutoffs["Stromal composite"], color="#444444", linestyle="--", linewidth=1)
    for sample in SAMPLES:
        ax.annotate(sample[:-1], (composites.loc["Immune composite", sample], composites.loc["Stromal composite", sample]), xytext=(5, 4), textcoords="offset points", fontsize=9)
    ax.set_title("C  Exploratory cohort-relative TME groups")
    ax.set_xlabel("Immune composite (mean feature z-score)"); ax.set_ylabel("Stromal composite (mean feature z-score)")
    ax.legend(loc="best", frameon=False, fontsize=8); ax.grid(color="#dddddd", linewidth=.5, zorder=0)
    save(fig, "Figure_DC_exploratory_TME_groups.png")


def quantiseq_composition() -> None:
    data = load("quantiseq_cell_fractions_11_pairs.csv")
    cell_types = [column for column in data.columns if column != "Sample"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    colors = sns.color_palette("tab20", len(cell_types))
    for ax, condition, label in zip(axes, ["N", "T"], ["Normal", "Tumour"]):
        block = data[data["Sample"].astype(str).str.endswith(condition)].copy()
        block["patient"] = block["Sample"].astype(str).str[:-1]
        block = block.set_index("patient").reindex(PATIENTS)
        bottom = np.zeros(len(PATIENTS))
        for cell_type, color in zip(cell_types, colors):
            values = pd.to_numeric(block[cell_type], errors="coerce").to_numpy() * 100
            ax.bar(PATIENTS, values, bottom=bottom, label=cell_type, color=color, width=.8)
            bottom += values
        ax.set_title(label); ax.set_xlabel("Patient"); ax.tick_params(axis="x", rotation=45)
    axes[0].set_ylabel("Estimated cell fraction (%)")
    axes[1].legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.suptitle("quanTIseq immune-cell composition: 11 paired PDAC samples")
    fig.tight_layout()
    save(fig, "Figure_immune_composition_quantiseq_11_pdac.png")


def quantiseq_changes() -> None:
    data = load("quantiseq_cell_fractions_11_pairs.csv")
    cell_types = [column for column in data.columns if column != "Sample"]
    fig, axes = plt.subplots(int(np.ceil(len(cell_types) / 4)), 4, figsize=(12, 8), squeeze=False)
    palette = {"Normal": "#4477AA", "Tumour": "#EE6677"}
    for index, cell_type in enumerate(cell_types):
        ax = axes.flat[index]
        for patient in PATIENTS:
            rows = data.set_index("Sample")
            values = [rows.loc[f"{patient}{condition}", cell_type] * 100 for condition in ["N", "T"]]
            ax.plot([0, 1], values, color="#999999", alpha=.65, linewidth=.8)
            ax.scatter([0, 1], values, c=[palette["Normal"], palette["Tumour"]], s=15, zorder=3)
        ax.set_title(cell_type.replace(".", " "), fontsize=8)
        ax.set_xticks([0, 1], ["Normal", "Tumour"])
        ax.set_ylabel("Estimated cell fraction (%)", fontsize=7)
        ax.grid(axis="y", alpha=.2, linestyle=":")
    for ax in axes.flat[len(cell_types):]: ax.axis("off")
    fig.suptitle("Paired tumour-normal quanTIseq estimates: 11 PDAC patients")
    fig.tight_layout()
    save(fig, "Additional /Figure_immune_paired_changes_quantiseq_11_pdac.png")


def quantiseq_table() -> None:
    paired = load("quantiseq_paired_comparison_11_pdac.csv")
    names = {"B.cells": "B cells", "Macrophages.M1": "Macrophages M1", "Macrophages.M2": "Macrophages M2", "Monocytes": "Monocytes", "Neutrophils": "Neutrophils", "NK.cells": "NK cells", "T.cells.CD4": "CD4 T cells", "T.cells.CD8": "CD8 T cells", "Tregs": "Regulatory T cells", "Dendritic.cells": "Dendritic cells", "Other": "Other (residual)"}
    display = pd.DataFrame({"Cell type": paired.cell_type.map(names), "Normal mean (%)": paired.normal_mean_percent.map(lambda x: f"{x:.2f}"), "Tumour mean (%)": paired.tumour_mean_percent.map(lambda x: f"{x:.2f}"), "Median paired delta (percentage points)": paired.median_delta_percentage_points.map(lambda x: f"{x:+.2f}"), "Paired p-value": paired.p_value.map(lambda x: f"{x:.3f}"), "BH FDR": paired.FDR.map(lambda x: f"{x:.3f}")})
    fig, ax = plt.subplots(figsize=(13, 4.5)); ax.axis("off")
    tab = ax.table(cellText=display.values, colLabels=display.columns, colWidths=[.23, .14, .14, .23, .14, .12], cellLoc="left", colLoc="left", bbox=[0, 0, 1, 1])
    tab.auto_set_font_size(False); tab.set_fontsize(10)
    for (row, _), cell in tab.get_celld().items():
        cell.PAD = .04; cell.set_edgecolor("#D0D0D0"); cell.set_linewidth(.4); cell.set_facecolor("#DDEBDD" if row == 0 else ("#F3F3F3" if row % 2 else "white"))
        if row == 0: cell.set_text_props(weight="bold"); cell.set_height(cell.get_height() * 1.6)
    save(fig, "Table_immune_paired_quantiseq_11_pdac.png")


def method_panels() -> None:
    for method, (score_file, test_file, eligible) in METHOD_FILES.items():
        if method == "quanTIseq":
            continue
        scores = load(score_file, index_col=0)
        scores.columns = scores.columns.astype(str)
        tests = load(test_file)
        feature_column = "cell_type" if "cell_type" in tests.columns else "feature"
        tests = tests.set_index(feature_column)
        features = list(scores.index)
        if method == "xCell":
            features = list(tests.sort_values(["FDR", "p_value"]).head(12).index.intersection(scores.index))
        if eligible is None:
            eligible = PATIENTS
        rows = [feature for feature in features]
        ncols = 4; nrows = int(np.ceil(len(rows) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(12, max(4.5, 2.2 * nrows)), squeeze=False)
        for index, feature in enumerate(rows):
            ax = axes.flat[index]
            for patient in eligible:
                normal, tumour = f"{patient}N", f"{patient}T"
                if normal not in scores.columns or tumour not in scores.columns: continue
                values = pd.to_numeric(scores.loc[feature, [normal, tumour]], errors="coerce")
                if values.notna().all():
                    ax.plot([0, 1], values.to_numpy(), color="#AAAAAA", linewidth=.8)
                    ax.scatter([0, 1], values.to_numpy(), c=["#4477AA", "#EE6677"], s=13, zorder=3)
            ax.set_title(str(feature), fontsize=8); ax.set_xticks([0, 1], ["Normal", "Tumour"]); ax.tick_params(axis="x", labelsize=7); ax.set_ylabel("Native score", fontsize=7); ax.grid(axis="y", alpha=.2, linestyle=":")
        for ax in axes.flat[len(rows):]: ax.axis("off")
        subtitle = "12 lowest-FDR features shown; all scores/tests in source tables" if method == "xCell" else ("Provisional: EPIC estimates for nonconverged samples shown" if method == "EPIC" else None)
        fig.suptitle(f"{method} paired tumour-normal inference: 11 PDAC pairs" + (f"\n{subtitle}" if subtitle else ""))
        fig.tight_layout()
        save(fig, f"Additional /Figure_TME_paired_{method.lower().replace('-', '_')}_11_pdac.png")


def main() -> None:
    immune_stromal()
    independent_programmes()
    quantiseq_composition()
    quantiseq_changes()
    quantiseq_table()
    method_panels()
    print("Exported all immune rerun panels as SVG and PDF")


if __name__ == "__main__":
    main()
