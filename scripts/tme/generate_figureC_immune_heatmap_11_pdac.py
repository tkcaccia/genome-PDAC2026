#!/usr/bin/env python3
# PDAC2026 analysis: Expression-to-deconvolution scores, paired comparisons and tumour-only Hallmark/TME groups.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Generate the 11-patient analogue of the paired immune/TME Figure C."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from scipy.stats import wilcoxon

ROOT = project_root()
RUN = ROOT / "results/immune_salmon_11_20261012"
OUT = ROOT / "results/figures/final_11_figures/Figure_C_immune_stromal_paired_11_pdac.png"
CSV = RUN / "tables/immune_stromal_paired_feature_summary_11_pdac.csv"
PATIENTS = configured_cohort("primary_pdac")
METHOD_ORDER = ["ESTIMATE", "MCP-counter", "EPIC", "xCell", "quanTIseq"]

FILES = {
    "ESTIMATE": ("estimate_native_scores_11_pairs.csv", "estimate_paired_comparison_11_pdac.csv", None),
    "MCP-counter": ("mcp_counter_native_scores_11_pairs.csv", "mcp_counter_paired_comparison_11_pdac.csv", None),
    "EPIC": ("epic_native_scores_11_pairs.csv", "epic_converged_pairs_sensitivity.csv", configured_cohort("epic_converged")),
    "xCell": ("xcell_native_scores_11_pairs.csv", "xcell_paired_comparison_11_pdac.csv", None),
    "quanTIseq": ("quantiseq_cell_fractions_11_pairs.csv", "quantiseq_paired_comparison_11_pdac.csv", None),
}


def load_method(method, score_file, test_file, eligible_patients):
    scores = pd.read_csv(RUN / "tables" / score_file)
    if method == "quanTIseq":
        scores = scores.set_index("Sample").T
    else:
        scores = scores.set_index(scores.columns[0])
    scores.columns = scores.columns.astype(str)
    tests = pd.read_csv(RUN / "tables" / test_file)
    feature_col = "cell_type" if "cell_type" in tests.columns else "feature"
    tests = tests.set_index(feature_col)
    if eligible_patients is None:
        eligible_patients = PATIENTS
    rows = []
    for feature in scores.index:
        deltas = []
        used_patients = []
        for patient in eligible_patients:
            normal, tumour = f"{patient}N", f"{patient}T"
            if normal not in scores.columns or tumour not in scores.columns:
                continue
            values = pd.to_numeric(scores.loc[feature, [normal, tumour]], errors="coerce")
            if values.notna().all():
                deltas.append(float(values[tumour] - values[normal]))
                used_patients.append(patient)
        if feature not in tests.index:
            continue
        test_row = tests.loc[feature]
        if isinstance(test_row, pd.DataFrame):
            test_row = test_row.iloc[0]
        delta_sd = float(np.std(deltas, ddof=1)) if len(deltas) > 1 else 0.0
        scaled = np.array(deltas) / delta_sd if delta_sd > 0 else np.zeros(len(deltas))
        patient_scaled = np.full(len(PATIENTS), np.nan)
        for patient, value in zip(used_patients, scaled):
            patient_scaled[PATIENTS.index(patient)] = value
        rows.append({
            "method": method,
            "feature": str(feature),
            "n_pairs": len(deltas),
            "patients_used": ";".join(used_patients),
            "median_paired_change_native": float(np.median(deltas)) if deltas else np.nan,
            "median_paired_change_sd": float(np.median(scaled)) if len(scaled) else np.nan,
            "p_value_nominal": float(test_row["p_value"]),
            "FDR_within_method": float(test_row["FDR"]),
            "standardized_patient_changes": patient_scaled,
        })
    return rows


def main():
    all_rows = []
    for method in METHOD_ORDER:
        score_file, test_file, eligible = FILES[method]
        all_rows.extend(load_method(method, score_file, test_file, eligible))
    table = pd.DataFrame(all_rows)
    table["feature_label"] = table["method"] + ": " + table["feature"]
    table.to_csv(CSV, index=False)

    heat = np.vstack(table["standardized_patient_changes"].to_list())
    row_count = len(table)
    fig = plt.figure(figsize=(17, max(17, .28 * row_count)))
    grid = GridSpec(1, 2, figure=fig, width_ratios=[3.4, 1.5], wspace=.04)
    ax_heat = fig.add_subplot(grid[0, 0])
    ax_summary = fig.add_subplot(grid[0, 1], sharey=ax_heat)
    cmap = sns.color_palette("vlag", as_cmap=True)
    cmap.set_bad("#d9d9d9")
    image = ax_heat.imshow(np.ma.masked_invalid(heat), aspect="auto", interpolation="nearest", cmap=cmap, norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5))
    ax_heat.set_xticks(np.arange(len(PATIENTS)))
    ax_heat.set_xticklabels(PATIENTS, rotation=0)
    ax_heat.set_yticks(np.arange(row_count))
    ax_heat.set_yticklabels(table["feature_label"], fontsize=7)
    ax_heat.set_xlabel("Patient")
    ax_heat.set_title("Patient-level tumour − normal change\n(row-standardized within each method-feature)")
    cbar = fig.colorbar(image, ax=ax_heat, fraction=.025, pad=.015)
    cbar.set_label("Within-feature SD units")

    summary_x = table["median_paired_change_sd"].to_numpy()
    nominal = table["p_value_nominal"].to_numpy() < .05
    ax_summary.axvline(0, color="black", linewidth=.7)
    ax_summary.scatter(summary_x[~nominal], np.where(~nominal)[0], color="#666666", s=17, label="Nominal p ≥ 0.05")
    ax_summary.scatter(summary_x[nominal], np.where(nominal)[0], color="#d1495b", s=24, label="Nominal p < 0.05")
    for index in np.where(nominal)[0]:
        ax_summary.text(summary_x[index], index, " *", va="center", fontsize=7, color="#8c1d2c")
    ax_summary.set_xlabel("Median paired change\n(within-feature SD units)")
    ax_summary.set_title("Paired-change summary\n(nominal Wilcoxon p-values)")
    ax_summary.tick_params(axis="y", left=False, labelleft=False)
    ax_summary.legend(loc="lower right", fontsize=8, frameon=False)

    # Separate method blocks and mark EPIC's converged-pair sensitivity subset.
    offsets = np.cumsum([sum(table.method.eq(method)) for method in METHOD_ORDER])
    for boundary in offsets[:-1]:
        ax_heat.axhline(boundary - .5, color="black", linewidth=1.0)
        ax_summary.axhline(boundary - .5, color="black", linewidth=1.0)
    epi_rows = table.method.eq("EPIC")
    if epi_rows.any():
        label = ax_heat.get_yticklabels()
        for index in np.where(epi_rows.to_numpy())[0]:
            label[index].set_text(label[index].get_text() + "†")
        ax_heat.set_yticklabels([item.get_text() for item in label], fontsize=7)
    fig.suptitle("Paired immune and stromal comparisons: revised 11-patient PDAC cohort", y=.995, fontsize=14, fontweight="bold")
    fig.text(.01, .002, "EPIC (†) uses the eight matched pairs for which both samples converged; other methods use 11 pairs. Red points denote nominal p<0.05; no feature survived method-wise BH FDR<0.05. Standardized changes aid display across native method scales.", fontsize=8, va="bottom")
    fig.subplots_adjust(left=.31, right=.98, top=.96, bottom=.04)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("features", row_count, "methods", table.groupby("method").size().to_dict(), "nominal", int(nominal.sum()), "FDR<0.05", int((table.FDR_within_method < .05).sum()))


if __name__ == "__main__":
    main()
