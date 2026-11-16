#!/usr/bin/env python3
# WES processing, variant summaries, coverage and paired MSI.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""
from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

ROOT = private_path("handoff_root")


def category(score):
    if not math.isfinite(score) or not 0 <= score <= 100:
        raise ValueError(f"Invalid MSI percentage: {score}")
    return "MSI-High" if score >= 20 else "MSI-Intermediate" if score >= 5 else "MSS"


def main():
    cohorts = pd.read_csv(ROOT / "manifests/patient_cohorts.tsv", sep="\t", dtype=str)
    rows = []
    for patient in cohorts.itertuples(index=False):
        source = ROOT / "wes/msi_2_5" / f"patient_{patient.patient_id}"
        table = pd.read_csv(source, sep=r"\s+")
        if len(table) != 1:
            raise ValueError(f"Expected one summary row: {source}")
        result = table.iloc[0]
        total = int(result["Total_Number_of_Sites"])
        unstable = int(result["Number_of_Unstable_Sites"])
        score = float(result["%"])
        if not 0 <= unstable <= total or total == 0:
            raise ValueError(f"Invalid site counts in {source}")
        if abs(score - 100 * unstable / total) > 0.011:
            raise ValueError(f"Score/count disagreement in {source}")
        rows.append({"patient": patient.patient_id, "cohort": patient.cohort,
                     "informative_sites": total, "unstable_sites": unstable,
                     "msi_score_pct": score, "category": category(score)})
    results = pd.DataFrame(rows)
    out = ROOT / "wes/msi_2_5"
    results.to_csv(out / "msi_scores_all_14.csv", index=False)
    for cohort, subset in results.groupby("cohort"):
        subset.to_csv(out / f"msi_scores_{cohort}.csv", index=False)
    colors = {"MSS": "#2a9d8f", "MSI-Intermediate": "#e9c46a", "MSI-High": "#d62728"}
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.stripplot(data=results, x="cohort", y="msi_score_pct", hue="category",
                  order=["PDAC", "periampullary"], hue_order=["MSS", "MSI-Intermediate", "MSI-High"],
                  dodge=False, jitter=0.15, palette=colors, size=9, ax=ax)
    ax.axhline(5, color="grey", linestyle="--")
    ax.axhline(20, color="grey", linestyle="--")
    ax.set(xlabel="Cohort", ylabel="MSI score (% unstable microsatellite sites)",
           title="Paired MSI scores: 11 PDAC and 3 exploratory periampullary cases")
    fig.tight_layout()
    for extension in ("png", "pdf", "svg"):
        fig.savefig(out / f"msi_scores_all_14.{extension}", dpi=300)
    plt.close(fig)
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()