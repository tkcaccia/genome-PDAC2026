#!/usr/bin/env python3
# Separate PDAC and periampullary comparisons.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""
from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import pandas as pd

ROOT = private_path("handoff_root")
OUT_DIR = ROOT / "comparisons/pdac_ampulla_2_5/figures"
OUT = OUT_DIR / "Figure_7A_kras_mutation.png"
STATUS_OUT = OUT_DIR / "pdac_ampulla_kras_status.png"
VARIANTS = ROOT / "wes/analysis/tables/kras_variants_all_cohorts.tsv"
PATIENTS = configured_cohort("primary_pdac") + configured_cohort("periampullary")
PDAC = set(PATIENTS[:11])
COLORS = {"G12V": "#4C78A8", "Q61L": "#F58518", "G12R": "#E45756", "G12D": "#54A24B", "No detected hotspot": "#EEEEEE"}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    variants = pd.read_csv(VARIANTS, sep="\t")
    variants["patient"] = variants["patient"].astype(str)
    status = {patient: "No detected hotspot" for patient in PATIENTS}
    for _, row in variants.iterrows():
        patient = str(row["patient"])
        if patient not in status:
            continue
        protein = {"A": "G12V", "T": "G12D", "G": "G12R"}.get(str(row["ALT"]), None) if int(row["POS"]) in {25245350, 25245351} else None
        if int(row["POS"]) == 25227342:
            protein = "Q61L"
        status[patient] = protein or "Mutated"

    fig, ax = plt.subplots(figsize=(12, 3.8))
    for index, patient in enumerate(PATIENTS):
        category = status[patient]
        color = COLORS.get(category, "#984EA3")
        ax.add_patch(Rectangle((index, 0), 1, 1, facecolor=color, edgecolor="white", linewidth=1.5))
        label = category if category != "No detected hotspot" else "No hotspot"
        if category != "No detected hotspot":
            ax.text(index + .5, .5, label, ha="center", va="center", fontsize=9, color="white", fontweight="bold")
    ax.axvline(len(PDAC), color="black", linewidth=2)
    ax.set_xlim(0, len(PATIENTS)); ax.set_ylim(0, 1)
    ax.set_xticks([i + .5 for i in range(len(PATIENTS))], [f"{p}T" for p in PATIENTS])
    ax.set_yticks([.5], ["KRAS"])
    ax.tick_params(axis="x", length=0); ax.tick_params(axis="y", length=0)
    ax.set_title("KRAS mutation status: PDAC versus ampulla\nExploratory patient-level comparison", pad=28)
    ax.text(len(PDAC) / 2, -0.16, "PDAC (n=11)", transform=ax.get_xaxis_transform(), ha="center", va="top", fontweight="bold")
    ax.text(len(PDAC) + len(PATIENTS[len(PDAC):]) / 2, -0.16, "Ampulla (n=3)", transform=ax.get_xaxis_transform(), ha="center", va="top", fontweight="bold")
    ax.legend(handles=[Patch(facecolor=color, edgecolor="white", label=label) for label, color in COLORS.items()], loc="upper center", bbox_to_anchor=(.5, -.31), ncol=5, frameon=False)
    for spine in ax.spines.values(): spine.set_visible(False)
    fig.tight_layout(rect=(0, .08, 1, .95))
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    fig.savefig(STATUS_OUT, dpi=300, bbox_inches="tight")
    for extension in ("pdf", "svg"):
        fig.savefig(OUT.with_suffix(f".{extension}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("status:", status)


if __name__ == "__main__":
    main()
