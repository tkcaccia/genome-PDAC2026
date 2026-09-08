#!/usr/bin/env python3
"""Create a publication figure from the WES depth and call-audit tables.

The plot keeps target depth, matched-callable territory and somatic call burden
as separate quantities.  It deliberately does not relabel the full-target call
burden as tumour mutational burden (TMB), because a clinical TMB calculation
requires a validated coding-territory denominator.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-qc", type=Path, required=True)
    parser.add_argument("--patient-summary", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument(
        "--label-map",
        type=Path,
        help="Optional TSV with patient_id and display_label columns.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric_patient_key(value: str) -> tuple[int, str]:
    try:
        return int(value), value
    except ValueError:
        return 10**9, value


def as_float(value: str | None) -> float:
    if value is None or value.strip().upper() in {"", "NA", "NAN"}:
        return float("nan")
    return float(value)


def as_bool(value: str | None) -> bool:
    return (value or "").strip().lower() in {"1", "true", "yes"}


def build_labels(patient_ids: list[str], label_map: Path | None) -> list[str]:
    if label_map is None:
        return [f"Tumour {index:02d}" for index in range(1, len(patient_ids) + 1)]

    mapping = {
        row["patient_id"]: row["display_label"]
        for row in read_tsv(label_map)
        if row.get("patient_id") and row.get("display_label")
    }
    missing = [patient_id for patient_id in patient_ids if patient_id not in mapping]
    if missing:
        raise ValueError(f"Missing display labels for patient IDs: {', '.join(missing)}")
    return [mapping[patient_id] for patient_id in patient_ids]


def main() -> None:
    args = parse_args()
    pair_rows = sorted(
        read_tsv(args.pair_qc),
        key=lambda row: numeric_patient_key(row["patient_id"]),
    )
    summary_by_patient = {
        row["patient_id"]: row for row in read_tsv(args.patient_summary)
    }

    patient_ids = [row["patient_id"] for row in pair_rows]
    labels = build_labels(patient_ids, args.label_map)
    x = np.arange(len(patient_ids))

    tumour_depth = np.array(
        [as_float(row["tumour_mean_target_depth_x"]) for row in pair_rows]
    )
    normal_depth = np.array(
        [as_float(row["normal_mean_target_depth_x"]) for row in pair_rows]
    )
    callable_pct = np.array(
        [as_float(row["pair_callable_pct_of_target"]) for row in pair_rows]
    )
    strict_burden = np.array(
        [
            as_float(
                summary_by_patient[patient_id].get(
                    "strict_read_evidence_variants_per_callable_target_mb"
                )
            )
            for patient_id in patient_ids
        ]
    )
    outlier = np.array(
        [
            as_bool(
                summary_by_patient[patient_id].get("strict_burden_cohort_outlier")
            )
            for patient_id in patient_ids
        ]
    )
    missing_calls = np.isnan(strict_burden)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 8.5), constrained_layout=True)

    for index in range(len(x)):
        axes[0].plot(
            [x[index], x[index]],
            [normal_depth[index], tumour_depth[index]],
            color="#CCD4D8",
            linewidth=1,
            zorder=1,
        )
    axes[0].scatter(x, tumour_depth, color="#B43C2F", s=34, label="Tumour", zorder=3)
    axes[0].scatter(x, normal_depth, color="#1F6585", s=34, label="Matched normal", zorder=3)
    axes[0].set_ylabel("Mean target depth (x)")
    axes[0].set_title("A  Targeted WES depth")
    axes[0].legend(frameon=False, ncol=2, loc="upper left")

    callable_colors = np.where(callable_pct < 70.0, "#D08B2D", "#3D7D6B")
    axes[1].bar(x, callable_pct, color=callable_colors, width=0.72)
    axes[1].axhline(70.0, color="#333333", linewidth=1, linestyle="--")
    axes[1].text(
        len(x) - 0.4,
        71.5,
        "Descriptive 70% review threshold",
        ha="right",
        va="bottom",
        fontsize=8,
    )
    axes[1].set_ylim(0, 105)
    axes[1].set_ylabel("Paired callable target (%)")
    axes[1].set_title("B  Territory callable in both tumour and normal")

    burden_colors = np.where(outlier, "#B43C2F", "#6E8791")
    burden_colors = np.where(missing_calls, "#D7D7D7", burden_colors)
    axes[2].bar(x, np.nan_to_num(strict_burden, nan=0.0), color=burden_colors, width=0.72)
    axes[2].set_yscale("symlog", linthresh=0.03, linscale=0.7)
    axes[2].set_ylabel("Strict variants per callable\ntarget Mb")
    axes[2].set_title("C  Strict read-evidence call burden")
    axes[2].text(
        0.01,
        0.98,
        "Full callable target denominator; not clinical TMB",
        transform=axes[2].transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    for index, is_outlier in enumerate(outlier):
        if is_outlier:
            axes[2].annotate(
                "Cohort outlier",
                (x[index], strict_burden[index]),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7.5,
                color="#8C2D24",
            )
    for index, is_missing in enumerate(missing_calls):
        if is_missing:
            axes[2].scatter(
                x[index],
                0.03,
                marker="x",
                s=42,
                linewidth=1.5,
                color="#555555",
                zorder=4,
            )
            axes[2].annotate(
                "Final calls\nnot recovered",
                (x[index], 0.03),
                xytext=(0, 7),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
                color="#555555",
            )

    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(labels, rotation=45, ha="right")
        axis.grid(axis="y", color="#E7E7E7", linewidth=0.7)
        axis.set_axisbelow(True)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix(".png"), dpi=320, bbox_inches="tight")
    fig.savefig(args.output_prefix.with_suffix(".pdf"), bbox_inches="tight")


if __name__ == "__main__":
    main()
