#!/usr/bin/env python3
"""Assign cohort-relative immune/stromal/EMT tumour phenotype groups.

This is a patient-data-safe example. It expects a tumour-level score table
containing immune, stromal/fibroblast/CAF/ECM and EMT/invasion score columns
generated upstream from normalized RNA-seq expression data.

The script standardizes each selected feature across tumours, averages related
features into meta-scores, and assigns broad phenotype labels from contrast
scores. This mirrors the logic used to define the exploratory PDAC2026 groups:

* StromalHigh_EMTHigh_ImmuneLow
* ImmuneHigh_StromalLow
* Intermediate
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


GROUP_STROMAL_EMT = "StromalHigh_EMTHigh_ImmuneLow"
GROUP_IMMUNE = "ImmuneHigh_StromalLow"
GROUP_INTERMEDIATE = "Intermediate"


def comma_list(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assign cohort-relative tumour immune/stromal/EMT phenotype groups."
    )
    parser.add_argument("--scores", required=True, type=Path, help="Tumour-level score table in TSV format.")
    parser.add_argument("--sample-column", default="sample_id", help="Column containing tumour/sample IDs.")
    parser.add_argument("--immune-columns", required=True, help="Comma-separated immune score columns.")
    parser.add_argument("--stromal-columns", required=True, help="Comma-separated stromal/CAF/ECM score columns.")
    parser.add_argument("--emt-columns", required=True, help="Comma-separated EMT/invasion score columns.")
    parser.add_argument(
        "--upper-quantile",
        type=float,
        default=0.75,
        help="Quantile used to define high contrast scores. Default: 0.75.",
    )
    parser.add_argument(
        "--target-per-extreme",
        type=int,
        default=None,
        help=(
            "Optional: assign exactly this many tumours to each extreme group by ranking contrast scores. "
            "For PDAC2026, the final reviewed split used three tumours per extreme group."
        ),
    )
    parser.add_argument("--out-prefix", required=True, type=Path, help="Output prefix.")
    return parser.parse_args()


def check_columns(df: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise SystemExit(f"Missing {label} columns: {', '.join(missing)}")


def zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    mean = values.mean(skipna=True)
    sd = values.std(skipna=True, ddof=0)
    if pd.isna(sd) or sd == 0:
        return pd.Series(np.zeros(len(values)), index=series.index)
    return (values - mean) / sd


def add_meta_score(df: pd.DataFrame, columns: list[str], prefix: str) -> pd.Series:
    z_cols = []
    for col in columns:
        z_col = f"z_{prefix}_{col}"
        df[z_col] = zscore(df[col])
        z_cols.append(z_col)
    return df[z_cols].mean(axis=1, skipna=True)


def assign_by_quantile(df: pd.DataFrame, upper_quantile: float) -> pd.Series:
    immune_threshold = df["immune_high_stromal_low_score"].quantile(upper_quantile)
    stromal_threshold = df["stromal_emt_high_immune_low_score"].quantile(upper_quantile)
    labels = pd.Series(GROUP_INTERMEDIATE, index=df.index)

    immune_mask = df["immune_high_stromal_low_score"] >= immune_threshold
    stromal_mask = df["stromal_emt_high_immune_low_score"] >= stromal_threshold

    # If a tumour meets both high-contrast rules, keep the stronger contrast.
    labels.loc[immune_mask] = GROUP_IMMUNE
    labels.loc[stromal_mask] = GROUP_STROMAL_EMT
    both = immune_mask & stromal_mask
    labels.loc[both & (df["immune_high_stromal_low_score"] > df["stromal_emt_high_immune_low_score"])] = GROUP_IMMUNE
    return labels


def assign_by_rank(df: pd.DataFrame, target_per_extreme: int) -> pd.Series:
    labels = pd.Series(GROUP_INTERMEDIATE, index=df.index)
    immune_order = df.sort_values("immune_high_stromal_low_score", ascending=False).index.tolist()
    stromal_order = df.sort_values("stromal_emt_high_immune_low_score", ascending=False).index.tolist()

    immune_selected: list[int] = []
    stromal_selected: list[int] = []

    for idx in immune_order:
        if len(immune_selected) >= target_per_extreme:
            break
        immune_selected.append(idx)

    for idx in stromal_order:
        if len(stromal_selected) >= target_per_extreme:
            break
        if idx in immune_selected:
            # Avoid assigning the same tumour to both extremes.
            continue
        stromal_selected.append(idx)

    labels.loc[immune_selected] = GROUP_IMMUNE
    labels.loc[stromal_selected] = GROUP_STROMAL_EMT
    return labels


def main() -> None:
    args = parse_args()
    immune_columns = comma_list(args.immune_columns)
    stromal_columns = comma_list(args.stromal_columns)
    emt_columns = comma_list(args.emt_columns)

    scores = pd.read_csv(args.scores, sep="\t")
    check_columns(scores, [args.sample_column], "sample ID")
    check_columns(scores, immune_columns, "immune")
    check_columns(scores, stromal_columns, "stromal")
    check_columns(scores, emt_columns, "EMT")

    out = scores[[args.sample_column]].copy()
    work = scores.copy()
    out["immune_meta_score"] = add_meta_score(work, immune_columns, "immune")
    out["stromal_meta_score"] = add_meta_score(work, stromal_columns, "stromal")
    out["emt_meta_score"] = add_meta_score(work, emt_columns, "emt")
    out["immune_high_stromal_low_score"] = out["immune_meta_score"] - out["stromal_meta_score"]
    out["stromal_emt_high_immune_low_score"] = (
        out[["stromal_meta_score", "emt_meta_score"]].mean(axis=1) - out["immune_meta_score"]
    )

    if args.target_per_extreme is not None:
        out["phenotype_group"] = assign_by_rank(out, args.target_per_extreme)
        out["assignment_method"] = f"ranked_top_{args.target_per_extreme}_per_extreme"
    else:
        out["phenotype_group"] = assign_by_quantile(out, args.upper_quantile)
        out["assignment_method"] = f"upper_quantile_{args.upper_quantile}"

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    assignment_path = args.out_prefix.with_suffix(".assignments.tsv")
    counts_path = args.out_prefix.with_suffix(".group_counts.tsv")

    out.to_csv(assignment_path, sep="\t", index=False)
    (
        out["phenotype_group"]
        .value_counts()
        .rename_axis("phenotype_group")
        .reset_index(name="n")
        .to_csv(counts_path, sep="\t", index=False)
    )

    print(f"Wrote assignments: {assignment_path}")
    print(f"Wrote group counts: {counts_path}")


if __name__ == "__main__":
    main()

