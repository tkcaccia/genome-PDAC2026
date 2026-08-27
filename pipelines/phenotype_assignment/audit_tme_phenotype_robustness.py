#!/usr/bin/env python3
"""Audit reproducibility and perturbation sensitivity of TME phenotype labels.

The script operates on private score/assignment tables but writes no results to
the code repository. It checks the recorded feature manifest, reconstructs the
ranked labels, compares a quantile rule, and measures label stability after
leaving out one feature or one method family at a time.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


GROUP_STROMAL = "StromalHigh_EMTHigh_ImmuneLow"
GROUP_IMMUNE = "ImmuneHigh_StromalLow"
GROUP_INTERMEDIATE = "Intermediate_or_mixed"
META_GROUPS = ("immune_meta_score", "stromal_meta_score", "emt_meta_score")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit TME phenotype assignment robustness.")
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--assignments", required=True, type=Path)
    parser.add_argument("--feature-manifest", required=True, type=Path)
    parser.add_argument("--sample-column", default="sample_id")
    parser.add_argument("--target-per-extreme", type=int, default=3)
    parser.add_argument("--upper-quantile", type=float, default=0.75)
    parser.add_argument("--out-prefix", required=True, type=Path)
    return parser.parse_args()


def zscore(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    standard_deviation = numeric.std(skipna=True, ddof=0)
    if pd.isna(standard_deviation) or standard_deviation == 0:
        return pd.Series(np.zeros(len(numeric)), index=numeric.index)
    return (numeric - numeric.mean(skipna=True)) / standard_deviation


def calculate_contrasts(
    scores: pd.DataFrame, feature_groups: dict[str, list[str]]
) -> pd.DataFrame:
    output = pd.DataFrame(index=scores.index)
    for meta_score, columns in feature_groups.items():
        if not columns:
            raise ValueError(f"No features remain for {meta_score}")
        standardized = pd.concat([zscore(scores[column]) for column in columns], axis=1)
        output[meta_score] = standardized.mean(axis=1, skipna=True)
    output["immune_high_stromal_low_score"] = (
        output["immune_meta_score"] - output["stromal_meta_score"]
    )
    output["stromal_emt_high_immune_low_score"] = (
        output[["stromal_meta_score", "emt_meta_score"]].mean(axis=1)
        - output["immune_meta_score"]
    )
    return output


def assign_by_rank(contrasts: pd.DataFrame, target: int) -> pd.Series:
    labels = pd.Series(GROUP_INTERMEDIATE, index=contrasts.index)
    immune_order = contrasts.sort_values(
        "immune_high_stromal_low_score", ascending=False
    ).index.tolist()
    stromal_order = contrasts.sort_values(
        "stromal_emt_high_immune_low_score", ascending=False
    ).index.tolist()
    immune_selected = immune_order[:target]
    stromal_selected = [index for index in stromal_order if index not in immune_selected][:target]
    labels.loc[immune_selected] = GROUP_IMMUNE
    labels.loc[stromal_selected] = GROUP_STROMAL
    return labels


def assign_by_quantile(contrasts: pd.DataFrame, quantile: float) -> pd.Series:
    labels = pd.Series(GROUP_INTERMEDIATE, index=contrasts.index)
    immune = contrasts["immune_high_stromal_low_score"]
    stromal = contrasts["stromal_emt_high_immune_low_score"]
    immune_mask = immune >= immune.quantile(quantile)
    stromal_mask = stromal >= stromal.quantile(quantile)
    labels.loc[immune_mask] = GROUP_IMMUNE
    labels.loc[stromal_mask] = GROUP_STROMAL
    both = immune_mask & stromal_mask
    labels.loc[both & (immune > stromal)] = GROUP_IMMUNE
    return labels


def method_family(column: str) -> str:
    for prefix in ("mcp_counter", "estimate", "epic", "xcell", "quantiseq", "programme"):
        if column.startswith(f"{prefix}_"):
            return prefix
    return column.split("_", maxsplit=1)[0]


def concordance(first: pd.Series, second: pd.Series) -> float:
    return float((first == second).mean())


def main() -> None:
    args = parse_args()
    scores = pd.read_csv(args.scores, sep="\t")
    assignments = pd.read_csv(args.assignments, sep="\t")
    manifest = pd.read_csv(args.feature_manifest, sep="\t")

    if set(manifest["meta_score"]) - set(META_GROUPS):
        raise SystemExit("Feature manifest contains an unknown meta-score name.")
    feature_groups = {
        meta_score: manifest.loc[manifest["meta_score"] == meta_score, "source_column"].tolist()
        for meta_score in META_GROUPS
    }
    selected = [column for columns in feature_groups.values() for column in columns]
    missing = [column for column in selected if column not in scores.columns]
    if missing:
        raise SystemExit(f"Selected features missing from score table: {', '.join(missing)}")
    if scores[selected].apply(pd.to_numeric, errors="coerce").isna().any().any():
        raise SystemExit("Selected phenotype features contain missing or non-numeric values.")

    assignment_columns = [args.sample_column, "phenotype_group"]
    recorded = assignments[assignment_columns].merge(
        scores[[args.sample_column]], on=args.sample_column, how="right", validate="one_to_one"
    )
    if recorded["phenotype_group"].isna().any():
        raise SystemExit("One or more score-table samples lack a recorded assignment.")
    recorded_labels = recorded.set_index(args.sample_column)["phenotype_group"].reindex(scores[args.sample_column])
    recorded_labels.index = scores.index

    baseline_contrasts = calculate_contrasts(scores, feature_groups)
    reconstructed = assign_by_rank(baseline_contrasts, args.target_per_extreme)
    quantile_labels = assign_by_quantile(baseline_contrasts, args.upper_quantile)

    perturbations: list[dict[str, object]] = []
    perturbation_labels: dict[str, pd.Series] = {}

    for meta_score, columns in feature_groups.items():
        for omitted in columns:
            reduced = {name: list(values) for name, values in feature_groups.items()}
            reduced[meta_score].remove(omitted)
            label = f"leave_feature:{omitted}"
            labels = assign_by_rank(calculate_contrasts(scores, reduced), args.target_per_extreme)
            perturbation_labels[label] = labels
            perturbations.append(
                {
                    "perturbation": label,
                    "type": "leave_one_feature_out",
                    "concordance_with_recorded": concordance(labels, recorded_labels),
                    "changed_samples": int((labels != recorded_labels).sum()),
                }
            )

    families = sorted({method_family(column) for column in selected})
    for family in families:
        reduced = {
            name: [column for column in columns if method_family(column) != family]
            for name, columns in feature_groups.items()
        }
        if any(not columns for columns in reduced.values()):
            perturbations.append(
                {
                    "perturbation": f"leave_family:{family}",
                    "type": "leave_one_family_out",
                    "concordance_with_recorded": np.nan,
                    "changed_samples": np.nan,
                }
            )
            continue
        label = f"leave_family:{family}"
        labels = assign_by_rank(calculate_contrasts(scores, reduced), args.target_per_extreme)
        perturbation_labels[label] = labels
        perturbations.append(
            {
                "perturbation": label,
                "type": "leave_one_family_out",
                "concordance_with_recorded": concordance(labels, recorded_labels),
                "changed_samples": int((labels != recorded_labels).sum()),
            }
        )

    perturbation_table = pd.DataFrame(perturbations)
    group_summary_input = baseline_contrasts.copy()
    group_summary_input["phenotype_group"] = recorded_labels
    group_summary = (
        group_summary_input.groupby("phenotype_group", sort=False)
        .agg(
            n=("immune_meta_score", "size"),
            mean_immune_meta_score=("immune_meta_score", "mean"),
            mean_stromal_meta_score=("stromal_meta_score", "mean"),
            mean_emt_meta_score=("emt_meta_score", "mean"),
            mean_immune_contrast=("immune_high_stromal_low_score", "mean"),
            mean_stromal_emt_contrast=("stromal_emt_high_immune_low_score", "mean"),
        )
        .reset_index()
    )
    stability = pd.DataFrame({args.sample_column: scores[args.sample_column]})
    if perturbation_labels:
        matches = pd.DataFrame(
            {name: labels == recorded_labels for name, labels in perturbation_labels.items()}
        )
        stability["stable_fraction"] = matches.mean(axis=1)
        stability["changed_in_n_perturbations"] = (~matches).sum(axis=1)
        stability["perturbations_tested"] = matches.shape[1]

    summary_rows = [
        ("samples", len(scores)),
        ("selected_features", len(selected)),
        ("selected_features_with_missing_values", int(scores[selected].isna().any().sum())),
        ("exact_ranked_label_reproduction", str(bool((reconstructed == recorded_labels).all())).lower()),
        ("ranked_label_concordance", f"{concordance(reconstructed, recorded_labels):.6f}"),
        ("quantile_label_concordance", f"{concordance(quantile_labels, recorded_labels):.6f}"),
        ("quantile_immune_n", int((quantile_labels == GROUP_IMMUNE).sum())),
        ("quantile_stromal_emt_n", int((quantile_labels == GROUP_STROMAL).sum())),
        ("quantile_intermediate_n", int((quantile_labels == GROUP_INTERMEDIATE).sum())),
        ("valid_perturbations", len(perturbation_labels)),
        ("minimum_sample_stability", f"{stability['stable_fraction'].min():.6f}"),
        ("median_sample_stability", f"{stability['stable_fraction'].median():.6f}"),
        ("minimum_perturbation_concordance", f"{perturbation_table['concordance_with_recorded'].min():.6f}"),
        ("median_perturbation_concordance", f"{perturbation_table['concordance_with_recorded'].median():.6f}"),
    ]
    summary = pd.DataFrame(summary_rows, columns=["metric", "value"])

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.out_prefix.with_suffix(".summary.tsv"), sep="\t", index=False)
    perturbation_table.to_csv(
        args.out_prefix.with_suffix(".perturbation_concordance.tsv"), sep="\t", index=False
    )
    group_summary.to_csv(
        args.out_prefix.with_suffix(".group_score_summary.tsv"), sep="\t", index=False
    )
    stability.to_csv(args.out_prefix.with_suffix(".sample_stability.tsv"), sep="\t", index=False)
    print(summary.to_csv(sep="\t", index=False), end="")


if __name__ == "__main__":
    main()
