#!/usr/bin/env python3
"""Merge method-specific RNA expression score matrices by sample.

Inputs may have features in rows and samples in columns, or already have one
sample per row. The resulting table always has one row per sample and prefixed,
machine-safe score-column names for phenotype assignment.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assemble immune, stromal and programme score matrices into one sample-level table."
    )
    parser.add_argument(
        "--score-table",
        action="append",
        default=[],
        metavar="PREFIX=PATH",
        help="Feature-by-sample score matrix. Repeat for each method.",
    )
    parser.add_argument(
        "--programme-table",
        type=Path,
        help="Optional GSVA/ssGSEA programme-by-sample score matrix.",
    )
    parser.add_argument("--metadata", type=Path, help="Optional sample metadata used to retain tumours only.")
    parser.add_argument("--sample-column", default="sample_id")
    parser.add_argument("--patient-column", default="patient_id")
    parser.add_argument("--condition-column", default="condition")
    parser.add_argument("--tumour-label", default="Tumour")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def slug(value: object) -> str:
    text = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", str(value).strip())
    text = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", text)
    text = text.lower().replace("+", "_plus_")
    text = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    return text or "unnamed_feature"


def parse_spec(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise SystemExit(f"Invalid --score-table value '{spec}'; expected PREFIX=PATH.")
    prefix, path = spec.split("=", 1)
    if not prefix.strip() or not path.strip():
        raise SystemExit(f"Invalid --score-table value '{spec}'; expected PREFIX=PATH.")
    return slug(prefix), Path(path)


def read_score_matrix(path: Path, prefix: str, sample_column: str) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    if table.shape[1] < 2:
        raise SystemExit(f"Score matrix has no sample columns: {path}")

    # Archived result tables may already be sample-by-feature. Keep numeric
    # feature columns and drop metadata such as patient_id and condition.
    if table.columns[0] == sample_column or sample_column in table.columns:
        sample_values = table[sample_column].astype(str)
        metadata_columns = {
            sample_column,
            "patient_id",
            "patient_num",
            "condition",
            "condition_order",
            "phenotype_group",
            "notes",
        }
        candidate_columns = [column for column in table.columns if column not in metadata_columns]
        feature_table = table[candidate_columns].apply(pd.to_numeric, errors="coerce")
        feature_table = feature_table.loc[:, feature_table.notna().any(axis=0)]
        feature_table.columns = [f"{prefix}_{slug(value)}" for value in feature_table.columns]
        if feature_table.columns.duplicated().any():
            duplicates = sorted(set(feature_table.columns[feature_table.columns.duplicated()]))
            raise SystemExit(f"Duplicate normalized feature names in {path}: {', '.join(duplicates)}")
        feature_table.insert(0, sample_column, sample_values)
        return feature_table

    # Reusable method runners write feature-by-sample matrices. Transpose these
    # so that all downstream joins use one row per sample.
    feature_column = table.columns[0]
    features = [f"{prefix}_{slug(value)}" for value in table[feature_column]]
    if len(features) != len(set(features)):
        duplicates = sorted({name for name in features if features.count(name) > 1})
        raise SystemExit(f"Duplicate normalized feature names in {path}: {', '.join(duplicates)}")

    matrix = table.drop(columns=[feature_column]).apply(pd.to_numeric, errors="coerce")
    matrix.index = features
    sample_table = matrix.transpose().reset_index().rename(columns={"index": sample_column})
    return sample_table


def main() -> None:
    args = parse_args()
    specifications = [parse_spec(spec) for spec in args.score_table]
    if args.programme_table:
        specifications.append(("programme", args.programme_table))
    if not specifications:
        raise SystemExit("Supply at least one --score-table or --programme-table.")

    merged: pd.DataFrame | None = None
    for prefix, path in specifications:
        current = read_score_matrix(path, prefix, args.sample_column)
        merged = current if merged is None else merged.merge(current, on=args.sample_column, how="outer")

    assert merged is not None
    if args.metadata:
        metadata = pd.read_csv(args.metadata, sep="\t", dtype=str)
        required = {args.sample_column, args.condition_column}
        missing = required.difference(metadata.columns)
        if missing:
            raise SystemExit(f"Metadata is missing columns: {', '.join(sorted(missing))}")
        metadata = metadata.loc[
            metadata[args.condition_column] == args.tumour_label
        ].copy()
        metadata_columns = [args.sample_column, args.condition_column]
        if args.patient_column in metadata.columns:
            metadata_columns.insert(1, args.patient_column)
        merged = metadata[metadata_columns].merge(
            merged,
            on=args.sample_column,
            how="inner",
            validate="one_to_one",
        )

    merged = merged.sort_values(args.sample_column)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(merged)} sample rows and {len(merged.columns) - 1} score columns: {args.output}")


if __name__ == "__main__":
    main()
