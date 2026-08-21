#!/usr/bin/env python3
"""Patient-data-safe example for identifying TP53 mutation evidence.

This script mirrors the logic used in the PDAC2026 downstream review:
start from a VEP-annotated somatic variant table, extract TP53 rows, classify
protein-impacting coding variants, and create a patient-level summary.

It deliberately does not contain real patient data. The input file is supplied
by the user and should remain outside GitHub if it contains patient results.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PROTEIN_IMPACTING_TERMS = {
    "missense_variant",
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost",
    "inframe_deletion",
    "inframe_insertion",
    "protein_altering_variant",
}

REQUIRED_COLUMNS = {
    "patient_id",
    "caller",
    "chrom",
    "pos",
    "ref",
    "alt",
    "filter",
    "gene",
    "impact",
    "consequence",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract and summarize TP53 mutation evidence from a VEP-annotated somatic variant table."
    )
    parser.add_argument(
        "--variants",
        required=True,
        type=Path,
        help="Tab-separated somatic variant table annotated with VEP gene/consequence columns.",
    )
    parser.add_argument(
        "--gene",
        default="TP53",
        help="Driver gene symbol to review. Default: TP53.",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        type=Path,
        help="Output prefix. The script writes <prefix>.tp53_variant_rows.tsv and <prefix>.tp53_patient_summary.tsv.",
    )
    return parser.parse_args()


def require_columns(df: pd.DataFrame) -> None:
    missing = sorted(REQUIRED_COLUMNS - set(df.columns))
    if missing:
        raise SystemExit(
            "Missing required columns: "
            + ", ".join(missing)
            + "\nPlease adapt the column names or convert your VEP export before running this example."
        )


def is_protein_impacting(consequence: object, impact: object) -> bool:
    """Return True for coding variants likely to affect the protein sequence."""
    consequence_terms = str(consequence).split("&")
    has_term = any(term in PROTEIN_IMPACTING_TERMS for term in consequence_terms)
    has_high_or_moderate_impact = str(impact).upper() in {"HIGH", "MODERATE"}
    return has_term and has_high_or_moderate_impact


def is_pass(filter_value: object) -> bool:
    return str(filter_value).upper() == "PASS"


def summarize_patient(patient_id: str, rows: pd.DataFrame) -> dict[str, object]:
    coding = rows[rows["is_protein_impacting"]]
    coding_pass = coding[coding["is_pass"]]

    if len(coding_pass) > 0:
        status = "TP53_somatic_coding_PASS"
        evidence = coding_pass
    elif len(coding) > 0:
        status = "TP53_somatic_coding_review"
        evidence = coding
    else:
        status = "TP53_no_coding_somatic_variant_detected"
        evidence = rows.head(0)

    def join_unique(column: str) -> str:
        if column not in evidence.columns or evidence.empty:
            return ""
        values = evidence[column].dropna().astype(str)
        values = values[values != ""].unique().tolist()
        return ";".join(values)

    return {
        "patient_id": patient_id,
        "tp53_somatic_mutation_status": status,
        "tp53_rows_total": len(rows),
        "tp53_protein_impacting_rows": len(coding),
        "tp53_pass_protein_impacting_rows": len(coding_pass),
        "callers": join_unique("caller"),
        "filters": join_unique("filter"),
        "consequences": join_unique("consequence"),
        "amino_acids": join_unique("amino_acids"),
        "protein_positions": join_unique("protein_position"),
        "clin_sig": join_unique("clin_sig"),
        "known_variant_ids": join_unique("existing_variation"),
    }


def main() -> None:
    args = parse_args()
    variants = pd.read_csv(args.variants, sep="\t", dtype=str).fillna("")
    require_columns(variants)

    gene = args.gene.upper()
    gene_rows = variants[variants["gene"].str.upper() == gene].copy()
    gene_rows["is_protein_impacting"] = [
        is_protein_impacting(cons, impact)
        for cons, impact in zip(gene_rows["consequence"], gene_rows["impact"])
    ]
    gene_rows["is_pass"] = [is_pass(value) for value in gene_rows["filter"]]

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    variant_out = args.out_prefix.with_suffix(".tp53_variant_rows.tsv")
    summary_out = args.out_prefix.with_suffix(".tp53_patient_summary.tsv")

    gene_rows.to_csv(variant_out, sep="\t", index=False)

    summaries = [
        summarize_patient(patient_id, rows)
        for patient_id, rows in gene_rows.groupby("patient_id", sort=True)
    ]
    summary = pd.DataFrame(summaries)
    if summary.empty:
        summary = pd.DataFrame(
            columns=[
                "patient_id",
                "tp53_somatic_mutation_status",
                "tp53_rows_total",
                "tp53_protein_impacting_rows",
                "tp53_pass_protein_impacting_rows",
                "callers",
                "filters",
                "consequences",
                "amino_acids",
                "protein_positions",
                "clin_sig",
                "known_variant_ids",
            ]
        )
    summary.to_csv(summary_out, sep="\t", index=False)

    print(f"Wrote TP53 variant rows: {variant_out}")
    print(f"Wrote TP53 patient summary: {summary_out}")


if __name__ == "__main__":
    main()

