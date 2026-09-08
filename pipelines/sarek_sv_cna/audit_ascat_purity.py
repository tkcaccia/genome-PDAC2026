#!/usr/bin/env python3
"""Audit ASCAT purity estimates against retained orthogonal metadata.

This script does not decide whether an ASCAT solution is biologically correct.
It identifies boundary and weakly supported solutions that require inspection of
the original ASCAT LogR/BAF profiles or an independent purity caller.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


PURITY_COLUMN = "ASCAT_aberrant_cell_fraction_purity_estimate"
PLOIDY_COLUMN = "ASCAT_ploidy"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit retained ASCAT purity/ploidy summaries."
    )
    parser.add_argument("--molecular-table", required=True, type=Path)
    parser.add_argument("--clinical-table", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value: str | None) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if not text or text.lower() in {"na", "nan", "none"}:
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def tumour_fraction_from_ratio(value: str | None) -> float | None:
    """Parse a tumour/stroma entry such as '40/60' as 0.40 tumour fraction."""
    if not value:
        return None
    text = value.strip().replace("%", "")
    for separator in ("/", ":"):
        if separator in text:
            left, right = text.split(separator, 1)
            tumour = as_float(left)
            stroma = as_float(right)
            if tumour is None or stroma is None or tumour + stroma <= 0:
                return None
            return tumour / (tumour + stroma)
    number = as_float(text)
    if number is None:
        return None
    return number / 100.0 if number > 1 else number


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def mean(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def rank(values: list[float]) -> list[float]:
    """Return average ranks for ties without requiring SciPy."""
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][1] == ordered[index][1]:
            end += 1
        average_rank = (index + 1 + end) / 2.0
        for position in range(index, end):
            ranks[ordered[position][0]] = average_rank
        index = end
    return ranks


def pearson(left: list[float], right: list[float]) -> float | None:
    if len(left) < 3 or len(left) != len(right):
        return None
    left_mean = mean(left)
    right_mean = mean(right)
    assert left_mean is not None and right_mean is not None
    numerator = sum((x - left_mean) * (y - right_mean) for x, y in zip(left, right))
    left_ss = sum((x - left_mean) ** 2 for x in left)
    right_ss = sum((y - right_mean) ** 2 for y in right)
    denominator = math.sqrt(left_ss * right_ss)
    return numerator / denominator if denominator else None


def spearman(left: list[float], right: list[float]) -> float | None:
    return pearson(rank(left), rank(right))


def fmt(value: float | None, digits: int = 3) -> str:
    return "NA" if value is None else f"{value:.{digits}f}"


def numeric_values(rows: list[dict[str, object]], key: str) -> list[float]:
    values = []
    for row in rows:
        value = row.get(key)
        if isinstance(value, (int, float)) and math.isfinite(float(value)):
            values.append(float(value))
    return values


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    molecular = read_tsv(args.molecular_table)
    clinical = read_tsv(args.clinical_table)
    clinical_by_patient = {row.get("patient_id", "").strip(): row for row in clinical}

    required = {"patient_id", PURITY_COLUMN, PLOIDY_COLUMN}
    if not molecular or not required.issubset(molecular[0]):
        missing = sorted(required - set(molecular[0] if molecular else {}))
        raise SystemExit(f"Molecular table is missing required columns: {missing}")

    audited: list[dict[str, object]] = []
    for row in molecular:
        patient_id = row.get("patient_id", "").strip()
        clinical_row = clinical_by_patient.get(patient_id, {})
        purity = as_float(row.get(PURITY_COLUMN))
        ploidy = as_float(row.get(PLOIDY_COLUMN))
        tumour_fraction = tumour_fraction_from_ratio(
            clinical_row.get("Tumor_stroma_ratio")
        )
        exact_upper_boundary = purity is not None and purity >= 0.9995
        near_diploid = ploidy is not None and abs(ploidy - 2.0) <= 0.10
        reasons = []
        if exact_upper_boundary:
            reasons.append("purity_at_upper_boundary")
        if near_diploid:
            reasons.append("near_diploid_solution")
        if tumour_fraction is None:
            reasons.append("no_pathology_tumour_fraction")
        elif purity is not None and abs(purity - tumour_fraction) >= 0.30:
            reasons.append("discordant_with_pathology_fraction")
        if not row.get("ASCAT_purity_ploidy_source", "").strip():
            reasons.append("missing_source_provenance")

        audited.append(
            {
                "patient_id": patient_id,
                "ascat_purity": purity,
                "ascat_ploidy": ploidy,
                "pathology_tumour_fraction": tumour_fraction,
                "purity_minus_pathology": (
                    purity - tumour_fraction
                    if purity is not None and tumour_fraction is not None
                    else None
                ),
                "exact_upper_boundary": exact_upper_boundary,
                "near_diploid": near_diploid,
                "large_CNV_event_count": as_float(row.get("large_CNV_event_count")),
                "cnv_burden_altered_segments": as_float(
                    row.get("cnv_burden_altered_segments")
                ),
                "segments": as_float(row.get("segments")),
                "SV_total_records": as_float(row.get("SV_total_records")),
                "tumour_pass_mutations": as_float(row.get("tumour_pass_mutations")),
                "tumour_wes_raw_total_reads": as_float(
                    row.get("tumour_wes_raw_total_reads")
                ),
                "normal_wes_raw_total_reads": as_float(
                    row.get("normal_wes_raw_total_reads")
                ),
                "tumour_wes_reads_duplicated_percent": as_float(
                    row.get("tumour_wes_reads_duplicated_percent")
                ),
                "normal_wes_reads_duplicated_percent": as_float(
                    row.get("normal_wes_reads_duplicated_percent")
                ),
                "tumour_wes_reads_mapped_percent": as_float(
                    row.get("tumour_wes_reads_mapped_percent")
                ),
                "normal_wes_reads_mapped_percent": as_float(
                    row.get("normal_wes_reads_mapped_percent")
                ),
                "diagnostic_flags": ";".join(reasons) or "none",
            }
        )

    patient_output = args.output_dir / "ascat_purity_patient_level_audit.tsv"
    with patient_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audited[0]), delimiter="\t")
        writer.writeheader()
        for row in audited:
            writer.writerow(
                {
                    key: (
                        ""
                        if value is None
                        else f"{value:.6g}"
                        if isinstance(value, float)
                        else str(value).lower()
                        if isinstance(value, bool)
                        else value
                    )
                    for key, value in row.items()
                }
            )

    boundary = [row for row in audited if row["exact_upper_boundary"]]
    interior = [row for row in audited if not row["exact_upper_boundary"]]
    paired = [
        row
        for row in audited
        if row["ascat_purity"] is not None
        and row["pathology_tumour_fraction"] is not None
    ]
    purities = numeric_values(audited, "ascat_purity")
    pathology = numeric_values(paired, "pathology_tumour_fraction")
    paired_purity = numeric_values(paired, "ascat_purity")
    absolute_errors = [
        abs(float(row["ascat_purity"]) - float(row["pathology_tumour_fraction"]))
        for row in paired
    ]

    metric_rows = []
    for metric in (
        "ascat_ploidy",
        "pathology_tumour_fraction",
        "large_CNV_event_count",
        "cnv_burden_altered_segments",
        "segments",
        "SV_total_records",
        "tumour_pass_mutations",
        "tumour_wes_raw_total_reads",
        "normal_wes_raw_total_reads",
        "tumour_wes_reads_duplicated_percent",
        "normal_wes_reads_duplicated_percent",
    ):
        metric_rows.append(
            {
                "metric": metric,
                "boundary_n": len(numeric_values(boundary, metric)),
                "boundary_median": median(numeric_values(boundary, metric)),
                "interior_n": len(numeric_values(interior, metric)),
                "interior_median": median(numeric_values(interior, metric)),
            }
        )

    metric_output = args.output_dir / "ascat_boundary_group_comparison.tsv"
    with metric_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metric_rows[0]), delimiter="\t")
        writer.writeheader()
        for row in metric_rows:
            writer.writerow(
                {
                    key: "" if value is None else fmt(value) if isinstance(value, float) else value
                    for key, value in row.items()
                }
            )

    summary = args.output_dir / "ascat_purity_audit_summary.md"
    exact_near_diploid = sum(
        bool(row["exact_upper_boundary"] and row["near_diploid"]) for row in audited
    )
    with summary.open("w", encoding="utf-8") as handle:
        handle.write("# ASCAT purity audit\n\n")
        handle.write("## Retained evidence\n\n")
        handle.write(f"- Patients with a retained ASCAT purity: {len(purities)}\n")
        handle.write(
            f"- Exact upper-boundary estimates (purity >= 0.9995): {len(boundary)}\n"
        )
        handle.write(f"- Distinct retained purity values: {len(set(purities))}\n")
        handle.write(
            f"- Boundary estimates with near-diploid ploidy (1.90-2.10): "
            f"{exact_near_diploid}/{len(boundary)}\n"
        )
        handle.write(
            f"- Patients with both ASCAT and pathology tumour fraction: {len(paired)}\n"
        )
        handle.write(
            f"- Spearman correlation, ASCAT versus pathology fraction: "
            f"{fmt(spearman(paired_purity, pathology))}\n"
        )
        handle.write(
            f"- Median absolute difference from pathology fraction: "
            f"{fmt(median(absolute_errors))}\n\n"
        )
        handle.write("## Interpretation rules\n\n")
        handle.write(
            "An exact 1.00 estimate is a model-boundary solution, not proof that a "
            "sample contains 100% tumour cells. Repeated boundary values, especially "
            "with near-diploid solutions or disagreement with pathology, indicate weak "
            "identifiability or an unsuitable optimum and require LogR/BAF profile review.\n\n"
        )
        handle.write(
            "The retained `mosdepth_mean_depth_genomewide_recal` field must not be used "
            "as exome target coverage: it was calculated over a genome-wide denominator. "
            "On-target depth and uniformity require the capture BED and original alignment "
            "or coverage reports.\n\n"
        )
        handle.write("## Evidence needed for resolution\n\n")
        handle.write(
            "1. Original ASCAT tumour/normal LogR and BAF files, segments, "
            "purity/ploidy table, sunrise/fit plots and process logs.\n"
        )
        handle.write(
            "2. The exact Sarek samplesheet, sex assignments, capture BED and ASCAT "
            "reference bundles used in the completed run.\n"
        )
        handle.write(
            "3. On-target WES depth/uniformity, contamination and tumour-normal "
            "identity checks.\n"
        )
        handle.write(
            "4. Independent purity/ploidy estimation from the original tumour-normal "
            "alignments, for example FACETS or PureCN, followed by concordance review.\n"
        )

    print(summary)
    print(patient_output)
    print(metric_output)


if __name__ == "__main__":
    main()
