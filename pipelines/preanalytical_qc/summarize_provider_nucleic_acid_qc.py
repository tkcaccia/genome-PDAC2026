#!/usr/bin/env python3
"""Extract and summarize the CPGR nucleic-acid QC report format.

The generated sample-level tables contain private sample identifiers and must be
written outside the public repository. The script reports extraction quality;
it does not estimate sequencing depth, library complexity, or tumour purity.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path
from typing import Any


DNA_FIELDS = [
    "sample_id",
    "concentration_ng_ul",
    "a260_280",
    "a260_230",
    "din",
    "quality_score",
    "provider_qc",
]
RNA_FIELDS = [
    "sample_id",
    "concentration_ng_ul",
    "a260_280",
    "a260_230",
    "rin",
    "dv200_percent",
    "provider_qc",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract aggregate DNA/RNA QC from the CPGR report format."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--pdf", type=Path)
    source.add_argument(
        "--dna-qc-tsv",
        type=Path,
        help="Previously extracted private DNA QC table.",
    )
    parser.add_argument(
        "--rna-qc-tsv",
        type=Path,
        help="Previously extracted private RNA QC table; required with --dna-qc-tsv.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--ascat-audit",
        type=Path,
        help="Optional private output from audit_ascat_purity.py.",
    )
    return parser.parse_args()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def number(value: Any) -> float | None:
    text = clean(value).replace(",", ".")
    if not text or text in {"---", "-"}:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def last_value(row: list[Any]) -> str:
    for value in reversed(row):
        text = clean(value)
        if text:
            return text
    return ""


def detect_table_type(table: list[list[Any]]) -> str | None:
    header = " ".join(clean(value) for row in table[:3] for value in row)
    if "DIN" in header and "Concentration" in header:
        return "DNA"
    if "RIN" in header and "DV200" in header:
        return "RNA"
    return None


def data_rows(table: list[list[Any]]) -> list[list[Any]]:
    return [row for row in table if clean(row[0]).isdigit()]


def parse_report(pdf_path: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    try:
        import pdfplumber
    except ImportError as error:
        raise RuntimeError(
            "PDF extraction requires pdfplumber. Install the pinned Python "
            "requirements before running this command."
        ) from error

    dna: list[dict[str, object]] = []
    rna: list[dict[str, object]] = []
    with pdfplumber.open(pdf_path) as report:
        for page in report.pages:
            for table in page.extract_tables():
                table_type = detect_table_type(table)
                if table_type is None:
                    continue
                for row in data_rows(table):
                    # The first DNA table spans more PDF extraction cells than
                    # the continued DNA and RNA tables.
                    if table_type == "DNA" and len(row) >= 29:
                        record = {
                            "sample_id": clean(row[7]),
                            "concentration_ng_ul": number(row[10]),
                            "a260_280": number(row[13]),
                            "a260_230": number(row[16]),
                            "din": number(row[22]),
                            "quality_score": number(row[25]),
                            "provider_qc": last_value(row).upper(),
                        }
                    elif table_type == "DNA" and len(row) >= 21:
                        record = {
                            "sample_id": clean(row[7]),
                            "concentration_ng_ul": number(row[9]),
                            "a260_280": number(row[12]),
                            "a260_230": number(row[14]),
                            "din": number(row[19]),
                            "quality_score": number(row[20]),
                            "provider_qc": last_value(row).upper(),
                        }
                    elif table_type == "RNA" and len(row) >= 21:
                        record = {
                            "sample_id": clean(row[7]),
                            "concentration_ng_ul": number(row[9]),
                            "a260_280": number(row[12]),
                            "a260_230": number(row[14]),
                            "rin": number(row[19]),
                            "dv200_percent": number(row[20]),
                            "provider_qc": last_value(row).upper(),
                        }
                    else:
                        continue
                    (dna if table_type == "DNA" else rna).append(record)

    validate_records(dna, rna)
    return dna, rna


def validate_records(
    dna: list[dict[str, object]], rna: list[dict[str, object]]
) -> None:
    for label, records in (("DNA", dna), ("RNA", rna)):
        sample_ids = [str(row["sample_id"]) for row in records]
        if len(records) != 28 or len(set(sample_ids)) != 28:
            raise ValueError(
                f"Expected 28 unique {label} records, found "
                f"{len(records)} rows and {len(set(sample_ids))} unique IDs."
            )
        invalid = sorted(
            str(row["provider_qc"])
            for row in records
            if row["provider_qc"] not in {"PASS", "FAIL"}
        )
        if invalid:
            raise ValueError(f"Unexpected {label} QC labels: {invalid}")


def values(records: list[dict[str, object]], field: str) -> list[float]:
    return [
        float(row[field])
        for row in records
        if isinstance(row.get(field), (int, float))
    ]


def median_range(records: list[dict[str, object]], field: str) -> str:
    observed = values(records, field)
    if not observed:
        return "NA"
    return (
        f"median {statistics.median(observed):.2f}; "
        f"range {min(observed):.2f}-{max(observed):.2f}; n={len(observed)}"
    )


def write_tsv(path: Path, records: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in records:
            writer.writerow({field: "" if row.get(field) is None else row.get(field) for field in fields})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def coerce_qc_rows(
    records: list[dict[str, str]], numeric_fields: list[str]
) -> list[dict[str, object]]:
    return [
        {
            **row,
            **{field: number(row.get(field)) for field in numeric_fields},
            "provider_qc": clean(row.get("provider_qc")).upper(),
        }
        for row in records
    ]


def ascat_comparison(
    dna: list[dict[str, object]], ascat_path: Path | None
) -> list[dict[str, object]]:
    if ascat_path is None:
        return []
    dna_by_sample = {str(row["sample_id"]): row for row in dna}
    result = []
    for row in read_tsv(ascat_path):
        sample = f"{row['patient_id']}T"
        dna_row = dna_by_sample.get(sample)
        if dna_row is None:
            raise ValueError(f"No provider DNA QC row found for ASCAT tumour {sample}")
        result.append(
            {
                **dna_row,
                "ascat_boundary": row.get("exact_upper_boundary", "").lower() == "true",
            }
        )
    return result


def comparison_summary(records: list[dict[str, object]]) -> list[dict[str, str]]:
    if not records:
        return []
    boundary = [row for row in records if row["ascat_boundary"]]
    interior = [row for row in records if not row["ascat_boundary"]]
    output = []
    for metric in ("concentration_ng_ul", "din", "a260_280", "a260_230"):
        output.append(
            {
                "metric": metric,
                "boundary": median_range(boundary, metric),
                "non_boundary": median_range(interior, metric),
            }
        )
    output.append(
        {
            "metric": "provider_qc_pass",
            "boundary": f"{sum(row['provider_qc'] == 'PASS' for row in boundary)}/{len(boundary)}",
            "non_boundary": f"{sum(row['provider_qc'] == 'PASS' for row in interior)}/{len(interior)}",
        }
    )
    return output


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.pdf:
        dna, rna = parse_report(args.pdf)
    else:
        if args.rna_qc_tsv is None:
            raise SystemExit("--rna-qc-tsv is required with --dna-qc-tsv")
        dna = coerce_qc_rows(
            read_tsv(args.dna_qc_tsv),
            ["concentration_ng_ul", "a260_280", "a260_230", "din", "quality_score"],
        )
        rna = coerce_qc_rows(
            read_tsv(args.rna_qc_tsv),
            ["concentration_ng_ul", "a260_280", "a260_230", "rin", "dv200_percent"],
        )
        validate_records(dna, rna)
    write_tsv(args.output_dir / "provider_dna_qc_private.tsv", dna, DNA_FIELDS)
    write_tsv(args.output_dir / "provider_rna_qc_private.tsv", rna, RNA_FIELDS)

    summary_rows = [
        {"material": "DNA", "metric": "provider_qc_pass", "result": f"{sum(row['provider_qc'] == 'PASS' for row in dna)}/28"},
        {"material": "DNA", "metric": "concentration_ng_ul", "result": median_range(dna, "concentration_ng_ul")},
        {"material": "DNA", "metric": "a260_280", "result": median_range(dna, "a260_280")},
        {"material": "DNA", "metric": "a260_230", "result": median_range(dna, "a260_230")},
        {"material": "DNA", "metric": "din", "result": median_range(dna, "din")},
        {"material": "RNA", "metric": "provider_qc_pass", "result": f"{sum(row['provider_qc'] == 'PASS' for row in rna)}/28"},
        {"material": "RNA", "metric": "concentration_ng_ul", "result": median_range(rna, "concentration_ng_ul")},
        {"material": "RNA", "metric": "a260_280", "result": median_range(rna, "a260_280")},
        {"material": "RNA", "metric": "a260_230", "result": median_range(rna, "a260_230")},
        {"material": "RNA", "metric": "rin", "result": median_range(rna, "rin")},
        {"material": "RNA", "metric": "dv200_percent", "result": median_range(rna, "dv200_percent")},
        {"material": "RNA", "metric": "dv200_at_least_30_percent", "result": f"{sum((row['dv200_percent'] or 0) >= 30 for row in rna)}/28"},
    ]
    write_tsv(
        args.output_dir / "provider_qc_aggregate_summary.tsv",
        summary_rows,
        ["material", "metric", "result"],
    )

    comparison = comparison_summary(ascat_comparison(dna, args.ascat_audit))
    if comparison:
        write_tsv(
            args.output_dir / "provider_dna_qc_by_ascat_boundary.tsv",
            comparison,
            ["metric", "boundary", "non_boundary"],
        )

    print(args.output_dir / "provider_qc_aggregate_summary.tsv")
    if comparison:
        print(args.output_dir / "provider_dna_qc_by_ascat_boundary.tsv")


if __name__ == "__main__":
    main()
