#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


CONTEXTS = {
    "SBS96": "SBS96",
    "DBS": "DBS",
    "ID": "ID",
}


def read_activity_table(path: Path):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise SystemExit(f"No rows found in {path}")
    return rows, [field for field in reader.fieldnames if field != "Samples"]


def write_tsv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize SigProfilerAssignment activities into cohort-friendly tables."
    )
    parser.add_argument("--input-root", required=True, help="SigProfiler output root containing SBS96/DBS/ID")
    parser.add_argument("--output-dir", required=True, help="Directory where summary tables will be written")
    parser.add_argument("--top-n", type=int, default=3, help="Top signatures per sample/context to report")
    args = parser.parse_args()

    input_root = Path(args.input_root)
    output_dir = Path(args.output_dir)
    long_rows = []
    top_rows = []
    cohort_totals = []

    for context_name in CONTEXTS:
        activity_path = (
            input_root
            / context_name
            / "Assignment_Solution"
            / "Activities"
            / "Assignment_Solution_Activities.txt"
        )
        if not activity_path.is_file():
            raise SystemExit(f"Missing activity file: {activity_path}")

        rows, signatures = read_activity_table(activity_path)
        signature_totals = {sig: 0.0 for sig in signatures}
        signature_nonzero_samples = {sig: 0 for sig in signatures}

        for row in rows:
            sample = row["Samples"]
            values = []
            total_activity = 0.0
            for sig in signatures:
                value = float(row[sig])
                values.append((sig, value))
                total_activity += value
                signature_totals[sig] += value

            for sig, value in values:
                if value > 0:
                    signature_nonzero_samples[sig] += 1
                    long_rows.append(
                        {
                            "sample_pair": sample,
                            "context_group": context_name,
                            "signature": sig,
                            "activity": int(value) if value.is_integer() else value,
                            "fraction_within_context": (value / total_activity) if total_activity else 0.0,
                        }
                    )

            ranked = sorted(values, key=lambda item: item[1], reverse=True)
            for rank, (sig, value) in enumerate(ranked[: args.top_n], start=1):
                if value <= 0:
                    continue
                top_rows.append(
                    {
                        "sample_pair": sample,
                        "context_group": context_name,
                        "rank": rank,
                        "signature": sig,
                        "activity": int(value) if value.is_integer() else value,
                        "fraction_within_context": (value / total_activity) if total_activity else 0.0,
                        "total_context_activity": int(total_activity) if total_activity.is_integer() else total_activity,
                    }
                )

        for sig in signatures:
            total = signature_totals[sig]
            if total <= 0:
                continue
            cohort_totals.append(
                {
                    "context_group": context_name,
                    "signature": sig,
                    "total_activity": int(total) if total.is_integer() else total,
                    "samples_with_nonzero_activity": signature_nonzero_samples[sig],
                }
            )

    write_tsv(
        output_dir / "sigprofiler_assignment_activity_long.tsv",
        long_rows,
        ["sample_pair", "context_group", "signature", "activity", "fraction_within_context"],
    )
    write_tsv(
        output_dir / "sigprofiler_assignment_top_signatures.tsv",
        top_rows,
        [
            "sample_pair",
            "context_group",
            "rank",
            "signature",
            "activity",
            "fraction_within_context",
            "total_context_activity",
        ],
    )
    write_tsv(
        output_dir / "sigprofiler_assignment_cohort_totals.tsv",
        cohort_totals,
        ["context_group", "signature", "total_activity", "samples_with_nonzero_activity"],
    )


if __name__ == "__main__":
    main()
