#!/usr/bin/env python3
import argparse
import csv
import sys
from pathlib import Path

REQUIRED = {
    "sarek": [("patient", "sample"), ("lane",), ("fastq_1", "bam", "cram")],
    "rnaseq": [("sample",), ("strandedness",), ("fastq_1",)],
    "rnafusion": [("sample",), ("strandedness",), ("fastq_1", "bam", "cram", "junctions", "splice_junctions")],
}


def load_rows(path: Path):
    lines = [line for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError("No non-comment rows found.")
    return list(csv.DictReader(lines))


def check_required(rows, pipeline):
    header = set(rows[0].keys())
    problems = []
    for group in REQUIRED[pipeline]:
        if not any(col in header for col in group):
            problems.append(f"missing one of required columns: {', '.join(group)}")
    if pipeline == "sarek" and {"status", "sex"} - header:
        problems.append("recommended columns missing for tumor/normal setups: status, sex")
    return problems


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline", choices=sorted(REQUIRED), required=True)
    parser.add_argument("--input", required=True)
    args = parser.parse_args()

    path = Path(args.input)
    if not path.exists():
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        return 1
    try:
        rows = load_rows(path)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    if not rows:
        print("ERROR: samplesheet has no data rows", file=sys.stderr)
        return 1
    problems = check_required(rows, args.pipeline)
    if problems:
        for problem in problems:
            print(f"ERROR: {problem}", file=sys.stderr)
        return 1
    print(f"OK: {args.pipeline} samplesheet looks structurally valid with {len(rows)} data row(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
