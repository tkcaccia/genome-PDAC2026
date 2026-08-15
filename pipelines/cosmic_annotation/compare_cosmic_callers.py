#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def truthy(value: str) -> bool:
    return str(value).strip().lower() == "yes"


def summarize(rows, caller_name: str):
    sample_stats = defaultdict(
        lambda: {
            "pass_variant_count": 0,
            "cgc_hit_count": 0,
            "cosmic_mutant_hit_count": 0,
            "cosmic_resistance_hit_count": 0,
            "unique_genes": set(),
            "cosmic_match_genes": set(),
        }
    )

    caller_stats = {
        "caller": caller_name,
        "pass_variant_count": 0,
        "cgc_hit_count": 0,
        "cosmic_mutant_hit_count": 0,
        "cosmic_resistance_hit_count": 0,
        "unique_genes": set(),
        "cosmic_match_genes": set(),
        "sample_pairs": set(),
    }

    for row in rows:
        sample = row["sample_pair"]
        gene = row.get("gene_symbol", "")
        caller_stats["sample_pairs"].add(sample)
        caller_stats["pass_variant_count"] += 1
        sample_stats[sample]["pass_variant_count"] += 1

        if gene:
            caller_stats["unique_genes"].add(gene)
            sample_stats[sample]["unique_genes"].add(gene)

        if truthy(row.get("cgc_hit", "")):
            caller_stats["cgc_hit_count"] += 1
            sample_stats[sample]["cgc_hit_count"] += 1

        if truthy(row.get("cosmic_mutant_hit", "")):
            caller_stats["cosmic_mutant_hit_count"] += 1
            sample_stats[sample]["cosmic_mutant_hit_count"] += 1
            if gene:
                caller_stats["cosmic_match_genes"].add(gene)
                sample_stats[sample]["cosmic_match_genes"].add(gene)

        if truthy(row.get("cosmic_resistance_hit", "")):
            caller_stats["cosmic_resistance_hit_count"] += 1
            sample_stats[sample]["cosmic_resistance_hit_count"] += 1

    caller_row = {
        "caller": caller_name,
        "sample_pair_count": len(caller_stats["sample_pairs"]),
        "pass_variant_count": caller_stats["pass_variant_count"],
        "cgc_hit_count": caller_stats["cgc_hit_count"],
        "cosmic_mutant_hit_count": caller_stats["cosmic_mutant_hit_count"],
        "cosmic_resistance_hit_count": caller_stats["cosmic_resistance_hit_count"],
        "unique_gene_count": len(caller_stats["unique_genes"]),
        "cosmic_match_gene_count": len(caller_stats["cosmic_match_genes"]),
    }

    sample_rows = []
    for sample, stats in sorted(sample_stats.items()):
        sample_rows.append(
            {
                "caller": caller_name,
                "sample_pair": sample,
                "pass_variant_count": stats["pass_variant_count"],
                "cgc_hit_count": stats["cgc_hit_count"],
                "cosmic_mutant_hit_count": stats["cosmic_mutant_hit_count"],
                "cosmic_resistance_hit_count": stats["cosmic_resistance_hit_count"],
                "unique_gene_count": len(stats["unique_genes"]),
                "cosmic_match_gene_count": len(stats["cosmic_match_genes"]),
            }
        )

    return caller_row, sample_rows


def main():
    parser = argparse.ArgumentParser(
        description="Build a caller-level comparison from Mutect2 and Strelka COSMIC annotation tables."
    )
    parser.add_argument("--mutect2", required=True)
    parser.add_argument("--strelka", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mutect_rows = read_tsv(Path(args.mutect2))
    strelka_rows = read_tsv(Path(args.strelka))

    caller_rows = []
    sample_rows = []
    for caller_name, rows in (("Mutect2", mutect_rows), ("Strelka", strelka_rows)):
        caller_row, sample_level = summarize(rows, caller_name)
        caller_rows.append(caller_row)
        sample_rows.extend(sample_level)

    with open(output_dir / "cosmic_caller_level_summary.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "caller",
                "sample_pair_count",
                "pass_variant_count",
                "cgc_hit_count",
                "cosmic_mutant_hit_count",
                "cosmic_resistance_hit_count",
                "unique_gene_count",
                "cosmic_match_gene_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(caller_rows)

    with open(output_dir / "cosmic_sample_level_summary.tsv", "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "caller",
                "sample_pair",
                "pass_variant_count",
                "cgc_hit_count",
                "cosmic_mutant_hit_count",
                "cosmic_resistance_hit_count",
                "unique_gene_count",
                "cosmic_match_gene_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(sample_rows)


if __name__ == "__main__":
    main()
