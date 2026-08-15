#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
from collections import Counter
from pathlib import Path


ALT_BREAKEND_RE = re.compile(r"[\[\]]([^:\[\]]+):([0-9]+)[\[\]]")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize interchromosomal translocations from Manta somatic SV VCFs."
    )
    parser.add_argument("--manta-root", required=True, help="Root directory containing Manta somatic SV VCFs")
    parser.add_argument("--output-dir", required=True, help="Directory where summary TSVs will be written")
    return parser.parse_args()


def strip_suffix(name: str) -> str:
    suffix = ".manta.somatic_sv.vcf.gz"
    return name[: -len(suffix)] if name.endswith(suffix) else name


def parse_info(info_str: str) -> dict[str, str]:
    result = {}
    for field in info_str.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = value
        else:
            result[field] = "True"
    return result


def parse_alt_breakend(alt: str):
    match = ALT_BREAKEND_RE.search(alt)
    if not match:
        return None, None
    chrom, pos = match.group(1), int(match.group(2))
    return chrom, pos


def parse_sample_counts(format_str: str, sample_str: str) -> dict[str, str]:
    keys = format_str.split(":")
    vals = sample_str.split(":")
    return {k: v for k, v in zip(keys, vals)}


def main():
    args = parse_args()
    manta_root = Path(args.manta_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    vcf_paths = sorted(manta_root.rglob("*.manta.somatic_sv.vcf.gz"))
    if not vcf_paths:
        raise SystemExit(f"No Manta somatic SV VCFs found under {manta_root}")

    event_rows = []
    sample_counts = []
    partner_counts = Counter()

    for vcf_path in vcf_paths:
        sample_pair = strip_suffix(vcf_path.name)
        seen_keys = set()
        total_interchrom = 0

        with gzip.open(vcf_path, "rt") as fh:
            header_samples = []
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header_samples = line.rstrip().split("\t")[9:]
                    continue

                chrom, pos, record_id, ref, alt, qual, filt, info_str, fmt, *samples = line.rstrip().split("\t")
                info = parse_info(info_str)
                if info.get("SVTYPE") != "BND":
                    continue

                mate_chrom, mate_pos = parse_alt_breakend(alt)
                if not mate_chrom or mate_chrom == chrom:
                    continue

                mate_id = info.get("MATEID", "")
                dedup_key = tuple(sorted([record_id, mate_id])) if mate_id else (record_id,)
                if dedup_key in seen_keys:
                    continue
                seen_keys.add(dedup_key)

                total_interchrom += 1
                chrom_pair = "--".join(sorted([chrom, mate_chrom]))
                partner_counts[(sample_pair, chrom_pair)] += 1

                tumor_counts = {}
                if samples:
                    tumor_sample = samples[-1]
                    tumor_counts = parse_sample_counts(fmt, tumor_sample)

                event_rows.append(
                    {
                        "sample_pair": sample_pair,
                        "chrom1": chrom,
                        "pos1": pos,
                        "chrom2": mate_chrom,
                        "pos2": mate_pos,
                        "chrom_pair": chrom_pair,
                        "record_id": record_id,
                        "mate_id": mate_id,
                        "filter": filt,
                        "somatic_score": info.get("SOMATICSCORE", ""),
                        "tumor_pr": tumor_counts.get("PR", ""),
                        "tumor_sr": tumor_counts.get("SR", ""),
                    }
                )

        sample_counts.append(
            {
                "sample_pair": sample_pair,
                "interchromosomal_translocation_events": total_interchrom,
            }
        )

    partner_rows = [
        {
            "sample_pair": sample_pair,
            "chrom_pair": chrom_pair,
            "event_count": count,
        }
        for (sample_pair, chrom_pair), count in sorted(
            partner_counts.items(), key=lambda item: (item[0][0], -item[1], item[0][1])
        )
    ]

    cohort_partner_totals = Counter()
    cohort_partner_samples = Counter()
    for (sample_pair, chrom_pair), count in partner_counts.items():
        cohort_partner_totals[chrom_pair] += count
        cohort_partner_samples[chrom_pair] += 1

    cohort_partner_rows = [
        {
            "chrom_pair": chrom_pair,
            "total_events": count,
            "samples_with_event": cohort_partner_samples[chrom_pair],
        }
        for chrom_pair, count in sorted(cohort_partner_totals.items(), key=lambda item: (-item[1], item[0]))
    ]

    def write_tsv(path: Path, rows: list[dict]):
        if not rows:
            rows = []
        with path.open("w", newline="", encoding="utf-8") as fh:
            if rows:
                writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            else:
                fh.write("")

    write_tsv(output_dir / "manta_interchromosomal_translocations_by_sample.tsv", sample_counts)
    write_tsv(output_dir / "manta_interchromosomal_translocation_events.tsv", event_rows)
    write_tsv(output_dir / "manta_interchromosomal_translocation_partner_counts.tsv", partner_rows)
    write_tsv(output_dir / "manta_interchromosomal_translocation_partner_counts_cohort.tsv", cohort_partner_rows)


if __name__ == "__main__":
    main()
