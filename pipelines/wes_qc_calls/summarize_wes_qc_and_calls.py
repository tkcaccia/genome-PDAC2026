#!/usr/bin/env python3
"""Summarize targeted WES depth and tumour-normal genomic call quality.

The script is designed for a provider DRAGEN enrichment-results directory plus
per-patient Hartwig Oncoanalyser output. It deliberately recalculates WES call
burden within regions that are callable in both tumour and matched normal. This
avoids treating whole-genome-style TMB fields as valid when Oncoanalyser was run
with ``targeted=false`` on exome data.

Only summary tables are written. Input BAM/CRAM/VCF files are never modified.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import math
import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence


SAMPLE_RE = re.compile(r"^(?P<patient>\d+)(?P<condition>[NT])WESFastq$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--enrichment-root",
        required=True,
        type=Path,
        help="Directory containing one DRAGEN enrichment result folder per WES sample.",
    )
    parser.add_argument(
        "--onco-root",
        required=True,
        type=Path,
        help="Directory containing PanC_<patient>/GRCh38_*_DNA_RNA_output results.",
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--min-tumour-depth",
        type=float,
        default=50.0,
        help="Descriptive review threshold for mean tumour target depth (default: 50x).",
    )
    parser.add_argument(
        "--min-normal-depth",
        type=float,
        default=30.0,
        help="Descriptive review threshold for mean normal target depth (default: 30x).",
    )
    parser.add_argument(
        "--min-pct-20x",
        type=float,
        default=80.0,
        help="Descriptive review threshold for target bases at >=20x (default: 80%%).",
    )
    parser.add_argument(
        "--min-pair-callable-pct",
        type=float,
        default=70.0,
        help="Descriptive review threshold for paired callable target territory (default: 70%%).",
    )
    parser.add_argument(
        "--outlier-robust-z",
        type=float,
        default=2.5,
        help="Descriptive log-burden robust-z threshold for a cohort outlier flag (default: 2.5).",
    )
    return parser.parse_args()


def clean_number(value: object) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NONE"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def format_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.6g}"
    return str(value)


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field)) for field in fields})


def read_two_column_csv(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            # DRAGEN summaries are comma-delimited, whereas PURPLE .qc files
            # use a two-column tab-delimited key/value format.
            row = line.rstrip("\n\r").split("\t") if "\t" in line else next(csv.reader([line]))
            if len(row) >= 2 and row[0].strip():
                result[row[0].strip()] = row[1].strip()
    return result


def read_metric_csv(path: Path) -> dict[str, tuple[str, str | None]]:
    result: dict[str, tuple[str, str | None]] = {}
    with path.open(newline="", errors="replace") as handle:
        for row in csv.reader(handle):
            if len(row) < 4:
                continue
            label = row[2].strip()
            if not label:
                continue
            value = row[3].strip()
            percent = row[4].strip() if len(row) > 4 and row[4].strip() else None
            result[label] = (value, percent)
    return result


def metric_number(metrics: Mapping[str, tuple[str, str | None]], label: str, percent: bool = False) -> float | None:
    value = metrics.get(label)
    if value is None:
        return None
    return clean_number(value[1] if percent else value[0])


@dataclass(frozen=True)
class SampleInput:
    patient_id: str
    condition: str
    source_sample: str
    sample_id: str
    result_dir: Path
    summary_csv: Path
    target_metrics_csv: Path
    mapping_metrics_csv: Path
    variant_metrics_csv: Path
    callable_bed: Path
    target_bed: Path
    bam: Path
    hard_filtered_vcf: Path


def discover_enrichment_samples(root: Path) -> list[SampleInput]:
    samples: list[SampleInput] = []
    for summary_path in sorted(root.glob("*/*/Additional Files/*.summary.csv")):
        source_sample = summary_path.name.removesuffix(".summary.csv")
        match = SAMPLE_RE.match(source_sample)
        if not match:
            continue
        result_dir = summary_path.parent.parent
        patient = str(int(match.group("patient")))
        condition = match.group("condition")
        samples.append(
            SampleInput(
                patient_id=patient,
                condition=condition,
                source_sample=source_sample,
                sample_id=f"{patient}{condition}_WES",
                result_dir=result_dir,
                summary_csv=summary_path,
                target_metrics_csv=result_dir / f"{source_sample}.target_bed_coverage_metrics.csv",
                mapping_metrics_csv=result_dir / f"{source_sample}.mapping_metrics.csv",
                variant_metrics_csv=result_dir / f"{source_sample}.vc_metrics.csv",
                callable_bed=result_dir / f"{source_sample}.target_bed_somatic_callable_regions.bed",
                target_bed=result_dir / "Twist_ILMN_Exome_2.5_Plus_Panel.hg38.bed",
                bam=result_dir / f"{source_sample}.bam",
                hard_filtered_vcf=result_dir / f"{source_sample}.hard-filtered.vcf.gz",
            )
        )
    return samples


def bed_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 3:
                fields = line.split()
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if end > start:
                intervals[fields[0]].append((start, end))
    return {chrom: merge_intervals(values) for chrom, values in intervals.items()}


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted(intervals)
    if not ordered:
        return []
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        previous_start, previous_end = merged[-1]
        if start <= previous_end:
            merged[-1] = (previous_start, max(previous_end, end))
        else:
            merged.append((start, end))
    return merged


def interval_bases(intervals: Mapping[str, Sequence[tuple[int, int]]]) -> int:
    return sum(end - start for values in intervals.values() for start, end in values)


def intersect_interval_sets(
    left: Mapping[str, Sequence[tuple[int, int]]],
    right: Mapping[str, Sequence[tuple[int, int]]],
) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = {}
    for chrom in sorted(set(left) & set(right)):
        a = left[chrom]
        b = right[chrom]
        i = j = 0
        overlaps: list[tuple[int, int]] = []
        while i < len(a) and j < len(b):
            start = max(a[i][0], b[j][0])
            end = min(a[i][1], b[j][1])
            if start < end:
                overlaps.append((start, end))
            if a[i][1] <= b[j][1]:
                i += 1
            else:
                j += 1
        if overlaps:
            result[chrom] = merge_intervals(overlaps)
    return result


class IntervalIndex:
    def __init__(self, intervals: Mapping[str, Sequence[tuple[int, int]]]):
        self.intervals = {chrom: list(values) for chrom, values in intervals.items()}
        self.starts = {chrom: [start for start, _ in values] for chrom, values in self.intervals.items()}

    def contains(self, chrom: str, zero_based_position: int) -> bool:
        values = self.intervals.get(chrom)
        starts = self.starts.get(chrom)
        if not values or not starts:
            return False
        index = bisect.bisect_right(starts, zero_based_position) - 1
        return index >= 0 and zero_based_position < values[index][1]


def count_bed_bases(path: Path) -> int | None:
    if not path.exists():
        return None
    return interval_bases(bed_intervals(path))


def collect_sample_qc(sample: SampleInput) -> dict[str, object]:
    summary = read_two_column_csv(sample.summary_csv)
    target = read_metric_csv(sample.target_metrics_csv)
    mapping = read_metric_csv(sample.mapping_metrics_csv)

    def summary_num(name: str) -> float | None:
        return clean_number(summary.get(name))

    return {
        "patient_id": sample.patient_id,
        "sample_id": sample.sample_id,
        "condition": "tumour" if sample.condition == "T" else "normal",
        "reference_genome": summary.get("Reference Genome", "NA"),
        "target_manifest": summary.get("Target manifest", "NA"),
        "target_bases": summary_num("Total length of targeted reference"),
        "target_regions": summary_num("Number of targeted regions"),
        "total_pf_reads": summary_num("Total PF reads"),
        "mapped_reads_pct": summary_num("Percent aligned reads"),
        "unique_mapped_reads_pct": summary_num("Percent unique aligned reads"),
        "duplicate_aligned_reads_pct": summary_num("Percent duplicate aligned reads"),
        "q30_bases_pct": summary_num("Percent Q30 bases"),
        "mismatch_bases_pct": summary_num("Percent mismatches"),
        "mean_target_depth_x": summary_num("Mean target coverage depth"),
        "quality_adjusted_target_depth_x": metric_number(target, "Average alignment coverage over target region"),
        "target_coverage_ge_1x_pct": summary_num("Percent Target coverage at 1X"),
        "target_coverage_ge_10x_pct": summary_num("Percent Target coverage at 10X"),
        "target_coverage_ge_15x_pct": summary_num("Percent Target coverage at 15X"),
        "target_coverage_ge_20x_pct": summary_num("Percent Target coverage at 20X"),
        "target_coverage_ge_50x_pct": summary_num("Percent Target coverage at 50X"),
        "target_coverage_ge_100x_pct": metric_number(
            target, "PCT of target region with coverage [ 100x: inf)"
        ),
        "target_uniformity_0_2_mean_pct": summary_num("Uniformity of coverage (Pct > 0.2*mean)"),
        "unique_read_enrichment_pct": summary_num("Unique read enrichment"),
        "unique_base_enrichment_pct": summary_num("Unique base enrichment"),
        "aligned_reads_in_target_pct": metric_number(target, "Aligned reads in target region", percent=True),
        "provider_callability_pct": summary_num("Percent callability"),
        "provider_autosome_callability_pct": summary_num("Percent autosome callability"),
        "somatic_callable_bases": count_bed_bases(sample.callable_bed),
        "estimated_chr_x_target_depth_x": metric_number(target, "Average chr X coverage over target region"),
        "estimated_chr_y_target_depth_x": metric_number(target, "Average chr Y coverage over target region"),
        "fragment_length_median_bp": summary_num("Fragment length median"),
        "dragen_germline_snvs": summary_num("SNVs"),
        "dragen_germline_indels": summary_num("Indels"),
        "dragen_germline_titv": summary_num("SNV Ts/Tv ratio"),
        "bam_present": sample.bam.exists(),
        "bam_bytes": sample.bam.stat().st_size if sample.bam.exists() else None,
        "hard_filtered_vcf_present": sample.hard_filtered_vcf.exists(),
        "source_summary_csv": str(sample.summary_csv),
        "source_callable_bed": str(sample.callable_bed),
    }


def parse_info(text: str) -> dict[str, str | bool]:
    info: dict[str, str | bool] = {}
    for item in text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def parse_format(format_text: str, sample_text: str) -> dict[str, str]:
    keys = format_text.split(":")
    values = sample_text.split(":")
    return dict(zip(keys, values))


def alt_depth(sample: Mapping[str, str]) -> int | None:
    ad = sample.get("AD", "")
    if not ad or ad == ".":
        return None
    values = ad.split(",")
    if len(values) < 2:
        return None
    try:
        return sum(int(value) for value in values[1:] if value != ".")
    except ValueError:
        return None


def variant_type(ref: str, alt: str) -> str:
    first_alt = alt.split(",", 1)[0]
    if len(ref) == len(first_alt) == 1:
        return "SNV"
    if len(ref) == len(first_alt):
        return "MNV"
    if len(ref) < len(first_alt):
        return "INS"
    return "DEL"


def impact_fields(info: Mapping[str, str | bool]) -> tuple[str, str, str, str]:
    impact = info.get("IMPACT")
    if not isinstance(impact, str):
        return "NA", "NA", "NA", "NA"
    fields = impact.split(",")
    gene = fields[0] if fields else "NA"
    transcript = fields[1] if len(fields) > 1 else "NA"
    consequence = fields[2] if len(fields) > 2 else "NA"
    coding_effect = fields[3] if len(fields) > 3 else "NA"
    return gene, transcript, consequence, coding_effect


def is_nonsynonymous(consequence: str, coding_effect: str) -> bool:
    coding = coding_effect.upper()
    if coding not in {"", "NA", "NONE", "SYNONYMOUS", "UNDEFINED"}:
        return True
    consequence_lower = consequence.lower()
    return any(
        token in consequence_lower
        for token in ("missense", "frameshift", "stop_gained", "stop_lost", "start_lost", "splice")
    )


def median_or_none(values: Sequence[float]) -> float | None:
    return statistics.median(values) if values else None


def read_vcf_records(path: Path) -> Iterator[tuple[list[str], list[str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    sample_names: list[str] = []
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n\r").split("\t")
                sample_names = header[9:]
                continue
            if not line.strip() or line.startswith("#"):
                continue
            yield line.rstrip("\n\r").split("\t"), sample_names


def summarize_vcf(
    path: Path,
    patient_id: str,
    caller: str,
    callable_index: IntervalIndex,
    collect_calls: bool,
) -> tuple[dict[str, object], Counter[str], list[dict[str, object]]]:
    counts = Counter()
    filters: Counter[str] = Counter()
    calls: list[dict[str, object]] = []
    tumour_depths: list[float] = []
    normal_depths: list[float] = []
    tumour_vafs: list[float] = []
    normal_vafs: list[float] = []

    for fields, sample_names in read_vcf_records(path):
        if len(fields) < 8:
            continue
        chrom, pos, _variant_id, ref, alt, qual, filter_value, info_text = fields[:8]
        try:
            position = int(pos)
        except ValueError:
            continue
        info = parse_info(info_text)
        in_callable = callable_index.contains(chrom, position - 1)
        passed = filter_value in {"PASS", "."}
        kind = variant_type(ref, alt)
        gene, transcript, consequence, coding_effect = impact_fields(info)
        nonsynonymous = is_nonsynonymous(consequence, coding_effect)

        counts["records"] += 1
        counts[f"type_{kind}"] += 1
        if passed:
            counts["pass"] += 1
        if in_callable:
            counts["pair_callable"] += 1
        if passed and in_callable:
            counts["pass_pair_callable"] += 1
            counts[f"pass_pair_callable_{kind}"] += 1
            if nonsynonymous:
                counts["pass_pair_callable_nonsynonymous"] += 1
            if info.get("TIER") == "HOTSPOT":
                counts["pass_pair_callable_hotspot"] += 1

        for filter_name in filter_value.split(";"):
            filters[filter_name or "NA"] += 1

        normal_name = next((name for name in sample_names if name.endswith("-N")), None)
        tumour_name = next((name for name in sample_names if name.endswith("-T")), None)
        sample_columns = dict(zip(sample_names, fields[9:]))
        format_text = fields[8] if len(fields) > 8 else ""
        normal = parse_format(format_text, sample_columns.get(normal_name, "")) if normal_name else {}
        tumour = parse_format(format_text, sample_columns.get(tumour_name, "")) if tumour_name else {}
        normal_dp = clean_number(normal.get("DP"))
        tumour_dp = clean_number(tumour.get("DP"))
        normal_af = clean_number(normal.get("AF"))
        tumour_af = clean_number(tumour.get("AF"))
        normal_alt = alt_depth(normal)
        tumour_alt = alt_depth(tumour)
        strict_evidence = bool(
            passed
            and in_callable
            and tumour_dp is not None
            and tumour_dp >= 20
            and normal_dp is not None
            and normal_dp >= 10
            and tumour_alt is not None
            and tumour_alt >= 3
            and tumour_af is not None
            and tumour_af >= 0.05
            and normal_alt is not None
            and normal_alt <= 1
            and normal_af is not None
            and normal_af <= 0.02
        )

        if passed and in_callable:
            if tumour_dp is not None:
                tumour_depths.append(tumour_dp)
            if normal_dp is not None:
                normal_depths.append(normal_dp)
            if tumour_af is not None:
                tumour_vafs.append(tumour_af)
            if normal_af is not None:
                normal_vafs.append(normal_af)
            if strict_evidence:
                counts["strict_evidence"] += 1
                counts[f"strict_evidence_{kind}"] += 1
                if nonsynonymous:
                    counts["strict_evidence_nonsynonymous"] += 1
                if info.get("TIER") == "HOTSPOT":
                    counts["strict_evidence_hotspot"] += 1

            if collect_calls:
                calls.append(
                    {
                        "patient_id": patient_id,
                        "caller": caller,
                        "chrom": chrom,
                        "position": position,
                        "ref": ref,
                        "alt": alt,
                        "variant_type": kind,
                        "filter": filter_value,
                        "qual": clean_number(qual),
                        "tier": info.get("TIER", "NA"),
                        "gene": gene,
                        "transcript": transcript,
                        "consequence": consequence,
                        "coding_effect": coding_effect,
                        "nonsynonymous": nonsynonymous,
                        "tumour_depth": tumour_dp,
                        "tumour_alt_depth": tumour_alt,
                        "tumour_vaf": tumour_af,
                        "normal_depth": normal_dp,
                        "normal_alt_depth": normal_alt,
                        "normal_vaf": normal_af,
                        "strict_read_evidence": strict_evidence,
                        "purple_copy_number": clean_number(info.get("PURPLE_CN")),
                        "purple_variant_copy_number": clean_number(info.get("PURPLE_VCN")),
                        "biallelic_flag": "BIALLELIC" in info,
                        "source_vcf": str(path),
                    }
                )

    summary = {
        "patient_id": patient_id,
        "caller": caller,
        "vcf_records": counts["records"],
        "vcf_pass_records": counts["pass"],
        "pair_callable_records": counts["pair_callable"],
        "pass_pair_callable_records": counts["pass_pair_callable"],
        "pass_pair_callable_snvs": counts["pass_pair_callable_SNV"],
        "pass_pair_callable_mnvs": counts["pass_pair_callable_MNV"],
        "pass_pair_callable_insertions": counts["pass_pair_callable_INS"],
        "pass_pair_callable_deletions": counts["pass_pair_callable_DEL"],
        "pass_pair_callable_nonsynonymous": counts["pass_pair_callable_nonsynonymous"],
        "pass_pair_callable_hotspots": counts["pass_pair_callable_hotspot"],
        "strict_read_evidence_records": counts["strict_evidence"],
        "strict_read_evidence_snvs": counts["strict_evidence_SNV"],
        "strict_read_evidence_mnvs": counts["strict_evidence_MNV"],
        "strict_read_evidence_insertions": counts["strict_evidence_INS"],
        "strict_read_evidence_deletions": counts["strict_evidence_DEL"],
        "strict_read_evidence_nonsynonymous": counts["strict_evidence_nonsynonymous"],
        "strict_read_evidence_hotspots": counts["strict_evidence_hotspot"],
        "median_tumour_depth_at_pass_calls": median_or_none(tumour_depths),
        "median_normal_depth_at_pass_calls": median_or_none(normal_depths),
        "median_tumour_vaf_at_pass_calls": median_or_none(tumour_vafs),
        "median_normal_vaf_at_pass_calls": median_or_none(normal_vafs),
        "source_vcf": str(path),
    }
    return summary, filters, calls


def read_single_tsv_row(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return next(reader, {})


def driver_catalog_counts(path: Path) -> tuple[int, int, int, int]:
    if not path.exists():
        return 0, 0, 0, 0
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    likelihood_ge_0_5 = sum(
        1 for row in rows if (clean_number(row.get("driverLikelihood")) or 0.0) >= 0.5
    )
    mutation_rows = sum("MUTATION" in str(row.get("driver", "")).upper() for row in rows)
    deletion_rows = sum("DELETION" in str(row.get("driver", "")).upper() for row in rows)
    return len(rows), likelihood_ge_0_5, mutation_rows, deletion_rows


def onco_patient_dir(root: Path, patient_id: str) -> Path:
    padded = f"{int(patient_id):02d}"
    return root / f"PanC_{int(patient_id)}" / f"GRCh38_PanC_{int(patient_id)}_DNA_RNA_output" / f"PATIENT{padded}"


def main() -> int:
    args = parse_args()
    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    samples = discover_enrichment_samples(args.enrichment_root)
    if not samples:
        raise SystemExit(f"No DRAGEN enrichment sample summaries found below {args.enrichment_root}")

    sample_inputs = {(sample.patient_id, sample.condition): sample for sample in samples}
    sample_qc = [collect_sample_qc(sample) for sample in samples]
    sample_qc_lookup = {
        (str(row["patient_id"]), "T" if row["condition"] == "tumour" else "N"): row for row in sample_qc
    }

    pair_rows: list[dict[str, object]] = []
    call_summaries: list[dict[str, object]] = []
    filter_rows: list[dict[str, object]] = []
    somatic_calls: list[dict[str, object]] = []
    germline_rows: list[dict[str, object]] = []
    input_rows: list[dict[str, object]] = []

    patient_ids = sorted({sample.patient_id for sample in samples}, key=int)
    for patient_id in patient_ids:
        tumour_input = sample_inputs.get((patient_id, "T"))
        normal_input = sample_inputs.get((patient_id, "N"))
        tumour_qc = sample_qc_lookup.get((patient_id, "T"), {})
        normal_qc = sample_qc_lookup.get((patient_id, "N"), {})
        onco_dir = onco_patient_dir(args.onco_root, patient_id)
        padded = f"{int(patient_id):02d}"
        prefix = f"PATIENT{padded}-T"
        purple_vcf = onco_dir / "purple" / f"{prefix}.purple.somatic.vcf.gz"
        sage_vcf = onco_dir / "sage" / "somatic" / f"{prefix}.sage.somatic.vcf.gz"
        purity_path = onco_dir / "purple" / f"{prefix}.purple.purity.tsv"
        purple_qc_path = onco_dir / "purple" / f"{prefix}.purple.qc"
        somatic_driver_path = onco_dir / "purple" / f"{prefix}.purple.driver.catalog.somatic.tsv"
        germline_driver_path = onco_dir / "purple" / f"{prefix}.purple.driver.catalog.germline.tsv"
        germline_vcf = onco_dir / "pave" / f"{prefix}.sage.germline.pave.vcf.gz"

        pair_callable: dict[str, list[tuple[int, int]]] = {}
        pair_callable_bases = 0
        if tumour_input and normal_input and tumour_input.callable_bed.exists() and normal_input.callable_bed.exists():
            pair_callable = intersect_interval_sets(
                bed_intervals(tumour_input.callable_bed), bed_intervals(normal_input.callable_bed)
            )
            pair_callable_bases = interval_bases(pair_callable)
        callable_index = IntervalIndex(pair_callable)
        target_bases = clean_number(tumour_qc.get("target_bases")) or clean_number(normal_qc.get("target_bases"))
        pair_callable_pct = (
            100.0 * pair_callable_bases / target_bases if target_bases and pair_callable_bases else None
        )

        reasons: list[str] = []
        tumour_depth = clean_number(tumour_qc.get("mean_target_depth_x"))
        normal_depth = clean_number(normal_qc.get("mean_target_depth_x"))
        tumour_20x = clean_number(tumour_qc.get("target_coverage_ge_20x_pct"))
        normal_20x = clean_number(normal_qc.get("target_coverage_ge_20x_pct"))
        if tumour_depth is None or tumour_depth < args.min_tumour_depth:
            reasons.append("tumour_mean_depth_below_review_threshold")
        if normal_depth is None or normal_depth < args.min_normal_depth:
            reasons.append("normal_mean_depth_below_review_threshold")
        if tumour_20x is None or tumour_20x < args.min_pct_20x:
            reasons.append("tumour_20x_coverage_below_review_threshold")
        if normal_20x is None or normal_20x < args.min_pct_20x:
            reasons.append("normal_20x_coverage_below_review_threshold")
        if pair_callable_pct is None or pair_callable_pct < args.min_pair_callable_pct:
            reasons.append("paired_callable_territory_below_review_threshold")

        pair_row: dict[str, object] = {
            "patient_id": patient_id,
            "tumour_sample_id": tumour_qc.get("sample_id"),
            "normal_sample_id": normal_qc.get("sample_id"),
            "tumour_mean_target_depth_x": tumour_depth,
            "normal_mean_target_depth_x": normal_depth,
            "tumour_target_20x_pct": tumour_20x,
            "normal_target_20x_pct": normal_20x,
            "tumour_duplicate_pct": tumour_qc.get("duplicate_aligned_reads_pct"),
            "normal_duplicate_pct": normal_qc.get("duplicate_aligned_reads_pct"),
            "tumour_provider_callability_pct": tumour_qc.get("provider_callability_pct"),
            "normal_provider_callability_pct": normal_qc.get("provider_callability_pct"),
            "pair_callable_bases": pair_callable_bases,
            "pair_callable_mb": pair_callable_bases / 1_000_000.0 if pair_callable_bases else None,
            "pair_callable_pct_of_target": pair_callable_pct,
            "descriptive_qc_review": "REVIEW" if reasons else "MEETS_THRESHOLDS",
            "review_reasons": ";".join(reasons) if reasons else "none",
        }

        purity = read_single_tsv_row(purity_path)
        purple_qc = read_two_column_csv(purple_qc_path) if purple_qc_path.exists() else {}
        pair_row.update(
            {
                "purple_purity": clean_number(purity.get("purity")),
                "purple_status": purity.get("status", "NA"),
                "purple_qc_status": purple_qc.get("QCStatus", "NA"),
                "onco_run_mode": purity.get("runMode", "NA"),
                "onco_targeted": purity.get("targeted", "NA"),
                "onco_reported_tmb_per_mb_not_used": clean_number(purity.get("tmbPerMb")),
                "onco_reported_tml": clean_number(purity.get("tml")),
            }
        )
        pair_rows.append(pair_row)

        for caller, vcf_path, collect_calls in (
            ("SAGE", sage_vcf, False),
            ("PURPLE", purple_vcf, True),
        ):
            if not vcf_path.exists() or pair_callable_bases == 0:
                continue
            summary, filters, calls = summarize_vcf(
                vcf_path, patient_id, caller, callable_index, collect_calls=collect_calls
            )
            summary["pair_callable_bases"] = pair_callable_bases
            summary["pair_callable_mb"] = pair_callable_bases / 1_000_000.0
            summary["pass_variants_per_pair_callable_mb"] = (
                float(summary["pass_pair_callable_records"]) / (pair_callable_bases / 1_000_000.0)
            )
            summary["pass_nonsynonymous_per_pair_callable_mb"] = (
                float(summary["pass_pair_callable_nonsynonymous"]) / (pair_callable_bases / 1_000_000.0)
            )
            summary["strict_read_evidence_per_pair_callable_mb"] = (
                float(summary["strict_read_evidence_records"]) / (pair_callable_bases / 1_000_000.0)
            )
            summary["strict_nonsynonymous_per_pair_callable_mb"] = (
                float(summary["strict_read_evidence_nonsynonymous"]) / (pair_callable_bases / 1_000_000.0)
            )
            summary["onco_targeted_mode"] = purity.get("targeted", "NA")
            summary["clinical_tmb_interpretation"] = (
                "NOT_CALCULATED: full paired callable target used, not validated coding territory"
            )
            call_summaries.append(summary)
            for filter_name, count in sorted(filters.items()):
                filter_rows.append(
                    {
                        "patient_id": patient_id,
                        "caller": caller,
                        "filter": filter_name,
                        "record_count": count,
                    }
                )
            somatic_calls.extend(calls)

        somatic_drivers, somatic_likely, somatic_mutations, somatic_deletions = driver_catalog_counts(
            somatic_driver_path
        )
        germline_drivers, germline_likely, germline_mutations, germline_deletions = driver_catalog_counts(
            germline_driver_path
        )
        germline_summary: dict[str, object] = {
            "patient_id": patient_id,
            "normal_sample_id": normal_qc.get("sample_id"),
            "dragen_normal_germline_snvs": normal_qc.get("dragen_germline_snvs"),
            "dragen_normal_germline_indels": normal_qc.get("dragen_germline_indels"),
            "dragen_normal_germline_titv": normal_qc.get("dragen_germline_titv"),
            "purple_somatic_driver_catalog_rows": somatic_drivers,
            "purple_somatic_driver_likelihood_ge_0_5_rows": somatic_likely,
            "purple_somatic_driver_mutation_rows": somatic_mutations,
            "purple_somatic_driver_deletion_rows": somatic_deletions,
            "purple_germline_driver_catalog_rows": germline_drivers,
            "purple_germline_driver_likelihood_ge_0_5_rows": germline_likely,
            "purple_germline_driver_mutation_rows": germline_mutations,
            "purple_germline_driver_deletion_rows": germline_deletions,
            "germline_candidate_vcf_present": germline_vcf.exists(),
            "interpretation": "research candidates; clinical germline confirmation required",
        }
        germline_rows.append(germline_summary)

        input_rows.append(
            {
                "patient_id": patient_id,
                "tumour_enrichment_metrics": bool(tumour_input),
                "normal_enrichment_metrics": bool(normal_input),
                "tumour_bam": bool(tumour_input and tumour_input.bam.exists()),
                "normal_bam": bool(normal_input and normal_input.bam.exists()),
                "sage_somatic_vcf": sage_vcf.exists(),
                "purple_somatic_vcf": purple_vcf.exists(),
                "purple_purity": purity_path.exists(),
                "somatic_driver_catalog": somatic_driver_path.exists(),
                "germline_driver_catalog": germline_driver_path.exists(),
            }
        )

    sample_fields = [
        "patient_id", "sample_id", "condition", "reference_genome", "target_manifest",
        "target_bases", "target_regions", "total_pf_reads", "mapped_reads_pct",
        "unique_mapped_reads_pct", "duplicate_aligned_reads_pct", "q30_bases_pct",
        "mismatch_bases_pct", "mean_target_depth_x", "quality_adjusted_target_depth_x",
        "target_coverage_ge_1x_pct", "target_coverage_ge_10x_pct", "target_coverage_ge_15x_pct",
        "target_coverage_ge_20x_pct", "target_coverage_ge_50x_pct", "target_coverage_ge_100x_pct",
        "target_uniformity_0_2_mean_pct", "unique_read_enrichment_pct", "unique_base_enrichment_pct",
        "aligned_reads_in_target_pct", "provider_callability_pct", "provider_autosome_callability_pct",
        "somatic_callable_bases", "estimated_chr_x_target_depth_x", "estimated_chr_y_target_depth_x",
        "fragment_length_median_bp", "dragen_germline_snvs", "dragen_germline_indels",
        "dragen_germline_titv", "bam_present", "bam_bytes", "hard_filtered_vcf_present",
        "source_summary_csv", "source_callable_bed",
    ]
    pair_fields = [
        "patient_id", "tumour_sample_id", "normal_sample_id", "tumour_mean_target_depth_x",
        "normal_mean_target_depth_x", "tumour_target_20x_pct", "normal_target_20x_pct",
        "tumour_duplicate_pct", "normal_duplicate_pct", "tumour_provider_callability_pct",
        "normal_provider_callability_pct", "pair_callable_bases", "pair_callable_mb",
        "pair_callable_pct_of_target", "descriptive_qc_review", "review_reasons", "purple_purity",
        "purple_status", "purple_qc_status", "onco_run_mode", "onco_targeted",
        "onco_reported_tmb_per_mb_not_used", "onco_reported_tml",
    ]
    call_fields = [
        "patient_id", "caller", "vcf_records", "vcf_pass_records", "pair_callable_records",
        "pass_pair_callable_records", "pass_pair_callable_snvs", "pass_pair_callable_mnvs",
        "pass_pair_callable_insertions", "pass_pair_callable_deletions",
        "pass_pair_callable_nonsynonymous", "pass_pair_callable_hotspots",
        "strict_read_evidence_records", "strict_read_evidence_snvs", "strict_read_evidence_mnvs",
        "strict_read_evidence_insertions", "strict_read_evidence_deletions",
        "strict_read_evidence_nonsynonymous", "strict_read_evidence_hotspots",
        "median_tumour_depth_at_pass_calls", "median_normal_depth_at_pass_calls",
        "median_tumour_vaf_at_pass_calls", "median_normal_vaf_at_pass_calls",
        "pair_callable_bases", "pair_callable_mb", "pass_variants_per_pair_callable_mb",
        "pass_nonsynonymous_per_pair_callable_mb", "strict_read_evidence_per_pair_callable_mb",
        "strict_nonsynonymous_per_pair_callable_mb", "onco_targeted_mode",
        "clinical_tmb_interpretation", "source_vcf",
    ]
    variant_fields = [
        "patient_id", "caller", "chrom", "position", "ref", "alt", "variant_type", "filter",
        "qual", "tier", "gene", "transcript", "consequence", "coding_effect", "nonsynonymous",
        "tumour_depth", "tumour_alt_depth", "tumour_vaf", "normal_depth", "normal_alt_depth",
        "normal_vaf", "strict_read_evidence", "purple_copy_number", "purple_variant_copy_number", "biallelic_flag",
        "source_vcf",
    ]
    germline_fields = [
        "patient_id", "normal_sample_id", "dragen_normal_germline_snvs",
        "dragen_normal_germline_indels", "dragen_normal_germline_titv",
        "purple_somatic_driver_catalog_rows", "purple_somatic_driver_likelihood_ge_0_5_rows",
        "purple_somatic_driver_mutation_rows", "purple_somatic_driver_deletion_rows",
        "purple_germline_driver_catalog_rows", "purple_germline_driver_likelihood_ge_0_5_rows",
        "purple_germline_driver_mutation_rows", "purple_germline_driver_deletion_rows",
        "germline_candidate_vcf_present", "interpretation",
    ]
    input_fields = [
        "patient_id", "tumour_enrichment_metrics", "normal_enrichment_metrics", "tumour_bam",
        "normal_bam", "sage_somatic_vcf", "purple_somatic_vcf", "purple_purity",
        "somatic_driver_catalog", "germline_driver_catalog",
    ]

    purple_by_patient = {
        str(row["patient_id"]): row for row in call_summaries if row["caller"] == "PURPLE"
    }
    driver_by_patient = {str(row["patient_id"]): row for row in germline_rows}
    burden_logs = [
        math.log10(float(row["strict_read_evidence_per_pair_callable_mb"]) + 0.01)
        for row in purple_by_patient.values()
    ]
    burden_log_median = statistics.median(burden_logs) if burden_logs else None
    burden_log_mad = (
        statistics.median(abs(value - burden_log_median) for value in burden_logs)
        if burden_log_median is not None
        else None
    )
    patient_summary_rows: list[dict[str, object]] = []
    for pair_row in pair_rows:
        patient_id = str(pair_row["patient_id"])
        calls = purple_by_patient.get(patient_id, {})
        drivers = driver_by_patient.get(patient_id, {})
        strict_burden = clean_number(calls.get("strict_read_evidence_per_pair_callable_mb"))
        robust_z = None
        if strict_burden is not None and burden_log_median is not None and burden_log_mad:
            robust_z = (math.log10(strict_burden + 0.01) - burden_log_median) / (1.4826 * burden_log_mad)
        limitations: list[str] = []
        if not calls:
            limitations.append("final_somatic_vcf_missing_after_failed_result_publication")
        if pair_row.get("descriptive_qc_review") == "REVIEW":
            limitations.append(str(pair_row.get("review_reasons", "qc_review")))
        if str(pair_row.get("onco_targeted", "NA")).lower() == "false":
            limitations.append("oncoanalyser_targeted_false")
        purple_status = str(pair_row.get("purple_status", "NA"))
        if purple_status not in {"SOMATIC", "NA"}:
            limitations.append(f"purple_status_{purple_status.lower()}")
        purple_qc_status = str(pair_row.get("purple_qc_status", "NA"))
        if purple_qc_status not in {"PASS", "NA"}:
            limitations.append("purple_qc_warning")
        patient_summary_rows.append(
            {
                "patient_id": patient_id,
                "tumour_sample_id": pair_row.get("tumour_sample_id"),
                "normal_sample_id": pair_row.get("normal_sample_id"),
                "tumour_mean_target_depth_x": pair_row.get("tumour_mean_target_depth_x"),
                "normal_mean_target_depth_x": pair_row.get("normal_mean_target_depth_x"),
                "tumour_target_20x_pct": pair_row.get("tumour_target_20x_pct"),
                "normal_target_20x_pct": pair_row.get("normal_target_20x_pct"),
                "pair_callable_mb": pair_row.get("pair_callable_mb"),
                "pair_callable_pct_of_target": pair_row.get("pair_callable_pct_of_target"),
                "purple_purity": pair_row.get("purple_purity"),
                "purple_status": pair_row.get("purple_status"),
                "purple_qc_status": pair_row.get("purple_qc_status"),
                "final_pass_pair_callable_somatic_records": calls.get("pass_pair_callable_records"),
                "final_strict_read_evidence_somatic_records": calls.get("strict_read_evidence_records"),
                "final_strict_read_evidence_nonsynonymous_records": calls.get(
                    "strict_read_evidence_nonsynonymous"
                ),
                "strict_read_evidence_variants_per_callable_target_mb": strict_burden,
                "strict_nonsynonymous_variants_per_callable_target_mb": calls.get(
                    "strict_nonsynonymous_per_pair_callable_mb"
                ),
                "strict_burden_robust_z_on_log_scale": robust_z,
                "strict_burden_cohort_outlier": bool(
                    robust_z is not None and robust_z >= args.outlier_robust_z
                ),
                "purple_somatic_driver_catalog_rows_unvalidated": drivers.get(
                    "purple_somatic_driver_catalog_rows"
                ),
                "purple_somatic_driver_likelihood_ge_0_5_rows_unvalidated": drivers.get(
                    "purple_somatic_driver_likelihood_ge_0_5_rows"
                ),
                "somatic_call_status": "COMPLETE_WITH_LIMITATIONS" if calls else "MISSING_FINAL_SOMATIC_VCF",
                "limitations": ";".join(limitations) if limitations else "none",
                "clinical_tmb": "NOT_CALCULATED",
            }
        )

    patient_summary_fields = [
        "patient_id", "tumour_sample_id", "normal_sample_id", "tumour_mean_target_depth_x",
        "normal_mean_target_depth_x", "tumour_target_20x_pct", "normal_target_20x_pct",
        "pair_callable_mb", "pair_callable_pct_of_target", "purple_purity", "purple_status",
        "purple_qc_status", "final_pass_pair_callable_somatic_records",
        "final_strict_read_evidence_somatic_records",
        "final_strict_read_evidence_nonsynonymous_records",
        "strict_read_evidence_variants_per_callable_target_mb",
        "strict_nonsynonymous_variants_per_callable_target_mb",
        "strict_burden_robust_z_on_log_scale", "strict_burden_cohort_outlier",
        "purple_somatic_driver_catalog_rows_unvalidated",
        "purple_somatic_driver_likelihood_ge_0_5_rows_unvalidated",
        "somatic_call_status", "limitations", "clinical_tmb",
    ]

    write_tsv(output_dir / "wes_sample_qc.tsv", sample_qc, sample_fields)
    write_tsv(output_dir / "wes_matched_pair_qc.tsv", pair_rows, pair_fields)
    write_tsv(output_dir / "somatic_call_summary.tsv", call_summaries, call_fields)
    write_tsv(output_dir / "somatic_filter_counts.tsv", filter_rows, ["patient_id", "caller", "filter", "record_count"])
    write_tsv(output_dir / "somatic_pass_pair_callable_variants.tsv", somatic_calls, variant_fields)
    write_tsv(output_dir / "germline_and_driver_call_summary.tsv", germline_rows, germline_fields)
    write_tsv(output_dir / "input_completeness.tsv", input_rows, input_fields)
    write_tsv(
        output_dir / "patient_level_wes_genomics_summary.tsv",
        patient_summary_rows,
        patient_summary_fields,
    )

    purple_summaries = [row for row in call_summaries if row["caller"] == "PURPLE"]
    depth_values = [clean_number(row["mean_target_depth_x"]) for row in sample_qc]
    depth_values = [value for value in depth_values if value is not None]
    pair_callable_values = [clean_number(row["pair_callable_mb"]) for row in pair_rows]
    pair_callable_values = [value for value in pair_callable_values if value is not None]
    pass_counts = [int(row["pass_pair_callable_records"]) for row in purple_summaries]
    strict_counts = [int(row["strict_read_evidence_records"]) for row in purple_summaries]
    review_count = sum(row["descriptive_qc_review"] == "REVIEW" for row in pair_rows)
    targeted_false = sum(str(row["onco_targeted"]).lower() == "false" for row in pair_rows)

    summary_lines = [
        "# WES depth and genomic call audit",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        f"Samples with enrichment metrics: {len(sample_qc)}",
        f"Matched tumour-normal pairs: {len(pair_rows)}",
        f"Patients with final PURPLE somatic summaries: {len(purple_summaries)}",
        f"Median sample mean target depth: {statistics.median(depth_values):.2f}x" if depth_values else "Median sample mean target depth: NA",
        f"Median matched-pair callable territory: {statistics.median(pair_callable_values):.2f} Mb" if pair_callable_values else "Median matched-pair callable territory: NA",
        f"Median final PASS pair-callable somatic records: {statistics.median(pass_counts):.1f}" if pass_counts else "Median final PASS pair-callable somatic records: NA",
        f"Median final strict-read-evidence somatic records: {statistics.median(strict_counts):.1f}" if strict_counts else "Median final strict-read-evidence somatic records: NA",
        f"Pairs flagged for descriptive QC review: {review_count}/{len(pair_rows)}",
        f"Oncoanalyser records marked targeted=false: {targeted_false}/{len(pair_rows)}",
        f"Descriptive log-burden outlier robust-z threshold: {args.outlier_robust_z}",
        "",
        "The reported Oncoanalyser TMB fields were retained only for provenance and were not used as valid WES TMB values when targeted=false. The recalculated burden is the number of PASS variants within territory callable in both tumour and matched normal divided by the full callable Twist target size. It is not a clinical TMB because the denominator is not restricted to a validated coding territory.",
        "The strict read-evidence sensitivity count additionally requires tumour depth >=20, normal depth >=10, tumour alternate depth >=3, tumour VAF >=0.05, normal alternate depth <=1, and normal VAF <=0.02.",
        "",
        "Review thresholds are descriptive flags, not clinical acceptance criteria.",
    ]
    (output_dir / "wes_qc_and_call_audit_summary.md").write_text("\n".join(summary_lines) + "\n")

    manifest_rows = []
    for path in sorted(output_dir.iterdir()):
        if path.is_file():
            manifest_rows.append({"file": path.name, "bytes": path.stat().st_size})
    write_tsv(output_dir / "output_manifest.tsv", manifest_rows, ["file", "bytes"])

    print("\n".join(summary_lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
