#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re


WES_RE = re.compile(
    r"^(?P<label>\d+)(?P<status>[NT])_\d+WESFastq_S\d+_(?P<lane>L\d{3})_R(?P<read>[12])_001\.fastq\.gz$"
)
RNA_RE = re.compile(
    r"^(?P<label>\d+)(?P<status>[NT])_\d+RNA_S\d+_(?P<lane>L\d{3})_R(?P<read>[12])_001\.fastq\.gz$"
)
WES_PLOIDY_RE = re.compile(r"^(?P<label>\d+)(?P<status>[NT])WESFastq\.ploidy_estimation_metrics\.csv$")
TWIST_BED_NAME = "Twist_ILMN_Exome_2.5_Plus_Panel.hg38.bed"
REPORT_SAMPLE_RE = re.compile(r"^(?P<label>\d+)(?P<status>[NT])RNA(?:_\d+)?$")


@dataclass(frozen=True)
class SampleKey:
    patient: str
    status_code: str
    lane: str


@dataclass(frozen=True)
class IntervalConsensus:
    source_path: str
    consensus_hash: str
    consensus_count: int
    total_count: int
    outlier_paths: tuple[str, ...]
    link_path: str | None = None


def backup_if_exists(path: Path, backup_root: Path) -> None:
    if not path.exists():
        return
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_dir = backup_root / stamp
    backup_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, backup_dir / path.name)


def write_text(path: Path, text: str, backup_root: Path) -> None:
    backup_if_exists(path, backup_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]], comments: list[str], backup_root: Path) -> None:
    backup_if_exists(path, backup_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        for line in comments:
            handle.write(f"# {line}\n")
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def discover_fastqs(input_root: Path, regex: re.Pattern[str], skip_dirs: set[str]) -> dict[SampleKey, dict[str, str]]:
    found: dict[SampleKey, dict[str, str]] = defaultdict(dict)
    for fastq in sorted(input_root.rglob("*.fastq.gz")):
        if any(part in skip_dirs for part in fastq.parts):
            continue
        match = regex.match(fastq.name)
        if not match:
            continue
        key = SampleKey(
            patient=match.group("label"),
            status_code=match.group("status"),
            lane=match.group("lane"),
        )
        found[key][match.group("read")] = str(fastq)
    return found


def count_matching_files(input_root: Path, suffixes: tuple[str, ...], include_parts: tuple[str, ...], skip_dirs: set[str]) -> list[Path]:
    hits: list[Path] = []
    for path in sorted(input_root.rglob("*")):
        if not path.is_file():
            continue
        if any(part in skip_dirs for part in path.parts):
            continue
        if suffixes and not path.name.endswith(suffixes):
            continue
        if include_parts and not any(part in path.parts for part in include_parts):
            continue
        hits.append(path)
    return hits


def infer_sex_from_ploidy_metrics(path: Path) -> str | None:
    values: dict[str, float] = {}
    with path.open() as handle:
        for line in handle:
            parts = [part.strip() for part in line.rstrip().split(",")]
            if len(parts) >= 4 and parts[0] == "PLOIDY ESTIMATION":
                try:
                    values[parts[2]] = float(parts[3])
                except ValueError:
                    continue
    autosomal = values.get("Autosomal median coverage")
    x_cov = values.get("X median coverage")
    y_cov = values.get("Y median coverage")
    if not autosomal or x_cov is None or y_cov is None:
        return None
    x_ratio = x_cov / autosomal
    y_ratio = y_cov / autosomal
    if y_ratio < 0.05 and x_ratio > 0.7:
        return "XX"
    if y_ratio > 0.10 and 0.25 < x_ratio < 0.7:
        return "XY"
    return None


def infer_patient_sex_map(input_root: Path, skip_dirs: set[str]) -> tuple[dict[str, str], list[str]]:
    normal_calls: dict[str, str] = {}
    tumor_calls: dict[str, str] = {}
    for path in sorted(input_root.rglob("*.ploidy_estimation_metrics.csv")):
        if any(part in skip_dirs for part in path.parts):
            continue
        match = WES_PLOIDY_RE.match(path.name)
        if not match:
            continue
        sex = infer_sex_from_ploidy_metrics(path)
        if sex is None:
            continue
        patient = match.group("label")
        status = match.group("status")
        if status == "N":
            normal_calls[patient] = sex
        else:
            tumor_calls[patient] = sex

    notes: list[str] = []
    discordant = sorted(
        patient for patient, sex in tumor_calls.items() if patient in normal_calls and normal_calls[patient] != sex
    )
    if discordant:
        notes.append(
            "Tumor ploidy was not used for patient sex inference because the following tumor calls differed from their normal: "
            + ", ".join(discordant)
        )
    return normal_calls, notes


def discover_interval_consensus(input_root: Path, skip_dirs: set[str]) -> IntervalConsensus | None:
    candidates = [
        path
        for path in sorted(input_root.rglob(TWIST_BED_NAME))
        if path.is_file() and not any(part in skip_dirs for part in path.parts)
    ]
    if not candidates:
        return None

    by_hash: dict[str, list[Path]] = defaultdict(list)
    for path in candidates:
        digest = hashlib.md5(path.read_bytes()).hexdigest()
        by_hash[digest].append(path)

    consensus_hash, consensus_paths = max(by_hash.items(), key=lambda item: (len(item[1]), str(item[1][0])))
    outlier_paths = tuple(
        str(path) for digest, paths in by_hash.items() if digest != consensus_hash for path in sorted(paths)
    )
    return IntervalConsensus(
        source_path=str(sorted(consensus_paths)[0]),
        consensus_hash=consensus_hash,
        consensus_count=len(consensus_paths),
        total_count=len(candidates),
        outlier_paths=outlier_paths,
    )


def ensure_symlink(path: Path, source: Path, backup_root: Path) -> None:
    if path.is_symlink():
        try:
            if path.resolve() == source.resolve():
                return
        except FileNotFoundError:
            pass
    backup_if_exists(path, backup_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        path.unlink()
    path.symlink_to(source)


def nfcore_strandedness_from_dragen_orientation(code: str) -> str | None:
    mapping = {
        "IU": "unstranded",
        "ISR": "reverse",
        "ISF": "forward",
        "U": "unstranded",
        "SR": "reverse",
        "SF": "forward",
    }
    return mapping.get(code)


def discover_rna_strandedness_map(input_root: Path, skip_dirs: set[str]) -> tuple[dict[str, str], list[str]]:
    html_candidates = [
        path
        for path in sorted(input_root.rglob("dragen-reports.html"))
        if path.is_file() and "DRAGEN_RNA_" in str(path) and not any(part in skip_dirs for part in path.parts)
    ]
    notes: list[str] = []
    if not html_candidates:
        return {}, notes

    sample_map: dict[str, str] = {}
    orientation_counts: dict[str, int] = defaultdict(int)
    for html_path in html_candidates:
        text = html_path.read_text(errors="ignore")
        for match in re.finditer(r"rowData:\s*(\[[^\n]*\])", text):
            block = match.group(1)
            if '"library_orientation"' not in block:
                continue
            try:
                rows = json.loads(block)
            except json.JSONDecodeError:
                continue
            for row in rows:
                report_sample = row.get("sample", "")
                orientation = row.get("library_orientation", "")
                strandedness = nfcore_strandedness_from_dragen_orientation(orientation)
                if not strandedness:
                    continue
                sample_match = REPORT_SAMPLE_RE.match(report_sample)
                if not sample_match:
                    continue
                sample_name = f"{sample_match.group('label')}{sample_match.group('status')}_RNA"
                sample_map[sample_name] = strandedness
                orientation_counts[orientation] += 1
            break

    if orientation_counts:
        summary = ", ".join(f"{code}={orientation_counts[code]}" for code in sorted(orientation_counts))
        notes.append(
            "RNA strandedness was populated from the DRAGEN RNA report library_orientation field "
            f"({summary})."
        )
    return sample_map, notes


def build_wes_rows(pairs: dict[SampleKey, dict[str, str]], patient_sex_map: dict[str, str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for key in sorted(pairs, key=lambda item: (int(item.patient), item.status_code, item.lane)):
        reads = pairs[key]
        if "1" not in reads or "2" not in reads:
            continue
        rows.append(
            {
                "patient": key.patient,
                "sex": patient_sex_map.get(key.patient, ""),
                "status": "0" if key.status_code == "N" else "1",
                "sample": f"{key.patient}{key.status_code}_WES",
                "lane": key.lane,
                "fastq_1": reads["1"],
                "fastq_2": reads["2"],
                "bam": "",
                "bai": "",
                "cram": "",
                "crai": "",
            }
        )
    return rows


def build_rna_rows(pairs: dict[SampleKey, dict[str, str]], strandedness_map: dict[str, str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for key in sorted(pairs, key=lambda item: (int(item.patient), item.status_code, item.lane)):
        reads = pairs[key]
        if "1" not in reads or "2" not in reads:
            continue
        sample_name = f"{key.patient}{key.status_code}_RNA"
        rows.append(
            {
                "sample": sample_name,
                "fastq_1": reads["1"],
                "fastq_2": reads["2"],
                "strandedness": strandedness_map.get(sample_name, "REVIEW_ME"),
            }
        )
    return rows


def build_rnafusion_rows(pairs: dict[SampleKey, dict[str, str]], strandedness_map: dict[str, str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for key in sorted(pairs, key=lambda item: (int(item.patient), item.status_code, item.lane)):
        reads = pairs[key]
        if "1" not in reads or "2" not in reads:
            continue
        sample_name = f"{key.patient}{key.status_code}_RNA"
        rows.append(
            {
                "sample": sample_name,
                "fastq_1": reads["1"],
                "fastq_2": reads["2"],
                "bam": "",
                "bai": "",
                "cram": "",
                "crai": "",
                "junctions": "",
                "splice_junctions": "",
                "strandedness": strandedness_map.get(sample_name, "REVIEW_ME"),
                "seq_platform": "ILLUMINA",
                "seq_center": "LOCAL",
            }
        )
    return rows


def launch_drafts_markdown(
    workspace_base: Path, interval_consensus: IntervalConsensus | None, unresolved_rna_rows: int
) -> str:
    scripts_dir = workspace_base / "scripts"
    samplesheets_dir = workspace_base / "samplesheets"
    bulk_results_root = Path(os.environ.get("PDAC2026_RESULTS_ROOT", workspace_base / "results"))
    interval_path = interval_consensus.link_path if interval_consensus and interval_consensus.link_path else (
        interval_consensus.source_path if interval_consensus else "/REVIEW_ME/EXOME_CAPTURE_INTERVALS.bed"
    )
    rna_note = (
        "Edit the RNA samplesheet first and replace every `REVIEW_ME` with `unstranded`, `forward`, or `reverse`."
        if unresolved_rna_rows
        else "RNA strandedness was auto-filled from the DRAGEN RNA report; spot-check the samplesheet before launch."
    )
    rnafusion_note = (
        "Edit the RNA fusion samplesheet first and replace every `REVIEW_ME` with `unstranded`, `forward`, or `reverse`."
        if unresolved_rna_rows
        else "RNA fusion strandedness was auto-filled from the DRAGEN RNA report; spot-check the samplesheet before launch."
    )
    return "\n".join(
        [
            "# PDAC Launch Drafts",
            "",
            "- Review the generated samplesheets before launching.",
            "- Keep patient data in place. These commands read directly from the source FASTQs.",
            f"- Suggested bulk-analysis output root: `{bulk_results_root}`",
            "- Keep every `--outdir` on a restricted Linux-native filesystem outside the public repository.",
            "",
            "## WES Tumor-Normal",
            "",
            "```bash",
            f"{scripts_dir}/run_sarek.sh --mode tumor-normal \\",
            f"  --samplesheet {samplesheets_dir}/sarek_samplesheet.PDAC_WES_fastq_autodraft.csv \\",
            f"  --intervals '{interval_path}' \\",
            f"  --outdir {bulk_results_root}/sarek_tumor_normal \\",
            "  -resume",
            "```",
            "",
            "## WES Germline",
            "",
            "```bash",
            f"{scripts_dir}/run_sarek.sh --mode germline \\",
            f"  --samplesheet {samplesheets_dir}/sarek_samplesheet.PDAC_WES_fastq_autodraft.csv \\",
            f"  --intervals '{interval_path}' \\",
            f"  --outdir {bulk_results_root}/sarek_germline \\",
            "  -resume",
            "```",
            "",
            "## RNA-Seq Expression",
            "",
            rna_note,
            "",
            "```bash",
            f"{scripts_dir}/run_rnaseq.sh \\",
            f"  --samplesheet {samplesheets_dir}/rnaseq_samplesheet.PDAC_RNA_fastq_autodraft.csv \\",
            f"  --outdir {bulk_results_root}/rnaseq_expression \\",
            "  -resume",
            "```",
            "",
            "## RNA Fusion",
            "",
            rnafusion_note,
            "",
            "```bash",
            f"{scripts_dir}/run_rnafusion.sh \\",
            f"  --samplesheet {samplesheets_dir}/rnafusion_samplesheet.PDAC_RNA_fastq_autodraft.csv \\",
            f"  --outdir {workspace_base}/results/rnafusion_pdac \\",
            "  -resume",
            "```",
        ]
    )


def inventory_markdown(
    workspace_base: Path,
    input_root: Path,
    wes_rows: list[dict[str, str]],
    rna_rows: list[dict[str, str]],
    rna_bams: list[Path],
    wes_bams: list[Path],
    wes_vcfs: list[Path],
    patient_sex_map: dict[str, str],
    sex_notes: list[str],
    rna_strandedness_map: dict[str, str],
    rna_notes: list[str],
    interval_consensus: IntervalConsensus | None,
) -> str:
    tumor_ids = sorted({row["patient"] for row in wes_rows if row["status"] == "1"}, key=int)
    normal_ids = sorted({row["patient"] for row in wes_rows if row["status"] == "0"}, key=int)
    generated_files = [
        workspace_base / "samplesheets" / "sarek_samplesheet.PDAC_WES_fastq_autodraft.csv",
        workspace_base / "samplesheets" / "rnaseq_samplesheet.PDAC_RNA_fastq_autodraft.csv",
        workspace_base / "samplesheets" / "rnafusion_samplesheet.PDAC_RNA_fastq_autodraft.csv",
        workspace_base / "docs" / "PDAC_LAUNCH_DRAFTS.md",
    ]
    lines = [
        "# PDAC Input Inventory",
        "",
        f"- Input root: `{input_root}`",
        f"- WES FASTQ rows drafted: {len(wes_rows)}",
        f"- RNA FASTQ rows drafted: {len(rna_rows)}",
        f"- RNA strandedness filled from DRAGEN report: {sum(1 for row in rna_rows if row['strandedness'] != 'REVIEW_ME')}/{len(rna_rows)}",
        f"- WES BAMs discovered: {len(wes_bams)}",
        f"- RNA BAMs discovered: {len(rna_bams)}",
        f"- WES hard-filtered VCFs discovered: {len(wes_vcfs)}",
        f"- WES patient sex inferred from normal ploidy metrics: {len(patient_sex_map)}/{len(normal_ids)}",
        "",
        "## Cohort IDs",
        "",
        f"- Normal labels: {', '.join(normal_ids)}",
        f"- Tumor labels: {', '.join(tumor_ids)}",
        "",
        "## Generated draft files",
        "",
    ]
    lines.extend(f"- `{path}`" for path in generated_files)
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- These drafts do not move or modify source data.",
            "- WES `sex` is filled only when a matching normal ploidy metric supports an `XX` or `XY` inference.",
            "- RNA `strandedness` prefers the DRAGEN RNA report `library_orientation` field when available; any remaining unknown rows stay as `REVIEW_ME`.",
            "- Existing DRAGEN BAM/VCF outputs were inventoried but not used as primary workflow inputs in these drafts.",
        ]
    )
    if interval_consensus is not None:
        interval_label = interval_consensus.link_path or interval_consensus.source_path
        lines.append(
            f"- Consensus exome interval candidate ({interval_consensus.consensus_count}/{interval_consensus.total_count} copies): `{interval_label}`"
        )
        if interval_consensus.outlier_paths:
            lines.append(f"- Interval outlier kept for review: `{interval_consensus.outlier_paths[0]}`")
    lines.extend(f"- {note}" for note in sex_notes)
    lines.extend(f"- {note}" for note in rna_notes)
    lines.extend(
        [
            "",
            "## Inferred WES Sex By Patient",
            "",
        ]
    )
    for patient in sorted(patient_sex_map, key=int):
        lines.append(f"- {patient}: {patient_sex_map[patient]}")
    lines.extend(
        [
            "",
            "## Example existing processed files",
            "",
        ]
    )
    for heading, files in (
        ("WES BAM examples", wes_bams[:5]),
        ("RNA BAM examples", rna_bams[:5]),
        ("WES VCF examples", wes_vcfs[:5]),
    ):
        lines.append(f"### {heading}")
        lines.append("")
        if files:
            lines.extend(f"- `{path}`" for path in files)
        else:
            lines.append("- None discovered")
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", required=True)
    parser.add_argument("--samplesheets-dir", required=True)
    parser.add_argument("--docs-dir", required=True)
    parser.add_argument("--backup-dir", required=True)
    parser.add_argument("--link-intervals-dir")
    args = parser.parse_args()

    input_root = Path(args.input_root)
    samplesheets_dir = Path(args.samplesheets_dir)
    docs_dir = Path(args.docs_dir)
    backup_dir = Path(args.backup_dir)
    link_intervals_dir = Path(args.link_intervals_dir) if args.link_intervals_dir else None
    skip_dirs = {"SEQ_analysis", "work", "results", "smoke_tests"}
    workspace_base = samplesheets_dir.parent

    wes_pairs = discover_fastqs(input_root, WES_RE, skip_dirs)
    rna_pairs = discover_fastqs(input_root, RNA_RE, skip_dirs)
    patient_sex_map, sex_notes = infer_patient_sex_map(input_root, skip_dirs)
    rna_strandedness_map, rna_notes = discover_rna_strandedness_map(input_root, skip_dirs)
    interval_consensus = discover_interval_consensus(input_root, skip_dirs)
    if interval_consensus is not None and link_intervals_dir is not None:
        link_path = link_intervals_dir / "PDAC_Twist_ILMN_Exome_2.5_Plus_Panel.hg38.majority.bed"
        ensure_symlink(link_path, Path(interval_consensus.source_path), backup_dir)
        interval_consensus = IntervalConsensus(
            source_path=interval_consensus.source_path,
            consensus_hash=interval_consensus.consensus_hash,
            consensus_count=interval_consensus.consensus_count,
            total_count=interval_consensus.total_count,
            outlier_paths=interval_consensus.outlier_paths,
            link_path=str(link_path),
        )

    wes_rows = build_wes_rows(wes_pairs, patient_sex_map)
    rna_rows = build_rna_rows(rna_pairs, rna_strandedness_map)
    rnafusion_rows = build_rnafusion_rows(rna_pairs, rna_strandedness_map)

    wes_bams = count_matching_files(input_root, (".bam",), ("1900NGS_WES_Enrichment_stats-10282273",), skip_dirs)
    rna_bams = count_matching_files(input_root, (".bam",), ("DRAGEN_RNA_03_17_2025_22_02_16-59174132",), skip_dirs)
    wes_vcfs = count_matching_files(input_root, (".hard-filtered.vcf.gz",), ("1900NGS_WES_Enrichment_stats-10282273",), skip_dirs)

    write_csv(
        samplesheets_dir / "sarek_samplesheet.PDAC_WES_fastq_autodraft.csv",
        ["patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2", "bam", "bai", "cram", "crai"],
        wes_rows,
        [
            "Auto-drafted from the PDAC input tree. Review before use.",
            "sex is inferred from normal-sample ploidy metrics when available; confirm before use.",
            "status uses 0 for normal and 1 for tumor.",
        ],
        backup_dir,
    )
    write_csv(
        samplesheets_dir / "rnaseq_samplesheet.PDAC_RNA_fastq_autodraft.csv",
        ["sample", "fastq_1", "fastq_2", "strandedness"],
        rna_rows,
        [
            "Auto-drafted from the PDAC input tree. Review before use.",
            "Strandedness is auto-filled from the DRAGEN RNA report when available; confirm any remaining REVIEW_ME rows before launching.",
        ],
        backup_dir,
    )
    write_csv(
        samplesheets_dir / "rnafusion_samplesheet.PDAC_RNA_fastq_autodraft.csv",
        ["sample", "fastq_1", "fastq_2", "bam", "bai", "cram", "crai", "junctions", "splice_junctions", "strandedness", "seq_platform", "seq_center"],
        rnafusion_rows,
        [
            "Auto-drafted from the PDAC input tree. Review before use.",
            "Strandedness is auto-filled from the DRAGEN RNA report when available; confirm any remaining REVIEW_ME rows before launching.",
        ],
        backup_dir,
    )
    write_text(
        docs_dir / "PDAC_INPUT_INVENTORY.md",
        inventory_markdown(
            workspace_base,
            input_root,
            wes_rows,
            rna_rows,
            rna_bams,
            wes_bams,
            wes_vcfs,
            patient_sex_map,
            sex_notes,
            rna_strandedness_map,
            rna_notes,
            interval_consensus,
        ),
        backup_dir,
    )
    unresolved_rna_rows = sum(1 for row in rna_rows if row["strandedness"] == "REVIEW_ME")
    write_text(
        docs_dir / "PDAC_LAUNCH_DRAFTS.md",
        launch_drafts_markdown(workspace_base, interval_consensus, unresolved_rna_rows),
        backup_dir,
    )

    print(f"Wrote {len(wes_rows)} WES rows, {len(rna_rows)} RNA rows, and inventory docs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
