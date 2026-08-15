#!/usr/bin/env python3
import argparse
import csv
import gzip
import shutil
from pathlib import Path

import pandas as pd
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerMatrixGenerator.scripts.CNVMatrixGenerator import generateCNVMatrix
from SigProfilerMatrixGenerator.scripts.SVMatrixGenerator import generateSVMatrix


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate CNV48 and SV32 COSMIC assignments from Sarek ASCAT and Manta outputs."
    )
    parser.add_argument("--sarek-root", required=True, help="Root of the Sarek SV/CNA result tree")
    parser.add_argument("--output", required=True, help="Output directory root")
    parser.add_argument("--project", default="PDAC2026")
    parser.add_argument("--genome-build", default="GRCh38")
    parser.add_argument("--cosmic-version", type=float, default=3.5)
    parser.add_argument("--cpu", type=int, default=8)
    parser.add_argument("--skip-completed", action="store_true")
    parser.add_argument("--make-plots", action="store_true")
    return parser.parse_args()


def normalize(name: str) -> str:
    return "".join(ch.lower() for ch in name if ch.isalnum())


def choose_column(columns, candidates):
    lookup = {normalize(col): col for col in columns}
    for candidate in candidates:
        if candidate in lookup:
            return lookup[candidate]
    return None


def strip_suffixes(name: str, suffixes):
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def is_completed(output_root: Path, context_name: str) -> bool:
    metadata = output_root / context_name / "JOB_METADATA_SPA.txt"
    if not metadata.is_file():
        return False
    return "completed successfully" in metadata.read_text(errors="ignore").lower()


def collect_ascat_segments(sarek_root: Path, merged_segments_path: Path):
    segment_files = sorted(sarek_root.rglob("*.segments.txt"))
    if not segment_files:
        raise SystemExit(f"No ASCAT segment files found under {sarek_root}")

    rows = []
    for segment_file in segment_files:
        sample = strip_suffixes(segment_file.name, [".segments.txt"])
        df = pd.read_csv(segment_file, sep="\t")

        sample_col = choose_column(df.columns, ["sample", "samples"])
        chr_col = choose_column(df.columns, ["chr", "chromosome", "chrom"])
        start_col = choose_column(df.columns, ["startpos", "startposition", "start"])
        end_col = choose_column(df.columns, ["endpos", "endposition", "end"])
        nmajor_col = choose_column(df.columns, ["nmajor", "majorcn", "majorcopynumber"])
        nminor_col = choose_column(df.columns, ["nminor", "minorcn", "minorcopynumber"])

        required = [chr_col, start_col, end_col, nmajor_col, nminor_col]
        if any(col is None for col in required):
            raise SystemExit(
                f"Could not map ASCAT segment columns in {segment_file}. "
                f"Found columns: {list(df.columns)}"
            )

        normalized = pd.DataFrame(
            {
                "sample": df[sample_col] if sample_col else sample,
                "chr": df[chr_col],
                "startpos": df[start_col],
                "endpos": df[end_col],
                "nMajor": df[nmajor_col],
                "nMinor": df[nminor_col],
            }
        )
        rows.append(normalized)

    merged = pd.concat(rows, ignore_index=True)
    merged_segments_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(merged_segments_path, sep="\t", index=False)
    return merged_segments_path


def collect_manta_vcfs(sarek_root: Path, output_dir: Path):
    patterns = [
        "*somaticSV.vcf.gz",
        "*.manta.somatic_sv.vcf.gz",
    ]
    manta_vcfs = []
    for pattern in patterns:
        manta_vcfs.extend(sarek_root.rglob(pattern))
    manta_vcfs = sorted(set(manta_vcfs))
    if not manta_vcfs:
        raise SystemExit(f"No Manta somatic SV VCFs found under {sarek_root}")

    output_dir.mkdir(parents=True, exist_ok=True)
    for vcf_gz in manta_vcfs:
        sample = strip_suffixes(
            vcf_gz.name,
            [".manta.somatic_sv.vcf.gz", ".somaticSV.vcf.gz", ".vcf.gz"],
        )
        dest = output_dir / f"{sample}.vcf"
        if dest.is_file() and dest.stat().st_size > 0:
            continue
        with gzip.open(vcf_gz, "rb") as src, open(dest, "wb") as dst:
            shutil.copyfileobj(src, dst)
    return output_dir


def run_cnv_assignment(sarek_root: Path, output_root: Path, project: str, genome_build: str, cosmic_version: float, cpu: int, make_plots: bool):
    matrix_root = output_root / "CNV48_matrix"
    merged_segments = matrix_root / f"{project}.ascat_segments.tsv"
    collect_ascat_segments(sarek_root, merged_segments)
    generateCNVMatrix("ASCAT", str(merged_segments), project, str(matrix_root))

    matrix_path = matrix_root / f"{project}.CNV48.matrix.tsv"
    Analyze.cosmic_fit(
        samples=str(matrix_path),
        output=str(output_root / "CNV48"),
        input_type="matrix",
        context_type="48",
        cosmic_version=cosmic_version,
        exome=False,
        genome_build=genome_build,
        signature_database=None,
        exclude_signature_subgroups=None,
        collapse_to_SBS96=False,
        export_probabilities=False,
        export_probabilities_per_mutation=False,
        make_plots=make_plots,
        sample_reconstruction_plots="none",
        verbose=True,
        cpu=cpu,
    )


def run_sv_assignment(sarek_root: Path, output_root: Path, project: str, genome_build: str, cosmic_version: float, cpu: int, make_plots: bool):
    matrix_root = output_root / "SV32_matrix"
    manta_input = collect_manta_vcfs(sarek_root, matrix_root / "manta_vcfs_plain")
    generateSVMatrix(str(manta_input), project, str(matrix_root))

    matrix_path = matrix_root / f"{project}.SV32.matrix.tsv"
    Analyze.cosmic_fit(
        samples=str(matrix_path),
        output=str(output_root / "SV32"),
        input_type="matrix",
        context_type="32",
        cosmic_version=cosmic_version,
        exome=False,
        genome_build=genome_build,
        signature_database=None,
        exclude_signature_subgroups=None,
        collapse_to_SBS96=False,
        export_probabilities=False,
        export_probabilities_per_mutation=False,
        make_plots=make_plots,
        sample_reconstruction_plots="none",
        verbose=True,
        cpu=cpu,
    )


def main():
    args = parse_args()
    sarek_root = Path(args.sarek_root)
    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=True)

    if not sarek_root.is_dir():
        raise SystemExit(f"Sarek root not found: {sarek_root}")

    if not (args.skip_completed and is_completed(output_root, "CNV48")):
        run_cnv_assignment(
            sarek_root=sarek_root,
            output_root=output_root,
            project=args.project,
            genome_build=args.genome_build,
            cosmic_version=args.cosmic_version,
            cpu=args.cpu,
            make_plots=args.make_plots,
        )

    if not (args.skip_completed and is_completed(output_root, "SV32")):
        run_sv_assignment(
            sarek_root=sarek_root,
            output_root=output_root,
            project=args.project,
            genome_build=args.genome_build,
            cosmic_version=args.cosmic_version,
            cpu=args.cpu,
            make_plots=args.make_plots,
        )


if __name__ == "__main__":
    main()
