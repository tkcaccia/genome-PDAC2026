#!/usr/bin/env python3
import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run COSMIC mutational signature assignment on a folder of PASS-only somatic VCFs."
    )
    parser.add_argument("--samples", required=True, help="Directory containing PASS-only VCFs")
    parser.add_argument("--output", required=True, help="Output directory root")
    parser.add_argument("--genome-build", default="GRCh38")
    parser.add_argument("--cosmic-version", type=float, default=3.5)
    parser.add_argument("--cpu", type=int, default=8)
    parser.add_argument(
        "--contexts",
        nargs="+",
        choices=("96", "DINUC", "ID"),
        default=("96", "DINUC", "ID"),
        help="One or more SigProfiler context types to run.",
    )
    parser.add_argument(
        "--skip-completed",
        action="store_true",
        help="Skip contexts whose JOB_METADATA_SPA.txt already reports successful completion.",
    )
    parser.add_argument("--make-plots", action="store_true")
    return parser.parse_args()


def run_context(samples: str, output_root: Path, context: str, cosmic_version: float, genome_build: str, cpu: int, make_plots: bool):
    try:
        from SigProfilerAssignment import Analyzer as Analyze
    except ImportError as exc:
        raise SystemExit(
            "SigProfilerAssignment is required to run signature assignment. "
            "Install the environment documented in env/environment.yml."
        ) from exc

    context_name = {
        "96": "SBS96",
        "DINUC": "DBS",
        "ID": "ID",
    }[context]
    output_dir = output_root / context_name
    output_dir.mkdir(parents=True, exist_ok=True)

    Analyze.cosmic_fit(
        samples=samples,
        output=str(output_dir),
        input_type="vcf",
        context_type=context,
        cosmic_version=cosmic_version,
        exome=True,
        genome_build=genome_build,
        signature_database=None,
        exclude_signature_subgroups=None,
        collapse_to_SBS96=(context == "96"),
        export_probabilities=False,
        export_probabilities_per_mutation=False,
        make_plots=make_plots,
        sample_reconstruction_plots="none",
        verbose=True,
        cpu=cpu,
    )


def is_completed(output_root: Path, context: str) -> bool:
    context_name = {
        "96": "SBS96",
        "DINUC": "DBS",
        "ID": "ID",
    }[context]
    metadata = output_root / context_name / "JOB_METADATA_SPA.txt"
    if not metadata.is_file():
        return False
    return "completed successfully" in metadata.read_text(errors="ignore").lower()


def main():
    args = parse_args()
    samples = Path(args.samples)
    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=True)

    if not samples.is_dir():
        raise SystemExit(f"Samples directory not found: {samples}")

    for context in args.contexts:
        if args.skip_completed and is_completed(output_root, context):
            print(f"Skipping completed context: {context}")
            continue
        run_context(
            samples=str(samples),
            output_root=output_root,
            context=context,
            cosmic_version=args.cosmic_version,
            genome_build=args.genome_build,
            cpu=args.cpu,
            make_plots=args.make_plots,
        )


if __name__ == "__main__":
    main()
