#!/usr/bin/env python3
# PDAC2026 analysis: VEP-annotated DeepSomatic calls, retained-variant summaries and paired MSI launcher.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Build a WES gene x patient mutation matrix from per-patient VEP-annotated
DeepSomatic VCFs, with AD/VAF/DP quality filtering and mutation-count plots.

Each input VCF represents one patient and contains a tumour sample column.
Matched-normal evidence can occur in FORMAT fields NDP/NAD/NAF.
Somatic call status is read from FILTER; no genotype-based subtraction
is performed by this downstream parser.

Inputs:
- A directory of per-patient VEP-annotated VCFs, one per patient, e.g.:
    vep_annotated/patient_<id>.somatic.vep.vcf.gz
    vep_annotated/patient_<id>.somatic.vep.vcf.gz
    ...
  (VEP annotation must be run first — this script does not annotate.)

Outputs:
- results/wes/wes_somatic_variants.tsv
- results/wes/wes_mutation_matrix.csv
- results/wes/wes_mutation_burden.csv
- results/wes/wes_oncoprint.png
- results/wes/wes_mutation_burden.png
"""

from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

import argparse
import gzip
from pathlib import Path

import pandas as pd


PROJECT_ROOT = project_root()
DEFAULT_VCF_DIR = PROJECT_ROOT / "vep_annotated"
DEFAULT_OUTDIR = PROJECT_ROOT / "results" / "wes"
DEFAULT_VCF_GLOB = "patient_*.vep.vcf.gz"

# Variant-level depth and allele-support thresholds for matrix construction.
DEFAULT_MIN_DP = 10
DEFAULT_MIN_ALT_READS = 3
DEFAULT_MIN_VAF = 0.05


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf-dir", type=Path, default=DEFAULT_VCF_DIR,
                         help="Directory containing one VEP-annotated VCF per patient")
    parser.add_argument("--vcf-glob", type=str, default=DEFAULT_VCF_GLOB,
                         help="Glob pattern to match per-patient VCFs within --vcf-dir")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--min-dp", type=int, default=DEFAULT_MIN_DP,
                         help="Minimum tumour read depth (DP) to keep a variant")
    parser.add_argument("--min-alt-reads", type=int, default=DEFAULT_MIN_ALT_READS,
                         help="Minimum alt-supporting reads (from AD) to keep a variant")
    parser.add_argument("--min-vaf", type=float, default=DEFAULT_MIN_VAF,
                         help="Minimum variant allele fraction (VAF) to keep a variant")
    return parser.parse_args()


def parse_info_field(info_text: str) -> dict[str, str]:
    info: dict[str, str] = {}
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = "True"
    return info


def parse_format_values(format_keys: list[str], sample_field: str) -> dict[str, str]:
    values = sample_field.split(":")
    return dict(zip(format_keys, values))


def passes_quality_filters(fmt: dict[str, str], min_dp: int, min_alt_reads: int, min_vaf: float) -> bool:
    """Apply depth, alternate-read and variant-allele-fraction criteria."""
    dp_str = fmt.get("DP", "")
    if dp_str in ("", "."):
        return False
    try:
        dp = int(dp_str)
    except ValueError:
        return False
    if dp < min_dp:
        return False

    vaf_str = fmt.get("VAF", "")
    if vaf_str not in ("", "."):
        try:
            vaf = float(vaf_str)
            if vaf < min_vaf:
                return False
        except ValueError:
            pass  # fall through to AD-based check if VAF unparsable

    ad_str = fmt.get("AD", "")
    if ad_str not in ("", "."):
        parts = ad_str.split(",")
        if len(parts) >= 2:
            try:
                alt_reads = int(parts[1])  # AD = ref,alt for biallelic records
                if alt_reads < min_alt_reads:
                    return False
            except ValueError:
                pass

    return True


def iter_patient_vcfs(vcf_dir: Path, pattern: str):
    files = sorted(vcf_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No VCFs matched {pattern} in {vcf_dir}")
    for f in files:
        yield f


def extract_patient_id(sample_name: str) -> str:
    # Strip the tissue suffix from the numeric patient identifier.
    return sample_name[:-1] if sample_name and sample_name[-1] in ("N", "T") else sample_name


def process_vcf(vcf_path: Path, min_dp: int, min_alt_reads: int, min_vaf: float) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    opener = gzip.open if vcf_path.suffix == ".gz" else open

    with opener(vcf_path, "rt") as handle:
        sample_name = None
        for line in handle:
            if line.startswith("#CHROM"):
                cols = line.rstrip().split("\t")
                sample_cols = cols[9:]
                if len(sample_cols) != 1:
                    raise ValueError(
                        f"{vcf_path.name}: expected exactly 1 sample column (tumour only), "
                        f"found {len(sample_cols)}: {sample_cols}"
                    )
                sample_name = sample_cols[0]
                break
        if sample_name is None:
            raise ValueError(f"{vcf_path.name}: could not find #CHROM header line")

        patient_id = extract_patient_id(sample_name)

        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 10:
                continue

            chrom, pos, _id, ref, alt, _qual, filt = fields[:7]
            if filt != "PASS":
                continue  # DeepSomatic already classifies somatic status via FILTER

            info = parse_info_field(fields[7])
            # VEP emits one entry per overlapping transcript, comma-joined
            # (e.g. "MODERATE,MODERATE,MODERATE" for a variant hitting 3
            # transcripts) — so these fields must be split per-transcript
            # before filtering, not compared as a single exact string.
            impacts = info.get("CSQ_IMPACT", info.get("IMPACT", "")).strip().upper().split(",")
            symbols = info.get("CSQ_SYMBOL", info.get("SYMBOL", "")).strip().split(",")
            consequences = info.get("CSQ_Consequence", info.get("Consequence", "")).split(",")
            sifts = info.get("CSQ_SIFT", info.get("SIFT", "")).split(",")
            polyphens = info.get("CSQ_PolyPhen", info.get("PolyPhen", "")).split(",")

            impact_rank = {"HIGH": 0, "MODERATE": 1}
            best_idx = None
            for i, imp in enumerate(impacts):
                if imp in impact_rank and (
                    best_idx is None or impact_rank[imp] < impact_rank[impacts[best_idx]]
                ):
                    best_idx = i
            if best_idx is None:
                continue  # no transcript for this variant reaches HIGH/MODERATE

            impact = impacts[best_idx]
            symbol = symbols[best_idx].strip() if best_idx < len(symbols) else ""
            if symbol in {"", "."}:
                continue

            format_keys = fields[8].split(":")
            fmt = parse_format_values(format_keys, fields[9])

            if not passes_quality_filters(fmt, min_dp, min_alt_reads, min_vaf):
                continue

            consequence = consequences[best_idx] if best_idx < len(consequences) else ""
            sift = sifts[best_idx] if best_idx < len(sifts) else ""
            polyphen = polyphens[best_idx] if best_idx < len(polyphens) else ""

            rows.append({
                "patient": patient_id,
                "sample": sample_name,
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "SYMBOL": symbol,
                "Consequence": consequence,
                "IMPACT": impact,
                "SIFT": sift,
                "PolyPhen": polyphen,
                "DP": fmt.get("DP", ""),
                "VAF": fmt.get("VAF", ""),
                "AD": fmt.get("AD", ""),
            })

    return rows


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, object]] = []
    for vcf_path in iter_patient_vcfs(args.vcf_dir, args.vcf_glob):
        print(f"Processing {vcf_path.name}...")
        rows = process_vcf(vcf_path, args.min_dp, args.min_alt_reads, args.min_vaf)
        print(f"  {len(rows)} PASS variants retained after impact + quality filters")
        all_rows.extend(rows)

    somatic_df = pd.DataFrame(all_rows)
    if somatic_df.empty:
        raise ValueError("No somatic variants were identified across any patient VCF")

    somatic_df["patient"] = somatic_df["patient"].astype(str)
    somatic_df["patient_id"] = somatic_df["patient"].str.zfill(2)
    somatic_df["tumour_sample"] = somatic_df["patient_id"] + "T"

    mutation_matrix = (
        somatic_df.groupby(["SYMBOL", "tumour_sample"])
        .size()
        .unstack(fill_value=0)
        .clip(upper=1)
    )
    mutation_matrix = mutation_matrix.reindex(sorted(mutation_matrix.columns), axis=1)

    burden = (
        somatic_df.groupby("tumour_sample")
        .size()
        .reset_index(name="somatic_variant_count")
        .sort_values("somatic_variant_count", ascending=False)
    )

    somatic_df = somatic_df.drop(columns=["patient_id", "tumour_sample"])

    somatic_df.to_csv(args.outdir / "wes_somatic_variants.tsv", sep="\t", index=False)
    mutation_matrix.to_csv(args.outdir / "wes_mutation_matrix.csv")
    burden.to_csv(args.outdir / "wes_mutation_burden.csv", index=False)

    print(f"\nSaved somatic variants: {args.outdir / 'wes_somatic_variants.tsv'}")
    print(f"Saved mutation matrix:  {args.outdir / 'wes_mutation_matrix.csv'}")
    print(f"Saved burden table:     {args.outdir / 'wes_mutation_burden.csv'}")
    print(f"\nSomatic variants (post-filter): {len(somatic_df)}")
    print(f"Matrix shape (genes x patients): {mutation_matrix.shape}")


if __name__ == "__main__":
    main()
