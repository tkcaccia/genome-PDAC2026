#!/usr/bin/env python3
# PDAC2026 analysis: VEP-annotated DeepSomatic calls, retained-variant summaries and paired MSI launcher.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Extract genuine DeepSomatic-native KRAS hotspot fields for the 14-patient cohort.

DeepSomatic VCFs never emit a TLOD field (that is Mutect2-only); this script reads
only the fields DeepSomatic actually produces (QUAL, GT, GQ, DP, AD, VAF, NDP, NAD,
NAF) directly from the per-patient VEP-annotated VCFs, so the resulting table cannot
be accidentally populated with Mutect2 values.

Input:
    wes_somatic/vep_annotated/patient_*.vep.vcf.gz

Output:
    results/figures_deepsomatic_corrected_20260828/deepsomatic_kras_native_fields.csv
"""

from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

import gzip
import glob
import re
from pathlib import Path

import pandas as pd

ROOT = project_root()
VCF_DIR = ROOT / "wes_somatic" / "vep_annotated"
OUT_CSV = ROOT / "results" / "figures_deepsomatic_corrected_20260828" / "deepsomatic_kras_native_fields.csv"

# GRCh38 KRAS codon 12/13 (chr12:25245274-25245384) and codon 61 (chr12:25227336-25227348) hotspot window.
KRAS_CONTIG = "NC_000012.12"
KRAS_START = 25225000
KRAS_END = 25250929


def parse_kras_records(vcf_path: Path) -> list[dict]:
    patient = re.search(r"patient_(\w+)\.vep", vcf_path.name).group(1)
    rows = []
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, _id, ref, alt, qual, filt = fields[:7]
            if chrom != KRAS_CONTIG:
                continue
            pos_i = int(pos)
            if not (KRAS_START <= pos_i <= KRAS_END):
                continue

            info_text = fields[7]
            symbol = ""
            for item in info_text.split(";"):
                if item.startswith("CSQ_SYMBOL="):
                    symbol = item.split("=", 1)[1].split(",")[0]
                    break
            if symbol != "KRAS":
                continue

            format_keys = fields[8].split(":")
            sample_vals = fields[9].split(":")
            fmt = dict(zip(format_keys, sample_vals))

            rows.append({
                "patient": patient,
                "CHROM": chrom,
                "POS": pos_i,
                "REF": ref,
                "ALT": alt,
                "QUAL": qual,
                "FILTER": filt,
                "GT": fmt.get("GT", ""),
                "GQ": fmt.get("GQ", ""),
                "DP": fmt.get("DP", ""),
                "AD": fmt.get("AD", ""),
                "VAF": fmt.get("VAF", ""),
                "NDP": fmt.get("NDP", ""),
                "NAD": fmt.get("NAD", ""),
                "NAF": fmt.get("NAF", ""),
            })
    return rows


def main() -> None:
    all_rows = []
    for vcf_path in sorted(Path(p) for p in glob.glob(str(VCF_DIR / "patient_*.vep.vcf.gz"))):
        all_rows.extend(parse_kras_records(vcf_path))

    df = pd.DataFrame(all_rows)
    df["patient"] = df["patient"].astype(int)
    df = df.sort_values(["patient", "POS"]).reset_index(drop=True)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    print(f"Wrote {len(df)} KRAS-region records across {df['patient'].nunique()} patients to {OUT_CSV}")
    print()
    pass_only = df[df["FILTER"] == "PASS"]
    print(f"PASS KRAS records ({len(pass_only)}):")
    print(pass_only.to_string(index=False))


if __name__ == "__main__":
    main()
