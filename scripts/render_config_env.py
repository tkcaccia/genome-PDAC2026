#!/usr/bin/env python3
"""Render the supported YAML configuration keys as shell-safe variables."""

from __future__ import annotations

import argparse
import shlex
from pathlib import Path

import yaml


KEYS = {
    "OUTPUT_ROOT": ("project", "output_root"),
    "LOG_DIR": ("project", "log_dir"),
    "STAR_ANALYSIS_DIR": ("inputs", "star_analysis_dir"),
    "GTF": ("inputs", "gtf"),
    "METADATA": ("inputs", "metadata"),
    "PHENOTYPE_ASSIGNMENT": ("inputs", "phenotype_assignment"),
    "TME_SCORE_TABLE": ("inputs", "tme_score_table"),
    "PROGRAMME_GMT": ("inputs", "programme_gmt"),
    "CIBERSORT_SCRIPT": ("inputs", "cibersort_script"),
    "CIBERSORT_LM22": ("inputs", "cibersort_lm22"),
    "MEMORY_GATE": ("runtime", "memory_gate"),
    "MINIMUM_AVAILABLE_GB": ("runtime", "minimum_available_gb"),
    "MAXIMUM_SWAP_PERCENT": ("runtime", "maximum_swap_percent"),
    "RUN_STANDARD_DE": ("steps", "standard_de"),
    "RUN_GENE_TPM": ("steps", "gene_tpm"),
    "RUN_PATHWAY_SCORING": ("steps", "pathway_scoring"),
    "RUN_IMMUNE_DECONVOLUTION": ("steps", "immune_deconvolution"),
    "RUN_PHENOTYPE_ASSIGNMENT": ("steps", "phenotype_assignment"),
    "IMMUNE_METHODS": ("immune", "methods"),
    "CIBERSORT_PERMUTATIONS": ("immune", "cibersort_permutations"),
    "PHENOTYPE_IMMUNE_COLUMNS": ("phenotype", "immune_columns"),
    "PHENOTYPE_STROMAL_COLUMNS": ("phenotype", "stromal_columns"),
    "PHENOTYPE_EMT_COLUMNS": ("phenotype", "emt_columns"),
    "PHENOTYPE_TARGET_PER_EXTREME": ("phenotype", "target_per_extreme"),
}


def get_nested(config: dict, path: tuple[str, ...]):
    value = config
    for part in path:
        if not isinstance(value, dict) or part not in value:
            raise KeyError(".".join(path))
        value = value[part]
    return value


def shell_value(value: object) -> str:
    if isinstance(value, bool):
        value = "true" if value else "false"
    elif value is None:
        value = ""
    return shlex.quote(str(value))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise SystemExit("Configuration root must be a YAML mapping")

    missing = []
    for env_name, path in KEYS.items():
        try:
            value = get_nested(config, path)
        except KeyError:
            missing.append(".".join(path))
            continue
        print(f"{env_name}={shell_value(value)}")
    if missing:
        raise SystemExit("Missing configuration keys: " + ", ".join(missing))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
