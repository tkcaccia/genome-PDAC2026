#!/usr/bin/env python3
# Current collaborator code export. Public-cohort preparation and descriptive comparisons; not ancestry-controlled validation.
# Source: Original Source and Code/scripts/19_prepare_tcga_paad_matrices.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Prepare TCGA PAAD matrices from downloaded Xena files.

Creates comparison-ready matrices in `external/tcga_paad/processed/`:
- `tcga_rna_matrix.csv`
- `tcga_rppa_matrix.csv`
- `tcga_mc3_gene_binary.csv`
"""

from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

from pathlib import Path
import pandas as pd


ROOT = project_root()
RAW = ROOT / 'external' / 'tcga_paad' / 'raw'
OUT = ROOT / 'external' / 'tcga_paad' / 'processed'
OUT.mkdir(parents=True, exist_ok=True)


def keep_primary_tumors(df: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in df.columns if str(c).startswith('TCGA-')]
    # primary tumor typically sample type code 01 in TCGA barcode
    keep = []
    for c in cols:
        parts = c.split('-')
        if len(parts) >= 4 and parts[3][:2] == '01':
            keep.append(c)
    return df[keep]


def main() -> None:
    rna = pd.read_csv(RAW / 'TCGA.PAAD.sampleMap_HiSeqV2.gz', sep='\t')
    rna = rna.rename(columns={rna.columns[0]: 'feature'})
    rna = rna.set_index('feature')
    rna_t = keep_primary_tumors(rna)
    rna_t.to_csv(OUT / 'tcga_rna_matrix.csv')

    rppa = pd.read_csv(RAW / 'TCGA.PAAD.sampleMap_RPPA.gz', sep='\t')
    rppa = rppa.rename(columns={rppa.columns[0]: 'feature'})
    rppa = rppa.set_index('feature')
    rppa_t = keep_primary_tumors(rppa)
    rppa_t.to_csv(OUT / 'tcga_rppa_matrix.csv')

    mc3 = pd.read_csv(RAW / 'mc3_gene_level_PAAD_mc3_gene_level.txt.gz', sep='\t')
    mc3 = mc3.rename(columns={mc3.columns[0]: 'feature'})
    mc3 = mc3.set_index('feature')
    mc3_t = keep_primary_tumors(mc3)
    mc3_t.to_csv(OUT / 'tcga_mc3_gene_binary.csv')

    summary = pd.DataFrame([
        {'matrix': 'rna', 'n_features': rna_t.shape[0], 'n_samples': rna_t.shape[1]},
        {'matrix': 'rppa', 'n_features': rppa_t.shape[0], 'n_samples': rppa_t.shape[1]},
        {'matrix': 'mc3_gene_binary', 'n_features': mc3_t.shape[0], 'n_samples': mc3_t.shape[1]},
    ])
    summary.to_csv(OUT / 'tcga_matrix_summary.csv', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
