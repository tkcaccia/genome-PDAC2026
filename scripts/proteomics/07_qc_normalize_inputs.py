#!/usr/bin/env python3
# PDAC2026 analysis: Protein-group abundance assembly and integration/QC normalization.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""QC and normalize the integration-ready matrices.

Inputs are read from `results/integration/` and cleaned outputs are written to
`results/qc/`.

This step prepares feature-filtered, analysis-ready matrices for MOFA+:
- WES: binary mutation matrix, filtered to recurrent genes
- RNA-seq: delta matrix, row-wise z-scored after removing zero-variance genes
- Proteomics: log2-abundance matrix, filtered by detection rate and row-wise z-scored
"""

from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

from pathlib import Path
import pandas as pd
import numpy as np


ROOT = project_root()
INTEGRATION = ROOT / 'results' / 'integration'
OUT = ROOT / 'results' / 'qc'
OUT.mkdir(parents=True, exist_ok=True)


def row_zscore(df: pd.DataFrame) -> pd.DataFrame:
    means = df.mean(axis=1)
    sds = df.std(axis=1).replace(0, np.nan)
    return df.sub(means, axis=0).div(sds, axis=0).fillna(0.0)


def summarize_matrix(name: str, df: pd.DataFrame) -> dict[str, float]:
    nonzero = (df != 0).mean(axis=1)
    variance = df.var(axis=1)
    return {
        'view': name,
        'n_features': float(df.shape[0]),
        'n_samples': float(df.shape[1]),
        'mean_nonzero_prevalence': float(nonzero.mean()),
        'median_nonzero_prevalence': float(nonzero.median()),
        'mean_variance': float(variance.mean()),
        'median_variance': float(variance.median()),
    }


def main() -> None:
    wes = pd.read_csv(INTEGRATION / 'wes_binary.csv', index_col=0).astype(float)
    rna = pd.read_csv(INTEGRATION / 'rnaseq_delta.csv', index_col=0).astype(float)
    prot = pd.read_csv(INTEGRATION / 'proteomics_log2.csv', index_col=0).astype(float)
    samples = pd.read_csv(INTEGRATION / 'samples.tsv', sep='\t')['patient'].tolist()

    # keep matrices aligned to the same sample order
    wes = wes[samples]
    rna = rna[samples]
    prot = prot[samples]

    # feature filtering
    wes_prevalence = (wes > 0).mean(axis=1)
    wes_keep = wes_prevalence >= 0.10  # recurrent mutations present in at least 10% of samples
    wes_f = wes.loc[wes_keep].copy()

    rna_variance = rna.var(axis=1)
    rna_keep = rna_variance > 0
    rna_f = rna.loc[rna_keep].copy()

    prot_detection = (prot > 0).mean(axis=1)
    prot_keep = prot_detection >= 0.20  # observed in at least 20% of samples
    prot_f = prot.loc[prot_keep].copy()

    # normalization
    rna_z = row_zscore(rna_f)
    prot_z = row_zscore(prot_f)

    # WES stays binary for downstream integration
    wes_qc = wes_f.fillna(0).astype(int)

    # write cleaned matrices
    wes_qc.to_csv(OUT / 'wes_binary_qc.csv')
    rna_z.to_csv(OUT / 'rnaseq_delta_zscore.csv')
    prot_z.to_csv(OUT / 'proteomics_log2_zscore.csv')

    # QC summary tables
    summary = pd.DataFrame([
        summarize_matrix('wes', wes),
        summarize_matrix('rna', rna),
        summarize_matrix('proteomics', prot),
    ])
    summary.to_csv(OUT / 'qc_summary.csv', index=False)

    filters = pd.DataFrame([
        {
            'view': 'wes',
            'rule': 'prevalence >= 0.10',
            'n_retained': int(wes_qc.shape[0]),
            'n_removed': int(wes.shape[0] - wes_qc.shape[0]),
        },
        {
            'view': 'rna',
            'rule': 'variance > 0',
            'n_retained': int(rna_z.shape[0]),
            'n_removed': int(rna.shape[0] - rna_z.shape[0]),
        },
        {
            'view': 'proteomics',
            'rule': 'detection >= 0.20',
            'n_retained': int(prot_z.shape[0]),
            'n_removed': int(prot.shape[0] - prot_z.shape[0]),
        },
    ])
    filters.to_csv(OUT / 'feature_filters.csv', index=False)

    report = [
        '# QC and normalization report',
        '',
        f'- Samples aligned: {len(samples)}',
        f'- WES features retained: {wes_qc.shape[0]} / {wes.shape[0]}',
        f'- RNA features retained: {rna_z.shape[0]} / {rna.shape[0]}',
        f'- Proteomics features retained: {prot_z.shape[0]} / {prot.shape[0]}',
        '',
        '## Outputs',
        f'- `{(OUT / "wes_binary_qc.csv").relative_to(ROOT)}`',
        f'- `{(OUT / "rnaseq_delta_zscore.csv").relative_to(ROOT)}`',
        f'- `{(OUT / "proteomics_log2_zscore.csv").relative_to(ROOT)}`',
        f'- `{(OUT / "qc_summary.csv").relative_to(ROOT)}`',
        f'- `{(OUT / "feature_filters.csv").relative_to(ROOT)}`',
    ]
    (OUT / 'QC_REPORT.md').write_text('\n'.join(report) + '\n')

    print('QC complete')
    print(summary.to_string(index=False))
    print(filters.to_string(index=False))


if __name__ == '__main__':
    main()
