#!/usr/bin/env python3
# PDAC2026 analysis: Public-cohort preparation and descriptive comparisons; not ancestry-controlled validation.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Download TCGA PAAD core datasets from UCSC Xena.

Datasets pulled:
- RNA expression: `TCGA.PAAD.sampleMap/HiSeqV2.gz`
- Protein RPPA: `TCGA.PAAD.sampleMap/RPPA.gz`
- Somatic mutation (gene-level MC3): `mc3_gene_level/PAAD_mc3_gene_level.txt.gz`
- Somatic mutation (variant-level MC3): `mc3/PAAD_mc3.txt.gz`

Outputs are written under `external/tcga_paad/raw/`.
"""

from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

from pathlib import Path
import requests


ROOT = project_root()
OUT = ROOT / 'external' / 'tcga_paad' / 'raw'
OUT.mkdir(parents=True, exist_ok=True)

BASE = 'https://tcga.xenahubs.net/download/'
FILES = {
    'TCGA.PAAD.sampleMap_HiSeqV2.gz': 'TCGA.PAAD.sampleMap/HiSeqV2.gz',
    'TCGA.PAAD.sampleMap_RPPA.gz': 'TCGA.PAAD.sampleMap/RPPA.gz',
    'mc3_gene_level_PAAD_mc3_gene_level.txt.gz': 'mc3_gene_level/PAAD_mc3_gene_level.txt.gz',
    'mc3_PAAD_mc3.txt.gz': 'mc3/PAAD_mc3.txt.gz',
}


def download(url: str, target: Path) -> None:
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(target, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)


def main() -> None:
    rows = []
    for out_name, remote in FILES.items():
        url = BASE + remote
        target = OUT / out_name
        download(url, target)
        rows.append((out_name, target.stat().st_size))
        print('Downloaded', out_name, 'size_bytes=', target.stat().st_size)

    manifest = OUT / 'manifest.tsv'
    with open(manifest, 'w', encoding='utf-8') as f:
        f.write('file\tsize_bytes\n')
        for name, size in rows:
            f.write(f'{name}\t{size}\n')
    print('Wrote manifest:', manifest)


if __name__ == '__main__':
    main()
