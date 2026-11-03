#!/usr/bin/env python3
# Current collaborator code export. Salmon gene counts/TPM, paired DESeq2, response PCA and enrichment.
# Source: Original Source and Code/scripts/generate_native_rna_pathway_table_11.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Generate exploratory GO/Reactome enrichment tables for native RNA DEGs."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path

import pandas as pd

ROOT = project_root()
DEG = ROOT / "results" / "rnaseq_revised_11_native_deseq2_20260917" / "native_deseq2_degs_FDR05_abslog2FC1_11_pdac.csv"
OUT = ROOT / "results" / "rnaseq_revised_11_native_deseq2_20260917"
LIBRARIES = ["GO_Biological_Process_2023", "Reactome_2022", "KEGG_2021_Human"]


def enrich(genes: list[str], direction: str) -> list[pd.DataFrame]:
    import gseapy as gp

    results = []
    for library in LIBRARIES:
        try:
            result = gp.enrichr(gene_list=genes, gene_sets=library, organism="human", outdir=None, no_plot=True)
            table = result.results.copy()
            table.insert(0, "direction", direction)
            table.insert(1, "gene_set", library)
            results.append(table)
        except Exception as error:
            print(f"Skipping {direction}/{library}: {error}")
    return results


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    deg = pd.read_csv(DEG)
    deg = deg.dropna(subset=["gene_symbol"])
    up = sorted(set(deg.loc[deg["direction"] == "up_in_tumour", "gene_symbol"].astype(str)))
    down = sorted(set(deg.loc[deg["direction"] == "down_in_tumour", "gene_symbol"].astype(str)))
    tables = enrich(up, "up_in_tumour") + enrich(down, "down_in_tumour")
    if tables:
        result = pd.concat(tables, ignore_index=True)
        result.to_csv(OUT / "native_rna_de_pathway_enrichment_11_pdac.csv", index=False)
    else:
        pd.DataFrame(columns=["direction", "gene_set", "Term", "P-value", "Adjusted P-value", "Genes"]).to_csv(OUT / "native_rna_de_pathway_enrichment_11_pdac.csv", index=False)
    pd.DataFrame({
        "metric": ["up_genes_submitted", "down_genes_submitted", "tested_genes_background", "background_note"],
        "value": [len(up), len(down), 15422, "Enrichr library background; custom tested-gene background unavailable locally"],
    }).to_csv(OUT / "native_rna_de_pathway_enrichment_11_pdac_summary.csv", index=False)


if __name__ == "__main__":
    main()