"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""


from __future__ import annotations

from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path


from pathlib import Path


import matplotlib


matplotlib.use("Agg")


import matplotlib.pyplot as plt


from matplotlib.patches import Circle


import numpy as np


import pandas as pd


from sklearn.manifold import TSNE


ROOT = private_path("handoff_root")


RUN = ROOT / "mofa" / "pdac_11"


WES = ROOT / "wes" / "analysis" / "outputs" / "primary_pdac" / "wes_mutation_matrix.blacklist_filtered.csv"


FORMAL_DE = project_root() / "results" / "rnaseq_revised_11_native_deseq2_20260917" / "native_deseq2_degs_FDR05_abslog2FC1_11_pdac.csv"


FIG = RUN / "figures"


TABLE = RUN / "tables"


def make_formal_de_mutation_venn() -> None:
    wes = pd.read_csv(WES, index_col=0).fillna(0)
    de = pd.read_csv(FORMAL_DE)
    formal_de = set(de["gene_symbol"].dropna().astype(str))
    mutated = set(wes.index[wes.sum(axis=1) > 0])
    overlap = formal_de & mutated
    categories = {
        "formal_de_only": formal_de - overlap,
        "overlap": overlap,
        "mutated_only": mutated - overlap,
    }
    pd.DataFrame({
        "category": list(categories),
        "count": [len(genes) for genes in categories.values()],
        "definition": [
            "Paired native R/DESeq2 DEG (FDR < 0.05 and |log2 fold change| >= 1); not mutated by WES",
            "Paired native R/DESeq2 DEG and WES mutation, each present in the revised 11-patient PDAC cohort",
            "WES mutation; not a paired native R/DESeq2 DEG at FDR < 0.05 and |log2 fold change| >= 1",
        ],
    }).to_csv(TABLE / "formal_de_mutation_venn_counts_11_pdac.csv", index=False)
    for category, genes in categories.items():
        pd.DataFrame({"gene": sorted(genes)}).to_csv(TABLE / f"formal_de_mutation_venn_{category}_genes_11_pdac.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.set(xlim=(0, 7), ylim=(0, 4.8))
    ax.add_patch(Circle((2.65, 2.45), 1.8, color="#457b9d", alpha=0.5))
    ax.add_patch(Circle((4.35, 2.45), 1.8, color="#e9c46a", alpha=0.5))
    ax.text(1.55, 2.45, f"DEG only\n{len(categories['formal_de_only']):,}", ha="center", va="center", fontsize=11)
    ax.text(3.5, 2.45, f"Overlap\n{len(overlap):,}", ha="center", va="center", fontsize=11)
    ax.text(5.45, 2.45, f"Mutated only\n{len(categories['mutated_only']):,}", ha="center", va="center", fontsize=11)
    ax.text(2.3, 4.32, "Paired native DESeq2 DEGs", ha="center", fontsize=10)
    ax.text(4.7, 4.32, "Mutated genes", ha="center", fontsize=10)
    ax.set_axis_off()
    ax.set_title("Formal DE-mutation overlap: revised 11-patient PDAC cohort")
    fig.tight_layout()
    fig.savefig(FIG / "Figure_formal_DE_mutation_venn_11_pdac.png", dpi=300)
    plt.close(fig)



if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    TABLE.mkdir(parents=True, exist_ok=True)
    make_formal_de_mutation_venn()
