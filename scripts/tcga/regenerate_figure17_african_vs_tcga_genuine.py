#!/usr/bin/env python3
# PDAC2026 analysis: Public-cohort preparation and descriptive comparisons; not ancestry-controlled validation.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""PDAC versus TCGA PAAD MC3 mutation-associated pathway-theme comparison.

Both cohorts use somatic mutation calls and the recurrence rule of
mutation in at least two samples. Enrichment results are aggregated
into broad pathway themes for descriptive comparison.

Outputs:
- results/mofa_revised_11_20260916/tables/revised_11_vs_tcga_pathway_theme_concordance.csv
- results/mofa_revised_11_20260916/tables/revised_11_wes_recurrent_gene_enrichment.csv
- results/mofa_revised_11_20260916/tables/tcga_paad_mc3_recurrent_gene_enrichment.csv
- results/figures/final_11_figures/Figure_17_pathway_theme_concordance.png
"""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

ROOT = project_root()
AFRICAN_GENE_SUMMARY = ROOT / "results/wes_revised_11_20260916/outputs/primary_pdac/gene_summary_table.csv"
TCGA_MC3 = ROOT / "external/tcga_paad/processed/tcga_mc3_gene_binary.csv"
TABLES = ROOT / "results/mofa_revised_11_20260916/tables"
OUT_FIG = ROOT / "results/figures/final_11_figures/Figure_17_pathway_theme_concordance.png"
RECURRENCE_MIN = 2
GENE_SETS = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]


def theme_for_term(term: str) -> str:
    text = term.lower()
    rules = [
        ("DNA repair / replication", ["dna", "repair", "replication", "recombination", "nhej", "mismatch", "telomere"]),
        ("Translation / ribosome", ["translation", "ribosome", "ribosomal", "translational", "elongation"]),
        ("Metabolism / lipid", ["metab", "lipid", "cholesterol", "fatty", "glycol", "sphing", "ceramide", "sterol"]),
        ("Immune / inflammation", ["immune", "interleukin", "cytokine", "inflamm", "leukocyte", "t cell", "b cell", "antigen"]),
        ("ECM / adhesion", ["extracellular matrix", "ecm", "adhesion", "collagen", "integrin", "fibronectin", "basement membrane"]),
        ("Cytoskeleton / motility", ["actin", "microtubule", "cytoskeleton", "motility", "migration", "myosin", "tubulin"]),
        ("Mitochondria / stress", ["mitochond", "oxidative", "stress", "ros", "respiratory", "electron transport"]),
        ("Vesicle / trafficking", ["vesicle", "traffick", "endocyt", "exocyt", "golgi", "lysosome", "autophagy"]),
        ("Cell cycle / signaling", ["cell cycle", "mitotic", "g2/m", "checkpoint", "erk", "mapk", "akt", "mtor", "signaling"]),
    ]
    for label, keywords in rules:
        if any(keyword in text for keyword in keywords):
            return label
    return "Other / mixed"


def enrich(genes: list[str]) -> pd.DataFrame:
    import gseapy as gp

    rows = []
    for gene_set in GENE_SETS:
        try:
            result = gp.enrichr(gene_list=genes, gene_sets=gene_set, organism="human", outdir=None, cutoff=1.0)
            if result is None or result.results is None or len(result.results) == 0:
                continue
            table = result.results.copy()
            table.insert(0, "gene_set", gene_set)
            rows.append(table)
        except Exception as error:
            print(f"Skipping {gene_set}: {error}")
    if rows:
        return pd.concat(rows, ignore_index=True)
    return pd.DataFrame(columns=["gene_set", "Term", "Adjusted P-value", "Combined Score", "Genes"])


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)

    african_summary = pd.read_csv(AFRICAN_GENE_SUMMARY)
    african_genes = sorted(african_summary.loc[african_summary["patients_affected"] >= RECURRENCE_MIN, "SYMBOL"].dropna().astype(str).unique())

    tcga = pd.read_csv(TCGA_MC3, index_col=0).fillna(0)
    tcga_freq = (tcga != 0).sum(axis=1)
    tcga_genes = sorted(tcga_freq[tcga_freq >= RECURRENCE_MIN].index.astype(str))

    print(f"African recurrent genes (>= {RECURRENCE_MIN}/11 patients): {len(african_genes)}")
    print(f"TCGA PAAD MC3 recurrent genes (>= {RECURRENCE_MIN}/{tcga.shape[1]} samples): {len(tcga_genes)}")

    african_enrichment = enrich(african_genes)
    african_enrichment.to_csv(TABLES / "revised_11_wes_recurrent_gene_enrichment.csv", index=False)
    tcga_enrichment = enrich(tcga_genes)
    tcga_enrichment.to_csv(TABLES / "tcga_paad_mc3_recurrent_gene_enrichment.csv", index=False)

    african_theme_counts = Counter(theme_for_term(term) for term in african_enrichment.get("Term", pd.Series(dtype=str)).dropna().astype(str))
    tcga_theme_counts = Counter(theme_for_term(term) for term in tcga_enrichment.get("Term", pd.Series(dtype=str)).dropna().astype(str))

    theme_rows = []
    for theme in sorted(set(african_theme_counts) | set(tcga_theme_counts)):
        theme_rows.append({
            "theme": theme,
            "african_terms": african_theme_counts.get(theme, 0),
            "tcga_terms": tcga_theme_counts.get(theme, 0),
        })
    theme_df = pd.DataFrame(theme_rows).sort_values(["african_terms", "tcga_terms"], ascending=False)
    theme_df.to_csv(TABLES / "revised_11_vs_tcga_pathway_theme_concordance.csv", index=False)

    african_presence = (theme_df["african_terms"] > 0).astype(int).values
    tcga_presence = (theme_df["tcga_terms"] > 0).astype(int).values
    union = int((african_presence | tcga_presence).sum())
    intersection = int((african_presence & tcga_presence).sum())
    jaccard = intersection / union if union else 0.0
    rho, pvalue = spearmanr(theme_df["african_terms"], theme_df["tcga_terms"])

    plt.figure(figsize=(6, 6))
    themes = theme_df["theme"].tolist()
    palette = dict(zip(themes, sns.color_palette("tab10", n_colors=len(themes))))
    for theme in themes:
        row = theme_df[theme_df["theme"] == theme]
        plt.scatter(row["african_terms"], row["tcga_terms"], s=100, color=palette[theme], edgecolor="k",
                    label=f"{theme} ({int(row['african_terms'].iloc[0])}/{int(row['tcga_terms'].iloc[0])})")
    plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8)
    min_value = min(theme_df["african_terms"].min(), theme_df["tcga_terms"].min())
    max_value = max(theme_df["african_terms"].max(), theme_df["tcga_terms"].max())
    plt.plot([min_value - 5, max_value + 5], [min_value - 5, max_value + 5], "--", color="gray", linewidth=0.8)
    plt.xlabel("African pathway term count (aggregated)")
    plt.ylabel("TCGA pathway term count (aggregated)")
    plt.title(f"Pathway theme concordance: revised 11-patient PDAC vs TCGA PAAD\nRecurrently mutated genes (>= {RECURRENCE_MIN} samples); "
              f"Jaccard = {jaccard:.2f}; Spearman rho = {rho:.3f} (p = {pvalue:.3g})")
    plt.tight_layout()
    plt.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {OUT_FIG}")
    print(theme_df.to_string(index=False))


if __name__ == "__main__":
    main()
