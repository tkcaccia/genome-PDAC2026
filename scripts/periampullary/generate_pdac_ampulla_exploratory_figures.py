#!/usr/bin/env python3
# PDAC2026 analysis: Separate exploratory PDAC/periampullary comparisons.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Generate exploratory PDAC versus periampullary pathway/KRAS comparison figures."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = project_root()
OUT = ROOT / "results/pdac_ampulla_comparison_20260917"
FIG = OUT / "figures"
TAB = OUT / "tables"
PDAC = configured_cohort("primary_pdac")
AMP = configured_cohort("periampullary")
LIBRARIES = ["GO_Biological_Process_2023", "Reactome_2022", "KEGG_2021_Human"]
RNA_MAP = ROOT / "reproducible_pipeline_share_20260713_with_outputs/raw_outputs/results/figures/categories/04_cross_cohort_concordance/rnaseq_ensembl_to_symbol.csv"
PROTEIN_MAP = ROOT / "reproducible_pipeline_share_20260713_with_outputs/raw_outputs/results/comparison_african_tcga/proteomics_index_to_symbol.csv"


def symbols(values, mapping_path):
    mapping = pd.read_csv(mapping_path, index_col=0)["symbol"]
    return sorted(set(str(mapping.get(str(value).split(".")[0], value)) for value in values if pd.notna(value)))


def enrich(genes, group, modality):
    import gseapy as gp
    rows = []
    for library in LIBRARIES:
        try:
            result = gp.enrichr(gene_list=genes, gene_sets=library, organism="human", outdir=None, cutoff=1.0)
            if result is None or result.results is None or result.results.empty:
                continue
            table = result.results.copy()
            table.insert(0, "group", group)
            table.insert(1, "modality", modality)
            table.insert(2, "gene_set", library)
            rows.append(table)
        except Exception:
            pass
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def build_gene_sets():
    wes = pd.read_csv(ROOT / "results/wes/wes_mutation_matrix.blacklist_filtered.csv", index_col=0).fillna(0)
    wes.columns = [str(c).lstrip("0").rstrip("T") for c in wes.columns]
    rna = pd.read_csv(ROOT / "results/rnaseq/rnaseq_delta_matrix.csv", index_col=0)
    rna.columns = [str(c).lstrip("0").rstrip("T") for c in rna.columns]
    prot = pd.read_csv(ROOT / "results/proteomics/protein_abundance_matrix.csv", index_col=0)
    prot_response = {}
    for patient in PDAC + AMP:
        prot_response[patient] = np.log2(pd.to_numeric(prot[f"{patient}T"], errors="coerce") + 1) - np.log2(pd.to_numeric(prot[f"{patient}N"], errors="coerce") + 1)
    prot_response = pd.DataFrame(prot_response, index=prot.index)
    sets = []
    for group, patients in [("PDAC", PDAC), ("Ampulla", AMP)]:
        mutated = sorted(wes[patients].index[wes[patients].sum(axis=1) >= 2].astype(str))
        sets.append((group, "WES", mutated))
        rna_effect = rna[patients].mean(axis=1).abs().nlargest(100).index
        sets.append((group, "RNA", symbols(rna_effect, RNA_MAP)))
        prot_effect = prot_response[patients].mean(axis=1).abs().nlargest(100).index
        sets.append((group, "Proteomics", symbols(prot_effect, PROTEIN_MAP)))
    return sets


def make_pathway_figures():
    results = []
    for group, modality, genes in build_gene_sets():
        if len(genes) >= 3:
            result = enrich(genes, group, modality)
            if not result.empty:
                results.append(result)
    enrichment = pd.concat(results, ignore_index=True)
    enrichment.to_csv(TAB / "pdac_ampulla_exploratory_pathway_enrichment.csv", index=False)
    enrichment["Adjusted P-value"] = pd.to_numeric(enrichment["Adjusted P-value"], errors="coerce")
    enrichment["minus_log10_fdr"] = -np.log10(enrichment["Adjusted P-value"].clip(lower=1e-300))
    enrichment["label"] = enrichment["Term"].astype(str).map(lambda x: textwrap.fill(x[:70], 36))

    top_shared = []
    for term in enrichment["Term"].dropna().unique():
        subset = enrichment[enrichment["Term"] == term]
        groups = set(zip(subset.group, subset.modality))
        if any(group == "PDAC" for group, _ in groups) and any(group == "Ampulla" for group, _ in groups):
            top_shared.append(subset.sort_values("Adjusted P-value").iloc[0])
    shared = pd.DataFrame(top_shared).sort_values("Adjusted P-value").head(20)
    shared.to_csv(TAB / "pdac_ampulla_shared_pathways.csv", index=False)

    if not shared.empty:
        plot = shared.sort_values("minus_log10_fdr")
        fig, ax = plt.subplots(figsize=(11, 8))
        colors = plot["modality"].map({"WES": "#264653", "RNA": "#e76f51", "Proteomics": "#2a9d8f"})
        ax.barh(plot["label"], plot["minus_log10_fdr"], color=colors)
        ax.set(xlabel="-log10 adjusted p-value", ylabel="Shared pathway term", title="Shared pathway mechanisms: PDAC versus ampulla\nExploratory enrichment across WES, RNA-seq and proteomics")
        handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in ["#264653", "#e76f51", "#2a9d8f"]]
        ax.legend(handles, ["WES", "RNA", "Proteomics"], title="Modality", frameon=False)
        ax.grid(axis="x", alpha=.25, linestyle=":")
        fig.tight_layout(); fig.savefig(FIG / "pdac_ampulla_shared_pathways.png", dpi=300, bbox_inches="tight"); plt.close(fig)

    enrichment["pathway_key"] = enrichment["Term"].astype(str)
    top_terms = enrichment.groupby("pathway_key")["minus_log10_fdr"].max().nlargest(20).index
    heat = enrichment[enrichment.pathway_key.isin(top_terms)].pivot_table(index="pathway_key", columns=["group", "modality"], values="minus_log10_fdr", aggfunc="max", fill_value=0)
    heat = heat.reindex(top_terms)
    heat.columns = [" / ".join(col) for col in heat.columns]
    heat.to_csv(TAB / "pdac_ampulla_pathway_comparison_heatmap.csv")
    fig, ax = plt.subplots(figsize=(12, 9)); sns.heatmap(heat, cmap="magma", linewidths=.4, linecolor="white", ax=ax, cbar_kws={"label":"-log10 adjusted p-value"}); ax.set_title("PDAC versus ampulla pathway enrichment\nTop exploratory pathway terms"); ax.set_xlabel(""); ax.set_ylabel("Pathway term"); fig.tight_layout(); fig.savefig(FIG / "pdac_ampulla_pathway_comparison_heatmap.png", dpi=300, bbox_inches="tight"); plt.close(fig)


def make_kras_figures():
    wes = pd.read_csv(ROOT / "results/wes/wes_mutation_matrix.blacklist_filtered.csv", index_col=0).fillna(0)
    wes.columns = [str(c).lstrip("0") for c in wes.columns]
    patients = PDAC + AMP
    status = pd.DataFrame(index=patients)
    status["group"] = ["PDAC"] * len(PDAC) + ["Ampulla"] * len(AMP)
    status["KRAS_mutated"] = [int(wes.get(f"{p}T", pd.Series(0, index=wes.index)).get("KRAS", 0) > 0) for p in patients]
    kras_rows = pd.read_csv(ROOT / "results/wes_revised_11_20260916/tables/kras_variants_all_cohorts.tsv", sep="\t")
    status["KRAS_change"] = ["; ".join(kras_rows.loc[kras_rows.patient.astype(str).eq(p), "ALT"].astype(str)) if p in set(kras_rows.patient.astype(str)) else "None" for p in patients]
    status.to_csv(TAB / "pdac_ampulla_kras_status.csv")
    fig, ax = plt.subplots(figsize=(9, 6))
    matrix = status[["KRAS_mutated"]].T
    sns.heatmap(matrix, cmap=["#eeeeee", "#d1495b"], cbar=False, linewidths=.6, linecolor="white", annot=False, ax=ax)
    ax.set_title("KRAS mutation status: PDAC versus ampulla\nExploratory comparison")
    ax.set_xlabel("Patient"); ax.set_ylabel("")
    ax.set_xticklabels([f"{p}\n{g}" for p, g in zip(patients, status.group)], rotation=0)
    fig.tight_layout(); fig.savefig(FIG / "pdac_ampulla_kras_status.png", dpi=300, bbox_inches="tight"); plt.close(fig)


def main():
    FIG.mkdir(parents=True, exist_ok=True); TAB.mkdir(parents=True, exist_ok=True)
    make_pathway_figures(); make_kras_figures(); print("Generated PDAC-ampulla exploratory pathway and KRAS figures")

if __name__ == "__main__": main()
