#!/usr/bin/env python3
# WES processing, variant summaries, coverage and paired MSI.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""
from __future__ import annotations

from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path

import csv
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage

ROOT = private_path("handoff_root")
WORKSPACE_ROOT = project_root()
RUN = ROOT / "wes" / "analysis"
VCF_DIR = ROOT / "wes" / "annotated"
OUT = RUN / "outputs"
FIG = RUN / "figures"
TABLE = RUN / "tables"
PRIMARY = configured_cohort("primary_pdac")
PERI = configured_cohort("periampullary")
ALL_PATIENTS = PRIMARY + PERI
DRIVER_GENES = {
    "KRAS", "TP53", "SMAD4", "CDKN2A", "ARID1A", "KDM6A", "RNF43",
    "GNAS", "BRAF", "PIK3CA", "ATM", "BRCA1", "BRCA2", "PALB2",
    "ERBB2", "MYC", "STK11", "NF1", "TGFBR2", "FBXW7",
}
PATHWAYS = {
    "KRAS_MAPK": {"KRAS", "BRAF", "NF1", "MAP2K1", "MAP2K2", "RAF1", "DUSP6"},
    "TP53_DNA_damage": {"TP53", "ATM", "BRCA1", "BRCA2", "PALB2", "CHEK2", "ATR"},
    "TGF_BETA": {"SMAD4", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SKI"},
    "PI3K_AKT_MTOR": {"PIK3CA", "PIK3CB", "AKT1", "AKT2", "MTOR", "PTEN", "TSC1", "TSC2"},
    "WNT_BETA_CATENIN": {"RNF43", "APC", "CTNNB1", "AXIN1", "FBXW7", "GNAS"},
    "CHROMATIN_REMODELING": {"ARID1A", "KDM6A", "SMARCA4", "PBRM1", "BAP1", "KMT2C", "KMT2D"},
}


def load_fixed_parser():
    path = code_path("scripts/wes/wes_somatic_01_wes_build_mutation_matrix_FIXED.py")
    spec = importlib.util.spec_from_file_location("fixed_wes_parser", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_calls() -> pd.DataFrame:
    parser = load_fixed_parser()
    rows = []
    for patient in ALL_PATIENTS:
        path = VCF_DIR / f"patient_{patient}.vep.vcf.gz"
        if not path.exists():
            raise FileNotFoundError(path)
        rows.extend(parser.process_vcf(path, parser.DEFAULT_MIN_DP, parser.DEFAULT_MIN_ALT_READS, parser.DEFAULT_MIN_VAF))
    calls = pd.DataFrame(rows)
    calls["patient"] = calls["patient"].astype(str)
    calls["cohort"] = calls["patient"].map(lambda p: "PDAC" if p in PRIMARY else "periampullary")
    return calls


def matrix(calls: pd.DataFrame) -> pd.DataFrame:
    return calls.assign(sample=calls["patient"].map(lambda p: f"{p.zfill(2)}T")).groupby(["SYMBOL", "sample"]).size().unstack(fill_value=0).clip(upper=1)


def write_cohort(calls: pd.DataFrame, cohort: str, patients: list[str]) -> dict[str, pd.DataFrame]:
    out = OUT / ("primary_pdac" if cohort == "PDAC" else "periampullary")
    out.mkdir(parents=True, exist_ok=True)
    sub = calls[calls["patient"].isin(patients)].copy()
    mat = matrix(sub).reindex(columns=[f"{p.zfill(2)}T" for p in patients], fill_value=0)
    burden = sub.groupby("patient").size().reindex(patients, fill_value=0).rename("somatic_variant_count").reset_index()
    freq = mat.sum(axis=1).rename("n_patients").to_frame()
    freq["frequency"] = freq["n_patients"] / len(patients)
    freq = freq.sort_values(["n_patients"], ascending=False).reset_index(names="SYMBOL")
    classes = sub.groupby(["IMPACT", "Consequence"], dropna=False).size().reset_index(name="n_variants").sort_values("n_variants", ascending=False)
    drivers = freq[freq["SYMBOL"].isin(DRIVER_GENES)].copy()
    pathway_rows = []
    for pathway, genes in PATHWAYS.items():
        present = sorted(set(freq.loc[freq["SYMBOL"].isin(genes), "SYMBOL"]))
        pathway_rows.append({"pathway": pathway, "n_genes_observed": len(present), "genes_observed": ";".join(present), "n_patients_affected": int(mat.loc[mat.index.intersection(genes)].sum(axis=0).gt(0).sum())})
    pathway = pd.DataFrame(pathway_rows).sort_values("n_patients_affected", ascending=False)
    sub.to_csv(out / "wes_somatic_variants.tsv", sep="\t", index=False)
    mat.to_csv(out / "wes_mutation_matrix.csv")
    burden.to_csv(out / "wes_mutation_burden.csv", index=False)
    freq.to_csv(out / "mutation_frequency.csv", index=False)
    classes.to_csv(out / "variant_classification_summary.csv", index=False)
    drivers.to_csv(out / "driver_gene_summary.csv", index=False)
    pathway.to_csv(out / "pathway_mutation_summary.csv", index=False)
    return {"calls": sub, "matrix": mat, "burden": burden, "frequency": freq, "classes": classes, "drivers": drivers, "pathway": pathway}


def save_figures(primary: dict[str, pd.DataFrame], peri: dict[str, pd.DataFrame]) -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    for label, data in [("pdac11", primary), ("periampullary3", peri)]:
        burden = data["burden"]
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.bar(burden["patient"], burden["somatic_variant_count"], color="#264653", edgecolor="#1d3557")
        ax.set(xlabel="Patient", ylabel="Filtered somatic variants", title=f"WES somatic mutation burden: {label}")
        fig.tight_layout(); fig.savefig(FIG / f"{label}_mutation_burden.png", dpi=300); plt.close(fig)
        freq = data["frequency"].head(20).sort_values("frequency")
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.barh(freq["SYMBOL"], freq["frequency"] * 100, color="#2a9d8f")
        ax.set(xlabel="Patients with mutation (%)", ylabel="Gene", title=f"Top recurrent mutated genes: {label}")
        fig.tight_layout(); fig.savefig(FIG / f"{label}_mutation_frequency_top20.png", dpi=300); plt.close(fig)
        fig, ax = plt.subplots(figsize=(10, 7))
        top = data["frequency"].head(25)["SYMBOL"]
        heat = data["matrix"].reindex(top, fill_value=0)
        sns.heatmap(heat, cmap=["#f2f2f2", "#6a4c93"], cbar=False, linewidths=.3, linecolor="white", ax=ax)
        ax.set(xlabel="Tumor sample", ylabel="Gene", title=f"Recurrent mutation profile: {label}")
        fig.tight_layout(); fig.savefig(FIG / f"{label}_oncoprint_top25.png", dpi=300); plt.close(fig)
    comparison = pd.DataFrame({"cohort": ["PDAC", "periampullary"], "n_patients": [len(PRIMARY), len(PERI)], "n_variants": [len(primary["calls"]), len(peri["calls"])], "n_genes": [primary["matrix"].shape[0], peri["matrix"].shape[0]]})
    comparison.to_csv(TABLE / "cohort_summary.csv", index=False)
    top_genes = sorted(set(primary["matrix"].index) | set(peri["matrix"].index))
    combined = pd.concat([primary["matrix"].reindex(top_genes, fill_value=0), peri["matrix"].reindex(top_genes, fill_value=0)], axis=1).fillna(0)
    if combined.shape[1] > 1:
        linkage_matrix = linkage(combined.T, method="average", metric="jaccard")
        fig, ax = plt.subplots(figsize=(10, 5)); dendrogram(linkage_matrix, labels=combined.columns, ax=ax, leaf_rotation=45)
        ax.set_title("Hierarchical clustering of mutation profiles"); ax.set_ylabel("Average Jaccard distance")
        fig.tight_layout(); fig.savefig(FIG / "pdac11_periampullary_mutation_clustering.png", dpi=300); plt.close(fig)
    pd.DataFrame({"gene": top_genes, "PDAC_n": primary["matrix"].reindex(top_genes, fill_value=0).sum(axis=1), "periampullary_n": peri["matrix"].reindex(top_genes, fill_value=0).sum(axis=1)}).sort_values(["PDAC_n", "periampullary_n"], ascending=False).to_csv(TABLE / "pdac_vs_periampullary_gene_comparison.csv", index=False)


def write_reports(calls: pd.DataFrame, primary: dict[str, pd.DataFrame], peri: dict[str, pd.DataFrame]) -> None:
    shared = sorted(set(primary["matrix"].index) & set(peri["matrix"].index))
    pdac_only = sorted(set(primary["matrix"].index) - set(peri["matrix"].index))
    peri_only = sorted(set(peri["matrix"].index) - set(primary["matrix"].index))
    kras = calls[calls["SYMBOL"].eq("KRAS")].sort_values(["cohort", "patient"])
    kras.to_csv(TABLE / "kras_variants_all_cohorts.tsv", sep="\t", index=False)
    lines = [
        "REVISED WES ANALYSIS REPORT", "===========================", "",
        "Primary cohort: patients specified in private configuration (n=11).",
        "Exploratory periampullary cohort: patients specified in private configuration (n=3).", "",
        "Inputs", "------",
        "- Per-patient DeepSomatic 1.10.0 FFPE_WES calls using the vendor Exome 2.5 BED, annotated with VEP 101.",
        "- RNA-seq and proteomics inputs are unchanged and are reused only for the primary 11-patient integration.",
        "- MSI is generated separately by the paired MSI workflow.", "",
        "Regenerated analyses", "---------------------",
        "- Somatic variant tables, gene-by-patient matrices, burden, gene frequency, variant classification, driver-gene, and curated PDAC pathway summaries.",
        "- Mutation burden, recurrent-gene, oncoprint, and combined cohort clustering figures.",
        "- KRAS table and PDAC-versus-periampullary gene comparison.", "",
        "Cohort differences", "-------------------",
        f"- PDAC: {len(primary['calls'])} filtered variants across {primary['matrix'].shape[0]} genes.",
        f"- Periampullary: {len(peri['calls'])} filtered variants across {peri['matrix'].shape[0]} genes.",
        f"- Shared mutated genes: {len(shared)}; PDAC-only genes: {len(pdac_only)}; periampullary-only genes: {len(peri_only)}.",
        f"- KRAS calls: {len(kras)} total; PDAC={int((kras['cohort'] == 'PDAC').sum())}; periampullary={int((kras['cohort'] == 'periampullary').sum())}.",
        "- Clustering is exploratory and descriptive because the periampullary group has only three cases.", "",
        "Downstream integration note", "-----------------------------",
        "The primary PDAC WES matrix contains only patients specified in private configuration. The three periampullary patients are summarized separately and are not included in primary integration.",
        "This is the WES-derived analysis stage; it does not by itself refit MOFA or regenerate RNA/proteomics-only figures.",
    ]
    (RUN / "WES_2_5_ANALYSIS_REPORT.txt").write_text("\n".join(lines) + "\n")


def main() -> None:
    RUN.mkdir(parents=True, exist_ok=True); OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True); TABLE.mkdir(parents=True, exist_ok=True)
    calls = load_calls()
    calls.to_csv(TABLE / "all_14_filtered_calls_used_for_cohort_split.tsv", sep="\t", index=False)
    primary = write_cohort(calls, "PDAC", PRIMARY)
    peri = write_cohort(calls, "periampullary", PERI)
    save_figures(primary, peri)
    write_reports(calls, primary, peri)
    print(f"Primary variants: {len(primary['calls'])}; genes: {primary['matrix'].shape[0]}")
    print(f"Periampullary variants: {len(peri['calls'])}; genes: {peri['matrix'].shape[0]}")
    print(f"Outputs: {RUN}")


if __name__ == "__main__":
    main()
