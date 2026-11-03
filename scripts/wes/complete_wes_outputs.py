#!/usr/bin/env python3
# PDAC2026 analysis: VEP-annotated DeepSomatic calls, retained-variant summaries and paired MSI launcher.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Summarise filtered WES calls, gene-level exclusions and paired MSI scores."""
from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

import re
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = project_root()
RUN = ROOT / "results" / "wes_revised_11_20260916"
FIG = RUN / "figures" / "extended"
TABLE = RUN / "tables"
MSI_DIR = RUN / "inputs" / "msi"
COHORTS = {"primary_pdac": configured_cohort("primary_pdac"), "periampullary": configured_cohort("periampullary")}
VALIDATED_MSI_FALLBACKS = {}
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
CODING_TERMS = {"missense_variant", "synonymous_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant", "coding_sequence_variant"}
REFSEQ_TO_CHROM = {f"NC_{i:06d}.{v}": f"chr{i}" for i, v in []}
REFSEQ_TO_CHROM.update({
    "NC_000001.11":"chr1", "NC_000002.12":"chr2", "NC_000003.12":"chr3", "NC_000004.12":"chr4", "NC_000005.10":"chr5", "NC_000006.12":"chr6", "NC_000007.14":"chr7", "NC_000008.11":"chr8", "NC_000009.12":"chr9", "NC_000010.11":"chr10", "NC_000011.10":"chr11", "NC_000012.12":"chr12", "NC_000013.11":"chr13", "NC_000014.9":"chr14", "NC_000015.10":"chr15", "NC_000016.10":"chr16", "NC_000017.11":"chr17", "NC_000018.10":"chr18", "NC_000019.10":"chr19", "NC_000020.11":"chr20", "NC_000021.9":"chr21", "NC_000022.11":"chr22", "NC_000023.11":"chrX", "NC_000024.10":"chrY", "NC_012920.1":"chrMT"
})
BLACKLIST_PATTERNS = {
    "Mucin (MUC*)": re.compile(r"^MUC\d"), "Zinc-finger paralog (ZNF*)": re.compile(r"^ZNF\d"), "Olfactory receptor (OR*)": re.compile(r"^OR\d"), "MAGE family": re.compile(r"^MAGE[A-Z]"), "NBPF family": re.compile(r"^NBPF\d"), "Uncurated/predicted locus (LOC*)": re.compile(r"^LOC\d")
}
BLACKLIST_GENES = {"TTN":"large structural gene", "OBSCN":"large structural gene", "AHNAK2":"large structural gene", "RP1L1":"large photoreceptor gene", "HRNR":"large repetitive gene", "HERC1":"large gene", "CUBN":"large gene", "PKD1L1":"large complex gene", "UBTFL1":"paralog-adjacent locus", "FCGBP":"large repetitive gene"}


def blacklist_reason(symbol: str) -> str | None:
    for label, pattern in BLACKLIST_PATTERNS.items():
        if pattern.match(str(symbol)):
            return label
    return BLACKLIST_GENES.get(str(symbol))


def variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1: return "SNP"
    if len(ref) < len(alt): return "INS"
    if len(ref) > len(alt): return "DEL"
    return "MNP"


def load_cohort(cohort: str) -> pd.DataFrame:
    return pd.read_csv(RUN / "outputs" / cohort / "wes_somatic_variants.tsv", sep="\t")


def savefig(name: str, fig: plt.Figure) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def complete_tables(cohort: str, df: pd.DataFrame) -> pd.DataFrame:
    out = RUN / "outputs" / cohort
    x = df.copy()
    x["blacklist_reason"] = x["SYMBOL"].map(blacklist_reason)
    removed = x[x["blacklist_reason"].notna()].copy()
    kept = x[x["blacklist_reason"].isna()].drop(columns="blacklist_reason")
    removed.groupby(["SYMBOL", "blacklist_reason"]).size().reset_index(name="n_variants_removed").sort_values("n_variants_removed", ascending=False).to_csv(out / "wes_blacklist_removed_genes.tsv", sep="\t", index=False)
    kept.to_csv(out / "wes_somatic_variants.blacklist_filtered.tsv", sep="\t", index=False)
    matrix = kept.assign(sample=kept["patient"].map(lambda p: f"{str(p).zfill(2)}T")).groupby(["SYMBOL", "sample"]).size().unstack(fill_value=0).clip(upper=1)
    matrix = matrix.reindex(columns=[f"{str(p).zfill(2)}T" for p in COHORTS[cohort]], fill_value=0)
    matrix.to_csv(out / "wes_mutation_matrix.blacklist_filtered.csv")
    genes = kept.groupby("SYMBOL").agg(patients_affected=("patient", "nunique"), variant_count=("SYMBOL", "size"), chroms=("CHROM", lambda s: ",".join(sorted(set(s)))), patient_list=("patient", lambda s: ",".join(sorted(set(map(str, s)))))).reset_index().sort_values(["patients_affected", "variant_count"], ascending=False)
    genes.to_csv(out / "gene_summary_table.csv", index=False)
    return kept


def extended_figures(cohort: str, df: pd.DataFrame) -> None:
    d = df.copy(); d["TYPE"] = [variant_type(r, a) for r, a in zip(d.REF, d.ALT)]
    d["chrom_name"] = d.CHROM.map(lambda x: REFSEQ_TO_CHROM.get(x, x))
    prefix = "pdac11" if cohort == "primary_pdac" else "periampullary3"
    counts = d.TYPE.value_counts(); fig, ax = plt.subplots(figsize=(5,4)); sns.barplot(x=counts.index, y=counts.values, ax=ax, hue=counts.index, palette="viridis", legend=False); ax.set(title=f"Variant types: {prefix}", xlabel="", ylabel="Count"); savefig(f"{prefix}_variant_types.png", fig)
    snps = d[d.TYPE == "SNP"]; ts = sum((r, a) in TRANSITIONS for r, a in zip(snps.REF, snps.ALT)); tv = len(snps) - ts; fig, ax = plt.subplots(figsize=(5,4)); sns.barplot(x=["Transitions", "Transversions"], y=[ts, tv], ax=ax, hue=["Transitions", "Transversions"], palette="Set2", legend=False); ax.set(title=f"Ts/Tv ratio: {prefix} ({ts/tv:.2f})", ylabel="Count"); savefig(f"{prefix}_ts_tv.png", fig)
    subs = Counter(f"{r}>{a}" for r, a in zip(snps.REF, snps.ALT)); fig, ax = plt.subplots(figsize=(7,4)); sns.barplot(x=list(subs), y=list(subs.values()), ax=ax, hue=list(subs), palette="magma", legend=False); ax.set(title=f"Substitution spectrum: {prefix}", ylabel="Count"); ax.tick_params(axis="x", rotation=45); savefig(f"{prefix}_substitutions.png", fig)
    impacts = d.IMPACT.value_counts().reindex(["HIGH", "MODERATE", "LOW", "MODIFIER"]).dropna(); fig, ax = plt.subplots(figsize=(5,4)); sns.barplot(x=impacts.index, y=impacts.values, ax=ax, hue=impacts.index, palette={"HIGH":"#d62728","MODERATE":"#e9c46a","LOW":"#2a9d8f","MODIFIER":"#264653"}, legend=False); ax.set(title=f"VEP impact distribution: {prefix}", ylabel="Count"); savefig(f"{prefix}_impact_distribution.png", fig)
    coding = d.Consequence.fillna("").map(lambda x: bool(set(str(x).split("&")) & CODING_TERMS)); fig, ax = plt.subplots(figsize=(4,4)); coding.value_counts().rename({True:"Coding",False:"Non-coding"}).plot.pie(autopct="%1.1f%%", ax=ax, colors=["#2a9d8f","#d9d9d9"]); ax.set(title=f"Coding summary: {prefix}", ylabel=""); savefig(f"{prefix}_coding_summary.png", fig)
    genes = d.groupby("SYMBOL").agg(patients=("patient","nunique"), variants=("SYMBOL","size")).sort_values(["patients","variants"], ascending=False).head(20).sort_values("patients"); fig, ax = plt.subplots(figsize=(6,7)); ax.barh(genes.index, genes.patients, color="#6a4c93"); ax.set(title=f"Top mutated genes: {prefix}", xlabel="Patients affected"); savefig(f"{prefix}_top_genes.png", fig)
    chrom = d.chrom_name.value_counts(); chrom = chrom.reindex(sorted(chrom.index, key=lambda x: int(x[3:]) if x.startswith("chr") and x[3:].isdigit() else 99)); fig, ax = plt.subplots(figsize=(9,4)); ax.bar(chrom.index, chrom.values, color=sns.color_palette("crest", len(chrom))); ax.set(title=f"Variants per chromosome: {prefix}", ylabel="Count"); ax.tick_params(axis="x", rotation=45); savefig(f"{prefix}_genomic_position.png", fig)
    indels = d[d.TYPE.isin(["INS", "DEL"])].copy(); indels["length"] = indels.ALT.str.len() - indels.REF.str.len(); fig, ax = plt.subplots(figsize=(6,4)); ax.hist(indels.length, bins=range(-30, 31), color="#e76f51"); ax.set(title=f"Indel length distribution: {prefix}", xlabel="Length (ALT - REF, bp)", ylabel="Count"); savefig(f"{prefix}_indel_length.png", fig)
    consequences = Counter(term for value in d.Consequence.dropna() for term in str(value).split("&")); top_consequences = dict(sorted(consequences.items(), key=lambda item: item[1], reverse=True)[:15]); fig, ax = plt.subplots(figsize=(8,5)); ax.barh(list(top_consequences)[::-1], list(top_consequences.values())[::-1], color="#457b9d"); ax.set(title=f"Top consequence types: {prefix}", xlabel="Count"); savefig(f"{prefix}_consequence_types.png", fig)
    for chromosome in ["chr2", "chr6", "chr7", "chr11"]:
        sub = d[d.chrom_name.eq(chromosome)].sort_values("POS")
        if sub.empty:
            continue
        fig, ax = plt.subplots(figsize=(8,3)); ax.scatter(sub.POS, range(len(sub)), s=5, alpha=.55, color="#2a9d8f"); ax.set(title=f"Manhattan-style plot: {prefix}, {chromosome}", xlabel="Genomic position", ylabel="Variant index"); savefig(f"{prefix}_manhattan_{chromosome}.png", fig)
    windows = d.assign(window=(d.POS // 1_000_000) * 1_000_000).groupby(["chrom_name","window"]).size().reset_index(name="n").sort_values("n", ascending=False).head(20); windows["label"] = windows.chrom_name + ":" + (windows.window // 1_000_000).astype(str) + "Mb"; fig, ax = plt.subplots(figsize=(8,5)); ax.barh(windows.label.iloc[::-1], windows.n.iloc[::-1], color="#e76f51"); ax.set(title=f"Top 1-Mb variant windows: {prefix}", xlabel="Variant count"); savefig(f"{prefix}_hotspots.png", fig)
    d.groupby("patient").size().rename("variant_count").to_csv(TABLE / f"{prefix}_mutations_per_patient.csv")
    d.groupby("patient").SYMBOL.nunique().rename("unique_genes").to_csv(TABLE / f"{prefix}_genes_per_patient.csv")


def msi_outputs() -> None:
    rows = []
    for cohort, patients in COHORTS.items():
        for patient in patients:
            path = MSI_DIR / f"patient_{patient}"
            if not path.exists() and patient not in VALIDATED_MSI_FALLBACKS:
                continue
            try:
                if patient in VALIDATED_MSI_FALLBACKS and not path.exists():
                    total, unstable, score = VALIDATED_MSI_FALLBACKS[patient]
                else:
                    x = pd.read_csv(path, sep=r"\s+", engine="python")
                    row = x.iloc[0].to_dict(); score = float(row.get("%", np.nan)); total = int(row.get("Total_Number_of_Sites", np.nan)); unstable = int(row.get("Number_of_Unstable_Sites", np.nan))
                category = "MSI-High" if score >= 20 else "MSI-Intermediate" if score >= 5 else "MSS"
                rows.append({"patient": patient, "cohort": cohort, "informative_sites": total, "unstable_sites": unstable, "msi_score_pct": score, "category": category})
            except Exception:
                pass
    msi = pd.DataFrame(rows).sort_values(["cohort", "msi_score_pct"], ascending=[True, False]); msi.to_csv(TABLE / "msisensor_pro_revised_cohorts.csv", index=False)
    if not msi.empty:
        fig, ax = plt.subplots(figsize=(8,5)); sns.stripplot(data=msi, x="cohort", y="msi_score_pct", hue="category", dodge=False, jitter=.15, palette={"MSS":"#2a9d8f","MSI-Intermediate":"#e9c46a","MSI-High":"#d62728"}, size=9, ax=ax); ax.axhline(5, ls="--", color="grey"); ax.axhline(20, ls="--", color="grey"); ax.set(xlabel="Cohort", ylabel="MSI score (% unstable microsatellite sites)", title="MSIsensor-pro MSI scores in revised WES cohorts"); fig.tight_layout(); fig.savefig(FIG / "revised_cohort_msi_scores.png", dpi=300); plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True); TABLE.mkdir(parents=True, exist_ok=True)
    for cohort in COHORTS:
        raw = pd.read_csv(RUN / "outputs" / cohort / "wes_somatic_variants.tsv", sep="\t")
        kept = complete_tables(cohort, raw)
        extended_figures(cohort, kept)
    msi_outputs()
    print("Completed blacklist-filtered tables, extended WES figures, gene summaries, and MSI outputs")


if __name__ == "__main__":
    main()
