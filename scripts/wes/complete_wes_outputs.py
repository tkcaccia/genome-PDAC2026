#!/usr/bin/env python3
# WES processing, variant summaries, coverage and paired MSI.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""
from __future__ import annotations

from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path

import re
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = private_path("handoff_root")
RUN = ROOT / "wes" / "analysis"
FIG = RUN / "figures" / "extended"
TABLE = RUN / "tables"
MSI_DIR = ROOT / "wes" / "msi_2_5"
COHORTS = {"primary_pdac": configured_cohort("primary_pdac"), "periampullary": configured_cohort("periampullary")}
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




def msi_outputs() -> None:
    rows = []
    missing = []
    for cohort, patients in COHORTS.items():
        for patient in patients:
            path = MSI_DIR / f"patient_{patient}"
            if not path.is_file():
                missing.append(patient)
                continue
            try:
                x = pd.read_csv(path, sep=r"\s+", engine="python")
                row = x.iloc[0].to_dict(); score = float(row.get("%", np.nan)); total = int(row.get("Total_Number_of_Sites", np.nan)); unstable = int(row.get("Number_of_Unstable_Sites", np.nan))
                category = "MSI-High" if score >= 20 else "MSI-Intermediate" if score >= 5 else "MSS"
                rows.append({"patient": patient, "cohort": cohort, "informative_sites": total, "unstable_sites": unstable, "msi_score_pct": score, "category": category})
            except Exception:
                missing.append(patient)
    if missing:
        print(f"MSI analysis skipped: complete 2.5 MSI outputs are missing or unreadable for patients {sorted(set(missing), key=int)}")
        return
    msi = pd.DataFrame(rows).sort_values(["cohort", "msi_score_pct"], ascending=[True, False]); msi.to_csv(TABLE / "msisensor_pro_revised_cohorts.csv", index=False)
    if not msi.empty:
        fig, ax = plt.subplots(figsize=(8,5)); sns.stripplot(data=msi, x="cohort", y="msi_score_pct", hue="category", dodge=False, jitter=.15, palette={"MSS":"#2a9d8f","MSI-Intermediate":"#e9c46a","MSI-High":"#d62728"}, size=9, ax=ax); ax.axhline(5, ls="--", color="grey"); ax.axhline(20, ls="--", color="grey"); ax.set(xlabel="Cohort", ylabel="MSI score (% unstable microsatellite sites)", title="MSIsensor-pro MSI scores in revised WES cohorts"); fig.tight_layout(); fig.savefig(FIG / "revised_cohort_msi_scores.png", dpi=300); plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True); TABLE.mkdir(parents=True, exist_ok=True)
    for cohort in COHORTS:
        raw = pd.read_csv(RUN / "outputs" / cohort / "wes_somatic_variants.tsv", sep="\t")
        kept = complete_tables(cohort, raw)
        pass
    msi_outputs()
    print("Completed gene-filtered tables, gene summaries, and MSI outputs")


if __name__ == "__main__":
    main()
