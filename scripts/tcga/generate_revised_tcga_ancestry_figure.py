#!/usr/bin/env python3
# Current collaborator code export. Public-cohort preparation and descriptive comparisons; not ancestry-controlled validation.
# Source: Original Source and Code/scripts/generate_revised_tcga_ancestry_figure.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Rebuild TCGA PAAD ancestry/modality availability figures from source metadata."""
from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = project_root()
RUN = ROOT / "results" / "mofa_revised_11_20260916"
FIG = RUN / "figures"
TABLE = RUN / "tables"
TCGA = ROOT / "external" / "tcga_paad"


def read_sample_ids(path: Path) -> set[str]:
    return set(pd.read_csv(path, header=None, sep="\t", dtype=str).iloc[:, 0].dropna().astype(str))


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TABLE.mkdir(parents=True, exist_ok=True)
    clinical = pd.read_csv(TCGA / "clinical_tmp/paad_tcga_pan_can_atlas_2018/data_clinical_patient.txt", sep="\t", comment="#", dtype=str)
    clinical = clinical[clinical["PATIENT_ID"].notna()].copy()
    ancestry = clinical[clinical["GENETIC_ANCESTRY_LABEL"].isin(["AFR", "AFR_ADMIX"])].copy()
    ancestry_ids = set(ancestry["PATIENT_ID"].astype(str))
    rna = pd.read_csv(TCGA / "processed/tcga_rna_matrix.csv", index_col=0, nrows=1).columns.astype(str)
    rppa = pd.read_csv(TCGA / "processed/tcga_rppa_matrix.csv", index_col=0, nrows=1).columns.astype(str)
    mc3 = pd.read_csv(TCGA / "processed/tcga_mc3_gene_binary.csv", index_col=0, nrows=1).columns.astype(str)
    modality = {"RNA": set(rna), "RPPA": set(rppa), "MC3": set(mc3)}
    # TCGA matrices use sample IDs; reduce to patient IDs for ancestry matching.
    modality_patients = {name: {sample[:12] for sample in samples} for name, samples in modality.items()}
    full_counts = {name: len(values) for name, values in modality_patients.items()}
    ancestry_counts = {name: len(values & ancestry_ids) for name, values in modality_patients.items()}
    black_overlap = set.intersection(*(values & ancestry_ids for values in modality_patients.values()))
    full_overlap = set.intersection(*modality_patients.values())
    african_n = 11
    summary = pd.DataFrame({"quantity": ["RNA", "RPPA", "MC3", "all_three_modalities", "ancestry_labels_AFR_or_AFR_ADMIX"], "full_tcga": [full_counts["RNA"], full_counts["RPPA"], full_counts["MC3"], len(full_overlap), "NA"], "african_ancestry_subset": [ancestry_counts["RNA"], ancestry_counts["RPPA"], ancestry_counts["MC3"], len(black_overlap), len(ancestry_ids)]})
    summary.to_csv(TABLE / "revised_tcga_ancestry_modality_summary.csv", index=False)
    pd.DataFrame({"patient_id": sorted(black_overlap), "ancestry_label": "AFR_or_AFR_ADMIX", "has_rna": True, "has_rppa": True, "has_mc3": True}).to_csv(TABLE / "revised_tcga_ancestry_all_three_patients.csv", index=False)
    labels = ["RNA", "RPPA", "MC3"]
    full = [full_counts[x] for x in labels]
    subset = [ancestry_counts[x] for x in labels]
    x = range(len(labels)); width = .35
    fig, ax = plt.subplots(figsize=(9, 5.5)); ax.bar([i-width/2 for i in x], full, width, label="TCGA full", color="#457b9d"); ax.bar([i+width/2 for i in x], subset, width, label="AFR/AFR_ADMIX subset", color="#e76f51")
    ax.set_xticks(list(x)); ax.set_xticklabels(labels); ax.set_ylabel("Number of patients"); ax.set_title("TCGA PAAD modality availability by genetic ancestry")
    for i,(f,s) in enumerate(zip(full,subset)): ax.text(i-width/2,f+1,str(f),ha="center",fontsize=9); ax.text(i+width/2,s+1,str(s),ha="center",fontsize=9)
    ax.legend(frameon=False); ax.grid(axis="y",linestyle=":",alpha=.35); fig.tight_layout(); fig.savefig(FIG/"revised_tcga_ancestry_counts_by_modality.png",dpi=300); plt.close(fig)
    labels2=["Revised African PDAC", "TCGA AFR/AFR_ADMIX RNA", "TCGA AFR/AFR_ADMIX RPPA", "TCGA AFR/AFR_ADMIX MC3", "TCGA AFR/AFR_ADMIX all 3"]
    values=[african_n,ancestry_counts["RNA"],ancestry_counts["RPPA"],ancestry_counts["MC3"],len(black_overlap)]
    fig,ax=plt.subplots(figsize=(10,5.5)); ax.barh(labels2,values,color=["#2a9d8f","#e76f51","#e76f51","#e76f51","#8c564b"]); ax.set_xlabel("Number of patients"); ax.set_title("Effective sample sizes for ancestry-restricted comparisons")
    for i,v in enumerate(values): ax.text(v+.3,i,str(v),va="center",fontsize=9)
    ax.grid(axis="x",linestyle=":",alpha=.35); fig.tight_layout(); fig.savefig(FIG/"revised_tcga_ancestry_overlap_limit.png",dpi=300); plt.close(fig)
    print('Ancestry labels:',len(ancestry_ids),'all-three AFR/AFR_ADMIX overlap:',len(black_overlap))


if __name__ == "__main__":
    main()
