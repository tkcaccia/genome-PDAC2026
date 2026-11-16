#!/usr/bin/env python3
# MOFA fitting, molecular integration and publication panels.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""
from __future__ import annotations

from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path

from pathlib import Path
import re
import importlib.util

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
from scipy.stats import spearmanr, fisher_exact

ROOT = private_path("handoff_root")
WORKSPACE_ROOT = project_root()
RUN = ROOT / "mofa" / "pdac_11"
INTEGRATION = RUN / "integration"
QC = RUN / "qc"
MOFA = RUN / "mofa"
FIG = RUN / "figures"
TABLE = RUN / "tables"
PATIENTS = [f"{p}T" for p in configured_cohort("primary_pdac")]


def normalize_sample(x: str) -> str:
    text = str(x)
    match = re.match(r"^0*(\d+)T$", text)
    return f"{int(match.group(1))}T" if match else text


def load_wes() -> pd.DataFrame:
    path = ROOT / "wes/analysis/outputs/primary_pdac/wes_mutation_matrix.blacklist_filtered.csv"
    df = pd.read_csv(path, index_col=0).fillna(0).astype(float)
    df.columns = [normalize_sample(c) for c in df.columns]
    return df.loc[:, PATIENTS]


def load_rna() -> pd.DataFrame:
    path = WORKSPACE_ROOT / "results/mofa_revised_11_20260916/integration/rnaseq_delta.csv"
    df = pd.read_csv(path, index_col=0).groupby(level=0).mean()
    df.columns = [f"{normalize_sample(c)}" if str(c).endswith("T") else f"{normalize_sample(c)}T" for c in df.columns]
    return df.loc[:, PATIENTS]


def load_proteomics() -> pd.DataFrame:
    path = WORKSPACE_ROOT / "results/mofa_revised_11_20260916/integration/proteomics_log2.csv"
    df = pd.read_csv(path, index_col=0)
    df.columns = [normalize_sample(c) for c in df.columns]
    return df.loc[:, PATIENTS]


def write_integration(wes: pd.DataFrame, rna: pd.DataFrame, prot: pd.DataFrame) -> None:
    INTEGRATION.mkdir(parents=True, exist_ok=True)
    wes.to_csv(INTEGRATION / "wes_binary.csv")
    rna.to_csv(INTEGRATION / "rnaseq_delta.csv")
    prot.to_csv(INTEGRATION / "proteomics_log2.csv")
    pd.Series(PATIENTS, name="patient").to_csv(INTEGRATION / "samples.tsv", sep="\t", index=False)


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, 1), axis=0).fillna(0)


def qc_views(wes: pd.DataFrame, rna: pd.DataFrame, prot: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    QC.mkdir(parents=True, exist_ok=True)
    wes_keep = wes.sum(axis=1) / wes.shape[1] >= 0.10
    rna_keep = rna.var(axis=1) > 0
    prot_keep = prot.notna().mean(axis=1) >= 0.20
    wes_qc = wes.loc[wes_keep]
    rna_qc = zscore_rows(rna.loc[rna_keep])
    prot_qc = zscore_rows(prot.loc[prot_keep])
    wes_qc.to_csv(QC / "wes_binary_qc.csv")
    rna_qc.to_csv(QC / "rnaseq_delta_zscore.csv")
    prot_qc.to_csv(QC / "proteomics_log2_zscore.csv")
    pd.DataFrame({"view": ["wes", "rnaseq", "proteomics"], "input_features": [len(wes), len(rna), len(prot)], "retained_features": [len(wes_qc), len(rna_qc), len(prot_qc)], "n_samples": [len(PATIENTS)] * 3}).to_csv(QC / "qc_summary.csv", index=False)
    pd.DataFrame({"wes_prevalence": wes.sum(axis=1) / wes.shape[1]}).to_csv(QC / "wes_binary_qc_feature_prevalence.csv")
    return wes_qc, rna_qc, prot_qc


def top_variable(df: pd.DataFrame, n: int) -> pd.DataFrame:
    return df.loc[df.var(axis=1).sort_values(ascending=False).head(min(n, len(df))).index]


def run_mofa(wes: pd.DataFrame, rna: pd.DataFrame, prot: pd.DataFrame) -> Path:
    from mofapy2.run.entry_point import entry_point
    MOFA.mkdir(parents=True, exist_ok=True)
    wes = top_variable(wes, 600); rna = top_variable(rna, 2000); prot = top_variable(prot, 1200)
    pd.DataFrame({"view": ["wes", "rnaseq", "proteomics"], "features_used": [len(wes), len(rna), len(prot)], "samples": [len(PATIENTS)] * 3}).to_csv(TABLE / "mofa_feature_selection.csv", index=False)
    data = [[wes.T.to_numpy()], [rna.T.to_numpy()], [prot.T.to_numpy()]]
    ent = entry_point()
    outfile = MOFA / "mofa_model.hdf5"
    ent.set_data_matrix(data=data, likelihoods=["bernoulli", "gaussian", "gaussian"], views_names=["wes", "rnaseq", "proteomics"], groups_names=["all"], samples_names=[PATIENTS], features_names=[wes.index.tolist(), rna.index.tolist(), prot.index.tolist()])
    ent.set_data_options(scale_views=False, scale_groups=False, center_groups=True, use_float32=False)
    ent.set_model_options(factors=5, spikeslab_factors=False, spikeslab_weights=True, ard_factors=False, ard_weights=True)
    ent.set_train_options(iter=1500, startELBO=1, freqELBO=5, startSparsity=50, tolerance=1e-5, convergence_mode="fast", startDrop=1, freqDrop=5, verbose=False, quiet=True, seed=42, outfile=str(outfile))
    ent.build(); ent.run(); ent.save(outfile=str(outfile), save_data=True, save_parameters=True)
    return outfile


def decode(values):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def export_mofa(model: Path) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    out = MOFA / "exports"; out.mkdir(parents=True, exist_ok=True)
    with h5py.File(model, "r") as f:
        samples = decode(f["samples/all"][()]); z = np.asarray(f["expectations/Z/all"][()]).T
        factors = pd.DataFrame(z, index=samples, columns=[f"factor_{i+1}" for i in range(z.shape[1])]); factors.to_csv(out / "mofa_factors.csv")
        loadings = {}
        for view in ["wes", "rnaseq", "proteomics"]:
            features = decode(f[f"features/{view}"][()]); weights = np.asarray(f[f"expectations/W/{view}"][()]).T
            loadings[view] = pd.DataFrame(weights, index=features, columns=factors.columns)
            loadings[view].to_csv(out / f"mofa_loadings_{view}.csv")
    return factors, loadings


def cluster_and_figures(factors: pd.DataFrame) -> pd.DataFrame:
    out = MOFA / "clusters"; out.mkdir(parents=True, exist_ok=True)
    scores = {}
    for k in range(2, min(5, len(factors) - 1) + 1):
        labels = KMeans(n_clusters=k, random_state=42, n_init=20).fit_predict(factors)
        scores[k] = silhouette_score(factors, labels)
    best_k = max(scores, key=scores.get)
    labels = KMeans(n_clusters=best_k, random_state=42, n_init=20).fit_predict(factors)
    clusters = pd.DataFrame({"cluster": labels}, index=factors.index); clusters.to_csv(out / "mofa_patient_clusters.csv")
    pd.DataFrame({"k": list(scores), "silhouette": list(scores.values())}).to_csv(out / "cluster_k_selection.csv", index=False)
    fig, ax = plt.subplots(figsize=(7, 5)); ax.bar(list(scores), list(scores.values()), color="#2a9d8f"); ax.set(xlabel="Number of clusters", ylabel="Silhouette score", title="MOFA cluster-number selection (11 PDAC patients)"); fig.tight_layout(); fig.savefig(FIG / "mofa_cluster_k_selection.png", dpi=300); plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5)); x = factors.iloc[:, 0]; y = factors.iloc[:, 1]; sc = ax.scatter(x, y, c=labels, cmap="viridis", s=90, edgecolor="black");
    for i in factors.index: ax.annotate(i, (x.loc[i], y.loc[i]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set(xlabel="Factor 1", ylabel="Factor 2", title=f"MOFA factor clusters: k={best_k}, silhouette={scores[best_k]:.3f}"); fig.colorbar(sc, ax=ax, label="Cluster"); fig.tight_layout(); fig.savefig(FIG / "mofa_factor_clusters.png", dpi=300); plt.close(fig)
    return clusters


def interpretation_figures(factors: pd.DataFrame, loadings: dict[str, pd.DataFrame], clusters: pd.DataFrame) -> pd.DataFrame:
    out = MOFA / "interpretation"; out.mkdir(parents=True, exist_ok=True)
    f = factors.copy(); f["cluster"] = clusters["cluster"]; f = f.sort_values("cluster"); palette = {c: col for c, col in zip(sorted(f.cluster.unique()), ["#264653", "#6A4C93", "#2A9D8F", "#E9C46A"])}
    g = sns.clustermap(f.drop(columns="cluster"), cmap="viridis", row_cluster=False, col_cluster=False, row_colors=f.cluster.map(palette), figsize=(8, 6), center=0, dendrogram_ratio=(0.12, 0.02), colors_ratio=0.03, cbar_pos=(0.02, 0.78, 0.03, 0.15)); g.ax_heatmap.set_title("MOFA patient factor heatmap: 11 PDAC patients", loc="left", pad=14); g.fig.savefig(FIG / "mofa_factor_heatmap.png", dpi=300, bbox_inches="tight"); plt.close("all")
    d = factors.copy(); d["cluster"] = clusters["cluster"].astype(str); g = sns.pairplot(d, hue="cluster", corner=True, diag_kind="hist", plot_kws={"s": 55, "alpha": .9}); g.fig.suptitle("MOFA factor space by cluster: 11 PDAC patients", y=1.02); g.fig.savefig(FIG / "mofa_factor_pairplot.png", dpi=300); plt.close("all")
    rows = []
    for view, df in loadings.items():
        for factor in df.columns:
            for direction, sub in [("positive", df[factor].nlargest(20)), ("negative", df[factor].nsmallest(20))]:
                rows.extend({"view": view, "factor": factor, "direction": direction, "feature": feat, "loading": float(val), "abs_loading": abs(float(val))} for feat, val in sub.items())
    top = pd.DataFrame(rows); top.to_csv(out / "mofa_top_features_signed.csv", index=False); shortlist = top.sort_values("abs_loading", ascending=False).groupby("view", as_index=False).head(30); shortlist.to_csv(out / "mofa_preliminary_biomarker_shortlist.csv", index=False)
    sep = []
    for factor in factors.columns:
        grouped = factors.assign(cluster=clusters.cluster).groupby("cluster")[factor].agg(["mean", "std"])
        smd = np.nan if len(grouped) < 2 else (grouped["mean"].iloc[1] - grouped["mean"].iloc[0]) / np.sqrt((grouped["std"].fillna(0).iloc[0] ** 2 + grouped["std"].fillna(0).iloc[1] ** 2) / 2)
        sep.append({"factor": factor, "cluster_separation_smd": smd, "variance_across_patients": factors[factor].var()})
    pd.DataFrame(sep).to_csv(out / "mofa_factor_separation.csv", index=False)
    return shortlist




def robustness(factors: pd.DataFrame, clusters: pd.DataFrame) -> None:
    out = MOFA / "robustness"; out.mkdir(parents=True, exist_ok=True); base = clusters.cluster; scenarios = private_value("mofa_exclusion_scenarios"); rows = []
    for name, drops in scenarios.items():
        sub = factors.drop(index=drops, errors="ignore"); k = int(base.nunique()); labels = KMeans(n_clusters=k, random_state=42, n_init=20).fit_predict(sub); sil = silhouette_score(sub, labels) if len(sub) > k else np.nan; other = pd.Series(labels, index=sub.index); common = sorted(set(base.index) & set(other.index)); ari = adjusted_rand_score(base.loc[common], other.loc[common]) if len(common) >= 3 else np.nan; rows.append({"scenario": name, "n_patients": len(sub), "silhouette": sil, "ari_vs_full": ari})
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "mofa_cluster_sensitivity_summary.csv", index=False)
    fig, ax = plt.subplots(figsize=(8,4))
    ax.bar(summary.scenario, summary.ari_vs_full, color="#2a9d8f")
    ax.set_ylim(0,1.05)
    ax.set_ylabel("ARI vs full clustering")
    ax.set_title("MOFA clustering robustness: 11 PDAC patients")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(FIG / "mofa_cluster_robustness.png", dpi=300)
    plt.close(fig)




def main() -> None:
    for p in [INTEGRATION, QC, MOFA, FIG, TABLE]: p.mkdir(parents=True, exist_ok=True)
    wes, rna, prot = load_wes(), load_rna(), load_proteomics(); write_integration(wes, rna, prot); wes_qc, rna_qc, prot_qc = qc_views(wes, rna, prot); model = run_mofa(wes_qc, rna_qc, prot_qc); factors, loadings = export_mofa(model); clusters = cluster_and_figures(factors); shortlist = interpretation_figures(factors, loadings, clusters); pass; robustness(factors, clusters)
    (RUN / "MOFA_REVISED_11_REPORT.txt").write_text("11-patient PDAC MOFA integration\n\nSamples: " + ", ".join(PATIENTS) + "\nViews: blacklist-filtered WES, tumour-normal RNA log2 delta, log2 proteomics abundance.\nMOFA settings: 5 factors, 1500 iterations, seed 42, fast convergence.\nAll outputs are isolated under this directory.\nThe 11-patient model is exploratory and requires robustness-aware interpretation.\n")
    pass
    print("Completed isolated 11-patient MOFA workflow:", RUN)


if __name__ == "__main__":
    main()
