#!/usr/bin/env python3
# Current collaborator code export. Three-view MOFA, exported factors/loadings, candidate ranking and sensitivity summaries.
# Source: Original Source and Code/scripts/run_mofa_revised_11.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Run the complete 11-patient WES/RNA/proteomics MOFA workflow in isolation."""
from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

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

ROOT = project_root()
RUN = ROOT / "results" / "mofa_revised_11_20260916"
INTEGRATION = RUN / "integration"
QC = RUN / "qc"
MOFA = RUN / "mofa"
FIG = RUN / "figures"
TABLE = RUN / "tables"
PATIENTS = [f"{patient}T" for patient in configured_cohort("primary_pdac")]


def normalize_sample(x: str) -> str:
    text = str(x)
    match = re.match(r"^0*(\d+)T$", text)
    return f"{int(match.group(1))}T" if match else text


def load_wes() -> pd.DataFrame:
    path = ROOT / "results/wes_revised_11_20260916/outputs/primary_pdac/wes_mutation_matrix.blacklist_filtered.csv"
    df = pd.read_csv(path, index_col=0).fillna(0).astype(float)
    df.columns = [normalize_sample(c) for c in df.columns]
    return df.loc[:, PATIENTS]


def load_rna() -> pd.DataFrame:
    path = ROOT / "results/rnaseq_revised_11_20260916/tables/tumour_normal_log2_delta_matrix_11_pdac.csv"
    df = pd.read_csv(path, index_col=0).groupby(level=0).mean()
    df.columns = [f"{normalize_sample(c)}" if str(c).endswith("T") else f"{normalize_sample(c)}T" for c in df.columns]
    return df.loc[:, PATIENTS]


def load_proteomics() -> pd.DataFrame:
    # The model uses tumour protein abundance, not a paired protein response.
    path = ROOT / "results/proteomics/protein_abundance_matrix.csv"
    df = pd.read_csv(path, index_col=0)
    cols = [c for c in df.columns if str(c).endswith("T")]
    df = df[cols].apply(pd.to_numeric, errors="coerce")
    df.columns = [normalize_sample(c) for c in df.columns]
    return np.log2(df.loc[:, PATIENTS].fillna(0) + 1)


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


def pathway_and_candidates(loadings: dict[str, pd.DataFrame], shortlist: pd.DataFrame) -> None:
    out = MOFA / "biology"; out.mkdir(parents=True, exist_ok=True)
    try:
        import gseapy as gp
        import mygene
        rna = loadings["rnaseq"]
        mg = mygene.MyGeneInfo(); ids = [str(x).split(".")[0] for x in rna.index]; mapping = {}
        try:
            res = mg.querymany(sorted(set(ids)), scopes="ensembl.gene", fields="symbol", species="human", as_dataframe=True)
            if isinstance(res, pd.DataFrame):
                for q in set(ids):
                    if q in res.index:
                        row = res.loc[q]
                        row = row.iloc[0] if isinstance(row, pd.DataFrame) else row
                        symbol = row.get("symbol", q)
                        mapping[q] = symbol if pd.notna(symbol) and str(symbol).strip() else q
        except Exception: pass
        rows = []
        for factor in rna.columns:
            genes = [mapping.get(str(x).split(".")[0], str(x)) for x in rna[factor].abs().nlargest(300).index]
            genes = list(dict.fromkeys(g for g in genes if not str(g).startswith("ENSG")))
            if len(genes) < 10: continue
            for gs in ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]:
                try:
                    enr = gp.enrichr(gene_list=genes, gene_sets=gs, organism="human", outdir=None, cutoff=1.0)
                    if enr is not None and enr.results is not None and len(enr.results):
                        x = enr.results.head(25).copy(); x.insert(0, "gene_set", gs); x.insert(0, "factor", factor); rows.append(x)
                except Exception: pass
        enrichment = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["factor", "gene_set", "Term", "Adjusted P-value", "Combined Score", "Genes"])
    except Exception:
        enrichment = pd.DataFrame(columns=["factor", "gene_set", "Term", "Adjusted P-value", "Combined Score", "Genes"])
    enrichment.to_csv(out / "mofa_pathway_enrichment.csv", index=False)
    scores = []
    for view, df in loadings.items():
        x = df.abs().max(axis=1).rename("score").to_frame(); x["feature"] = x.index; x["view"] = view; scores.append(x.reset_index(drop=True))
    rank = pd.concat(scores, ignore_index=True); rank["score_percentile"] = rank.groupby("view")["score"].rank(pct=True); rank = rank.sort_values("score_percentile", ascending=False).reset_index(drop=True); rank.insert(0, "rank", range(1, len(rank) + 1)); rank.to_csv(out / "mofa_cross_modal_target_ranking.csv", index=False); rank.head(200).to_csv(out / "mofa_cross_modal_target_top200.csv", index=False)
    if not enrichment.empty:
        top = enrichment.sort_values("Combined Score", ascending=False).groupby("factor").head(5); fig, ax = plt.subplots(figsize=(10, max(5, len(top) * .3))); ax.barh(range(len(top)), top["Combined Score"].fillna(0).iloc[::-1], color="#457b9d"); ax.set_yticks(range(len(top))); ax.set_yticklabels(top["Term"].iloc[::-1], fontsize=7); ax.set_title("MOFA RNA pathway enrichment: 11 PDAC patients"); fig.tight_layout(); fig.savefig(FIG / "mofa_pathway_enrichment.png", dpi=300); plt.close(fig)


def robustness(factors: pd.DataFrame, clusters: pd.DataFrame) -> None:
    # Recluster retained rows of fixed factor scores; this does not retrain MOFA.
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


def downstream_figures(factors: pd.DataFrame, loadings: dict[str, pd.DataFrame], shortlist: pd.DataFrame, rna: pd.DataFrame) -> None:
    ranked = []
    for view, df in loadings.items():
        score = df.abs().max(axis=1).rename("mofa_score").to_frame()
        score["feature"] = score.index.astype(str)
        score["view"] = view
        ranked.append(score.reset_index(drop=True))
    ranked = pd.concat(ranked, ignore_index=True).sort_values("mofa_score", ascending=False).reset_index(drop=True)
    ranked.insert(0, "rank", range(1, len(ranked) + 1))
    top = ranked.head(100)
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.barplot(data=top.sort_values("mofa_score"), x="mofa_score", y="feature", hue="view", dodge=False, ax=ax, palette={"wes":"#264653", "rnaseq":"#e76f51", "proteomics":"#2a9d8f"})
    ax.set_title("Top 100 cross-modal MOFA features: 11 PDAC patients")
    ax.set_xlabel("Maximum absolute MOFA loading")
    ax.set_ylabel("Feature")
    ax.legend(title="View", frameon=False)
    fig.tight_layout(); fig.savefig(FIG / "mofa_top100_cross_modal_features.png", dpi=300); plt.close(fig)

    long = []
    for view, df in loadings.items():
        x = df.abs().max(axis=1).rename("abs_loading").to_frame()
        x["view"] = view
        long.append(x.reset_index(names="feature"))
    load_long = pd.concat(long, ignore_index=True)
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.boxplot(data=load_long, x="view", y="abs_loading", hue="view", palette={"wes":"#264653", "rnaseq":"#e76f51", "proteomics":"#2a9d8f"}, legend=False, ax=ax)
    ax.set_title("MOFA loading distributions by view")
    ax.set_xlabel("")
    ax.set_ylabel("Maximum absolute loading")
    fig.tight_layout(); fig.savefig(FIG / "mofa_loading_distributions.png", dpi=300); plt.close(fig)

    rna_top = shortlist[shortlist.view.eq("rnaseq")].sort_values("abs_loading", ascending=False).head(30)
    present = [g for g in rna_top.feature if g in rna.index]
    if present:
        plot = pd.DataFrame({"gene": present, "mean_log2FC": rna.loc[present].mean(axis=1).values, "sd": rna.loc[present].std(axis=1).values}).sort_values("mean_log2FC")
        fig, ax = plt.subplots(figsize=(9, max(6, len(plot) * .24)))
        ax.barh(plot.gene, plot.mean_log2FC, xerr=plot.sd, color=np.where(plot.mean_log2FC >= 0, "#e76f51", "#457b9d"), edgecolor="black", capsize=2)
        ax.axvline(0, color="black", lw=.7)
        ax.set(title="RNA expression response for top MOFA RNA features", xlabel="Mean log2(Tumour / Normal)", ylabel="Ensembl gene")
        fig.tight_layout(); fig.savefig(FIG / "rna_expression_by_mofa_feature.png", dpi=300); plt.close(fig)

    mofa_load = loadings["rnaseq"].abs().max(axis=1)
    mean_fc = rna.mean(axis=1)
    common = sorted(set(mofa_load.index) & set(mean_fc.index))
    joined = pd.DataFrame({"abs_loading": mofa_load.loc[common], "mean_log2FC": mean_fc.loc[common]})
    rho, pval = spearmanr(joined.abs_loading, joined.mean_log2FC.abs())
    joined.to_csv(TABLE / "rna_mofa_loading_vs_log2fc.csv")
    fig, ax = plt.subplots(figsize=(6, 6)); ax.scatter(joined.abs_loading, joined.mean_log2FC, s=10, alpha=.35, color="#457b9d")
    ax.set(xlabel="Maximum absolute MOFA RNA loading", ylabel="Mean log2(Tumour / Normal)", title=f"MOFA RNA loading vs expression response\nSpearman rho={rho:.2f}, p={pval:.2g}")
    fig.tight_layout(); fig.savefig(FIG / "rna_mofa_loading_vs_log2fc.png", dpi=300); plt.close(fig)

    threshold = 1.0
    de = set(rna.index[rna.mean(axis=1).abs() >= threshold])
    mofa_top_genes = set(rna_top.feature)
    overlap = len(de & mofa_top_genes)
    overlap_table = pd.DataFrame({"category": ["overlap", "mofa_only", "de_only", "neither"], "count": [overlap, len(mofa_top_genes - de), len(de - mofa_top_genes), len(rna.index) - len(de | mofa_top_genes)]})
    overlap_table.to_csv(TABLE / "rna_mofa_overlap_table.csv", index=False)
    fig, ax = plt.subplots(figsize=(6, 4)); ax.bar(["Overlap", "MOFA only"], [overlap, len(mofa_top_genes - de)], color=["#2a9d8f", "#d9d9d9"])
    ax.set_ylabel("Number of top RNA MOFA features"); ax.set_title("MOFA RNA feature overlap with |mean log2FC| >= 1")
    fig.tight_layout(); fig.savefig(FIG / "rna_mofa_concordance_bar.png", dpi=300); plt.close(fig)


def main() -> None:
    for p in [INTEGRATION, QC, MOFA, FIG, TABLE]: p.mkdir(parents=True, exist_ok=True)
    wes, rna, prot = load_wes(), load_rna(), load_proteomics(); write_integration(wes, rna, prot); wes_qc, rna_qc, prot_qc = qc_views(wes, rna, prot); model = run_mofa(wes_qc, rna_qc, prot_qc); factors, loadings = export_mofa(model); clusters = cluster_and_figures(factors); shortlist = interpretation_figures(factors, loadings, clusters); pathway_and_candidates(loadings, shortlist); robustness(factors, clusters)
    (RUN / "MOFA_REVISED_11_REPORT.txt").write_text("11-patient PDAC MOFA integration\n\nSamples: " + ", ".join(PATIENTS) + "\nViews: blacklist-filtered WES, tumour-normal RNA log2 delta, log2 proteomics abundance.\nMOFA settings retained: 5 factors, 1500 iterations, seed 42, fast convergence.\nAll outputs are isolated under this directory; legacy 14-patient integration outputs were not overwritten.\nThe 11-patient model is exploratory and requires robustness-aware interpretation.\n")
    downstream_figures(factors, loadings, shortlist, rna)
    print("Completed isolated 11-patient MOFA workflow:", RUN)


if __name__ == "__main__":
    main()
