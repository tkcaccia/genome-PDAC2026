#!/usr/bin/env python3
# PDAC2026 analysis: Expression-to-deconvolution scores, paired comparisons and tumour-only Hallmark/TME groups.
# Input formats and private configuration: docs/inputs.md.
# Methods and code identity: docs/workflow.md and docs/source_manifest.json.
"""Create an independently specified programme/TME profile for 11 PDAC tumours."""

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path

import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = project_root()
RUN = ROOT / "results/immune_salmon_11_20261012"
INPUT = RUN / "inputs/salmon_gene_symbol_tpm_11_pairs.csv"
HALLMARK_GMT = RUN / "inputs/figureD_msigdb_hallmark_2020_selected.gmt"
CACHED_HALLMARK_GMT = private_path("hallmark_gmt")
OUT_FIG = RUN / "figures/Figure_D_independent_programmes_TME_11_pdac.png"
TABLES = RUN / "tables"
PATIENTS = configured_cohort("primary_pdac")
TUMOUR_SAMPLES = [f"{patient}T" for patient in PATIENTS]

PROGRAMMES = [
    "Epithelial Mesenchymal Transition",
    "Hypoxia",
    "TGF-beta Signaling",
    "Angiogenesis",
    "KRAS Signaling Up",
    "KRAS Signaling Dn",
    "E2F Targets",
    "G2-M Checkpoint",
    "Myc Targets V1",
    "p53 Pathway",
    "Glycolysis",
    "Oxidative Phosphorylation",
    "Fatty Acid Metabolism",
    "Cholesterol Homeostasis",
    "Unfolded Protein Response",
    "Interferon Alpha Response",
    "Interferon Gamma Response",
    "Inflammatory Response",
    "Complement",
]

TME_FEATURES = [
    ("Immune", "ESTIMATE", "ImmuneScore"),
    ("Immune", "MCP-counter", "T cells"),
    ("Immune", "MCP-counter", "CD8 T cells"),
    ("Immune", "MCP-counter", "Cytotoxic lymphocytes"),
    ("Immune", "MCP-counter", "B lineage"),
    ("Immune", "MCP-counter", "NK cells"),
    ("Immune", "MCP-counter", "Monocytic lineage"),
    ("Immune", "MCP-counter", "Myeloid dendritic cells"),
    ("Immune", "MCP-counter", "Neutrophils"),
    ("Stromal", "ESTIMATE", "StromalScore"),
    ("Stromal", "MCP-counter", "Endothelial cells"),
    ("Stromal", "MCP-counter", "Fibroblasts"),
]

SCORE_FILES = {
    "ESTIMATE": "estimate_native_scores_11_pairs.csv",
    "MCP-counter": "mcp_counter_native_scores_11_pairs.csv",
}


def load_hallmark_sets():
    # Gene membership comes from Hallmark 2020, not from the patient's matrix.
    if HALLMARK_GMT.exists():
        gene_sets = {}
        with HALLMARK_GMT.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                gene_sets[fields[0]] = fields[2:]
    else:
        try:
            library = gp.get_library("MSigDB_Hallmark_2020", organism="Human")
        except Exception:
            if not CACHED_HALLMARK_GMT.exists():
                raise
            library = {}
            with CACHED_HALLMARK_GMT.open() as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    library[fields[0]] = fields[1:]
        missing = sorted(set(PROGRAMMES) - set(library))
        if missing:
            raise ValueError(f"Hallmark collection is missing selected sets: {missing}")
        gene_sets = {name: list(library[name]) for name in PROGRAMMES}
        HALLMARK_GMT.parent.mkdir(parents=True, exist_ok=True)
        with HALLMARK_GMT.open("w") as handle:
            for name in PROGRAMMES:
                handle.write("\t".join([name, "MSigDB_Hallmark_2020", *gene_sets[name]]) + "\n")
    missing = sorted(set(PROGRAMMES) - set(gene_sets))
    if missing:
        raise ValueError(f"Cached Hallmark file is missing selected sets: {missing}")
    return {name: gene_sets[name] for name in PROGRAMMES}


def zscore_rows(frame):
    means = frame.mean(axis=1, skipna=True)
    standard_deviations = frame.std(axis=1, skipna=True, ddof=1).replace(0, np.nan)
    return frame.sub(means, axis=0).div(standard_deviations, axis=0).fillna(0)


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    (RUN / "figures").mkdir(parents=True, exist_ok=True)
    tpm = pd.read_csv(INPUT, index_col=0)
    tpm.index = tpm.index.astype(str).str.upper()
    tpm = tpm.apply(pd.to_numeric, errors="coerce").fillna(0)
    missing_samples = sorted(set(TUMOUR_SAMPLES) - set(tpm.columns))
    if missing_samples:
        raise ValueError(f"Missing tumour samples: {missing_samples}")
    tumour_expression = np.log2(tpm[TUMOUR_SAMPLES] + 1)
    gene_sets = load_hallmark_sets()
    expression_genes = set(tumour_expression.index)
    overlap = pd.DataFrame([
        {"programme": name, "source_set_size": len(genes), "genes_in_expression_matrix": len(set(genes) & expression_genes)}
        for name, genes in gene_sets.items()
    ])
    if (overlap["genes_in_expression_matrix"] < 15).any():
        raise ValueError("At least one selected Hallmark set has fewer than 15 mapped genes")
    overlap.to_csv(TABLES / "figureD_programme_gene_set_overlap_11_pdac.csv", index=False)

    ssgsea = gp.ssgsea(
        data=tumour_expression,
        gene_sets=gene_sets,
        min_size=15,
        max_size=500,
        sample_norm_method="rank",
        correl_norm_type="rank",
        weight=0.25,
        threads=1,
        outdir=None,
        no_plot=True,
        seed=42,
        verbose=False,
    )
    score_result = ssgsea.res2d
    programme_scores = score_result.pivot(index="Term", columns="Name", values="ES").reindex(PROGRAMMES)[TUMOUR_SAMPLES]
    programme_scores.to_csv(TABLES / "figureD_hallmark_ssgsea_scores_11_pdac.csv", index_label="programme")
    programme_z = zscore_rows(programme_scores)
    programme_z.to_csv(TABLES / "figureD_hallmark_ssgsea_row_zscores_11_pdac.csv", index_label="programme")

    score_matrices = {}
    for method, filename in SCORE_FILES.items():
        frame = pd.read_csv(RUN / "tables" / filename, index_col=0)
        frame.columns = frame.columns.astype(str)
        score_matrices[method] = frame.apply(pd.to_numeric, errors="coerce")
    raw_tme = pd.DataFrame(index=[f"{method}: {feature}" for _, method, feature in TME_FEATURES], columns=TUMOUR_SAMPLES, dtype=float)
    definitions = []
    for domain, method, feature in TME_FEATURES:
        row_label = f"{method}: {feature}"
        matrix = score_matrices[method]
        if feature not in matrix.index:
            raise ValueError(f"Missing TME feature: {row_label}")
        for sample in TUMOUR_SAMPLES:
            if sample not in matrix.columns or pd.isna(matrix.loc[feature, sample]):
                raise ValueError(f"Missing complete-cohort value for {row_label}/{sample}")
            else:
                raw_tme.loc[row_label, sample] = matrix.loc[feature, sample]
        definitions.append({
            "domain": domain,
            "method": method,
            "feature": feature,
            "display_label": row_label,
            "scale_note": "Native score/fraction; row-standardized across observed PDAC tumours for display only.",
        })
    raw_tme.to_csv(TABLES / "figureD_tme_native_feature_values_11_pdac.csv", index_label="feature")
    pd.DataFrame(definitions).to_csv(TABLES / "figureD_tme_feature_definitions_11_pdac.csv", index=False)
    tme_z = zscore_rows(raw_tme)
    tme_z.to_csv(TABLES / "figureD_tme_feature_row_zscores_11_pdac.csv", index_label="feature")

    composite_scores = pd.DataFrame(index=["Immune composite", "Stromal composite", "Hallmark EMT"], columns=TUMOUR_SAMPLES, dtype=float)
    # Equal feature weights: nine immune rows, three stromal rows; EMT is separate.
    for domain in ["Immune", "Stromal"]:
        labels = [f"{method}: {feature}" for feature_domain, method, feature in TME_FEATURES if feature_domain == domain]
        composite_scores.loc[f"{domain} composite"] = tme_z.loc[labels].mean(axis=0)
    composite_scores.loc["Hallmark EMT"] = programme_z.loc["Epithelial Mesenchymal Transition"]
    composite_z = zscore_rows(composite_scores)
    composite_scores.to_csv(TABLES / "figureD_tme_composite_scores_11_pdac.csv", index_label="score")
    composite_z.to_csv(TABLES / "figureD_tme_composite_row_zscores_11_pdac.csv", index_label="score")

    thresholds = composite_scores.median(axis=1)
    # Cohort-relative tumour-only labels; ties are high and EMT does not restrict
    # the immune-high/stromal-low branch. Paired tests are a separate analysis.
    patient_groups = {}
    for sample in TUMOUR_SAMPLES:
        immune_high = composite_scores.loc["Immune composite", sample] >= thresholds["Immune composite"]
        stromal_high = composite_scores.loc["Stromal composite", sample] >= thresholds["Stromal composite"]
        emt_high = composite_scores.loc["Hallmark EMT", sample] >= thresholds["Hallmark EMT"]
        if stromal_high and emt_high and not immune_high:
            group = "Stromal/EMT-high, immune-low"
        elif immune_high and not stromal_high:
            group = "Immune-high, stromal-low"
        else:
            group = "Mixed/intermediate"
        patient_groups[sample] = group
    pd.DataFrame([
        {"score": name, "median_cutoff": float(thresholds[name]), "rule_high": ">= cohort median", "rule_low": "< cohort median"}
        for name in composite_scores.index
    ]).to_csv(TABLES / "figureD_tme_median_cutoffs_11_pdac.csv", index=False)
    pd.DataFrame([
        {"sample": sample, "patient": sample[:-1], "immune_composite": composite_scores.loc["Immune composite", sample],
         "stromal_composite": composite_scores.loc["Stromal composite", sample], "hallmark_emt": composite_scores.loc["Hallmark EMT", sample],
         "exploratory_group": patient_groups[sample]}
        for sample in TUMOUR_SAMPLES
    ]).to_csv(TABLES / "figureD_tme_exploratory_groups_11_pdac.csv", index=False)

    sns.set_theme(style="white", font="DejaVu Sans")
    fig = plt.figure(figsize=(15, 22), constrained_layout=False)
    grid = fig.add_gridspec(3, 2, width_ratios=[40, 1], height_ratios=[19, 12, 10], hspace=0.34, wspace=0.08)
    ax_programmes = fig.add_subplot(grid[0, 0])
    cax_programmes = fig.add_subplot(grid[0, 1])
    ax_tme = fig.add_subplot(grid[1, 0])
    cax_tme = fig.add_subplot(grid[1, 1])
    ax_composite = fig.add_subplot(grid[2, 0])
    cax_composite = fig.add_subplot(grid[2, 1])
    cmap = sns.color_palette("vlag", as_cmap=True)
    mask_programmes = programme_z.isna()
    mask_tme = tme_z.isna()
    sns.heatmap(programme_z, ax=ax_programmes, cmap=cmap, center=0, vmin=-2, vmax=2,
                mask=mask_programmes, linewidths=.25, linecolor="white", cbar=True, cbar_ax=cax_programmes,
                cbar_kws={"label": "Within-feature SD units"})
    ax_programmes.set_title("A  Hallmark programme scores (ssGSEA)", loc="left", fontsize=12, fontweight="bold", pad=10)
    ax_programmes.set_ylabel("")
    ax_programmes.set_xlabel("")
    ax_programmes.tick_params(axis="y", labelsize=9, length=0)
    ax_programmes.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_programmes.set_yticklabels(ax_programmes.get_yticklabels(), rotation=0)

    sns.heatmap(tme_z, ax=ax_tme, cmap=cmap, center=0, vmin=-2, vmax=2,
                mask=mask_tme, linewidths=.25, linecolor="white", cbar=True, cbar_ax=cax_tme,
                cbar_kws={"label": "Within-feature SD units"})
    ax_tme.set_title("B  Immune and stromal estimates", loc="left", fontsize=12, fontweight="bold", pad=10)
    ax_tme.set_ylabel("")
    ax_tme.set_xlabel("PDAC tumour sample")
    ax_tme.tick_params(axis="y", labelsize=9, length=0)
    ax_tme.tick_params(axis="x", labelrotation=45, labelsize=9)
    ax_tme.set_yticklabels(ax_tme.get_yticklabels(), rotation=0)
    ax_tme.axhline(9, color="black", linewidth=1.2)
    for axis in (ax_programmes, ax_tme):
        axis.set_xticklabels(TUMOUR_SAMPLES if axis is ax_tme else [], rotation=45 if axis is ax_tme else 0, ha="right")

    cax_composite.set_axis_off()
    group_colors = {
        "Stromal/EMT-high, immune-low": "#d1495b",
        "Immune-high, stromal-low": "#00798c",
        "Mixed/intermediate": "#edae49",
    }
    for group, color in group_colors.items():
        samples = [sample for sample in TUMOUR_SAMPLES if patient_groups[sample] == group]
        if samples:
            ax_composite.scatter(
                [composite_scores.loc["Immune composite", sample] for sample in samples],
                [composite_scores.loc["Stromal composite", sample] for sample in samples],
                color=color, s=65, label=group, zorder=3,
            )
    ax_composite.axvline(thresholds["Immune composite"], color="#444444", linestyle="--", linewidth=1)
    ax_composite.axhline(thresholds["Stromal composite"], color="#444444", linestyle="--", linewidth=1)
    for sample in TUMOUR_SAMPLES:
        ax_composite.annotate(sample[:-1], (composite_scores.loc["Immune composite", sample], composite_scores.loc["Stromal composite", sample]),
                              xytext=(4, 3), textcoords="offset points", fontsize=8)
    ax_composite.set_title("C  Exploratory cohort-relative TME groups", loc="left", fontsize=12, fontweight="bold", pad=10)
    ax_composite.set_xlabel("Immune composite (mean feature z-score)")
    ax_composite.set_ylabel("Stromal composite (mean feature z-score)")
    ax_composite.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8, frameon=False)
    ax_composite.grid(color="#dddddd", linewidth=.5, zorder=0)

    fig.suptitle("Tumour programme and microenvironment profiles in 11 PDAC cases", fontsize=15, fontweight="bold", y=.99)
    fig.text(.01, .005,
             "Tumour programme and microenvironment profiles. Panel A shows 19 selected Hallmark 2020 ssGSEA programmes. Panel B shows 12 ESTIMATE/MCP-counter features, row-standardized for display. Panel C classifies continuous immune, stromal and Hallmark EMT scores using each score's 11-tumour median; these cohort-relative groups are descriptive, not validated subtypes. MCP-counter NK signature coverage is limited (4/9 markers).",
             fontsize=8, va="bottom", wrap=True)
    fig.subplots_adjust(left=.30, right=.82, top=.965, bottom=.055)
    fig.savefig(OUT_FIG, dpi=300, bbox_inches="tight")
    plt.close(fig)

    methods = [
        "Cohort: tumour samples in the configured primary PDAC cohort; periampullary cases excluded.",
        "Expression: summed gene-symbol Salmon TPM from the immune-analysis input; log2(TPM+1) used for ssGSEA. The transform preserves within-sample gene ranks.",
        "Pathways: 19 selected sets from MSigDB_Hallmark_2020 retrieved with GSEApy 1.3.0; exact selected gene membership is archived in figureD_msigdb_hallmark_2020_selected.gmt. Sets were selected to represent EMT/stromal signalling, growth, metabolism and immune/inflammatory processes.",
        "Scoring: GSEApy ssGSEA, rank normalization, rank correlation normalization, weight 0.25, min_size 15, max_size 500, seed 42. Heatmap cells are row z-scores across observed tumours, not inferential p-values.",
        "TME: 12 features selected from ESTIMATE and MCP-counter, with complete values for all 11 tumours. Each feature is z-scored across the 11 tumours, then the immune and stromal composites are the unweighted means of their respective feature z-scores. Feature directions are positive for higher immune/stromal signal. Definitions are in figureD_tme_feature_definitions_11_pdac.csv.",
        "EMT: the selected Hallmark Epithelial Mesenchymal Transition ssGSEA score is used as a continuous programme score. For the descriptive group labels, each of the immune, stromal and EMT composite scores is classified high at or above its cohort median, and low below its cohort median. Stromal/EMT-high immune-low requires stromal and EMT high and immune low; immune-high stromal-low requires immune high and stromal low; all others are mixed/intermediate. Median cutoffs are saved in figureD_tme_median_cutoffs_11_pdac.csv.",
        "The 12 features comprise ESTIMATE ImmuneScore and StromalScore; MCP-counter T cells, CD8 T cells, Cytotoxic lymphocytes, B lineage, NK cells, Monocytic lineage, Myeloid dendritic cells, Neutrophils, Endothelial cells and Fibroblasts. MCP-counter NK coverage is limited (4/9 markers), so groupings involving the composite should be interpreted cautiously.",
        "Interpretation: programme scores and cohort-relative TME groups are descriptive and exploratory, not validated molecular subtypes.",
    ]
    (RUN / "METHODS_FIGURE_D_INDEPENDENT_11_PDac.txt").write_text("\n".join(methods) + "\n")
    print(f"figure={OUT_FIG}")
    print(f"programmes={len(programme_scores)} patients={len(TUMOUR_SAMPLES)} TME_features={len(raw_tme)} groups={pd.Series(patient_groups).value_counts().to_dict()}")
    print(f"overlap_min={overlap.genes_in_expression_matrix.min()} overlap_max={overlap.genes_in_expression_matrix.max()}")


if __name__ == "__main__":
    main()
