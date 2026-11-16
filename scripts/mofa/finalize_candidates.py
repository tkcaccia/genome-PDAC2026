#!/usr/bin/env python3
# MOFA fitting, molecular integration and publication panels.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""

from __future__ import annotations

from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path

from pathlib import Path
import math
import textwrap
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


ROOT = private_path("handoff_root")
BIO = ROOT / 'mofa' / 'pdac_11' / 'mofa' / 'biology'
OUT = ROOT / 'mofa' / 'pdac_11' / 'mofa' / 'final'
OUT.mkdir(parents=True, exist_ok=True)


def parse_genes_field(s: str) -> set[str]:
    if pd.isna(s):
        return set()
    return {x.strip() for x in str(s).split(';') if x.strip()}


def select_pathways(df: pd.DataFrame) -> pd.DataFrame:
    # Primary threshold
    sig = df[df['Adjusted P-value'] < 0.10].copy()
    if len(sig) == 0:
        # Fallback for small cohorts / weak signal: top terms by combined score per factor+gene_set
        sig = (
            df.sort_values('Combined Score', ascending=False)
              .groupby(['factor', 'gene_set'], as_index=False)
              .head(8)
              .copy()
        )
        sig['selection_rule'] = 'fallback_top_combined_score'
    else:
        sig['selection_rule'] = 'adjusted_pvalue_lt_0.10'

    sig = sig.sort_values(['factor', 'gene_set', 'Adjusted P-value', 'Combined Score'], ascending=[True, True, True, False])
    return sig


def save_pathway_figure(pathways: pd.DataFrame):
    if len(pathways) == 0:
        return

    # Keep the same top-six terms per factor, but give each factor a separate axis
    # so long pathway names cannot spill into neighboring panels.
    plot_df = pathways.sort_values('Combined Score', ascending=False).groupby('factor', as_index=False).head(6).copy()
    plot_df['neglog10_adj_p'] = -plot_df['Adjusted P-value'].clip(lower=1e-300).map(lambda x: __import__('math').log10(x))
    factors = sorted(plot_df['factor'].unique())
    ncols = 2
    nrows = math.ceil(len(factors) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 5.2 * nrows), squeeze=False)
    colors = {'GO_Biological_Process_2023': '#457b9d', 'Reactome_2022': '#e58b2a', 'KEGG_2021_Human': '#2a9d8f'}
    for index, factor in enumerate(factors):
        ax = axes[index // ncols][index % ncols]
        subset = plot_df[plot_df['factor'].eq(factor)].sort_values('Combined Score', ascending=True)
        labels = [textwrap.fill(str(term), width=34, break_long_words=False, break_on_hyphens=False) for term in subset['Term']]
        y = range(len(subset))
        ax.barh(y, subset['neglog10_adj_p'], color=[colors.get(name, '#777777') for name in subset['gene_set']])
        ax.set_yticks(list(y), labels, fontsize=9)
        ax.set_title(factor.replace('_', ' ').title())
        ax.set_xlabel('-log10(adjusted p-value)')
        ax.grid(axis='x', alpha=.25, linestyle=':')
        ax.set_axisbelow(True)
    for index in range(len(factors), nrows * ncols):
        axes[index // ncols][index % ncols].set_visible(False)
    handles = [Patch(color=color, label=gene_set) for gene_set, color in colors.items() if gene_set in set(plot_df['gene_set'])]
    fig.legend(handles=handles, title='Gene set', loc='lower center', ncol=max(1, len(handles)), frameon=False)
    fig.suptitle('Top pathway terms per MOFA factor', fontsize=16)
    fig.subplots_adjust(left=.24, right=.98, top=.90, bottom=.12, wspace=.85, hspace=.55)
    fig.savefig(OUT / 'mofa_pathway_summary.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def build_candidate_table(rank: pd.DataFrame, pathways: pd.DataFrame) -> pd.DataFrame:
    # Candidate pool: at least 2 modalities and in the top half by composite score
    cut = rank['composite_score'].quantile(0.50)
    cand = rank[(rank['n_modalities'] >= 2) & (rank['composite_score'] >= cut)].copy()

    pathway_gene_set = set()
    for s in pathways.get('Genes', pd.Series([], dtype=str)):
        pathway_gene_set.update(parse_genes_field(s))

    cand['pathway_hit'] = cand['feature_key'].astype(str).isin(pathway_gene_set)
    cand['evidence_modalities'] = cand[['wes_score', 'rnaseq_score', 'proteomics_score']].notna().sum(axis=1)

    # Prioritize cross-modal consistency and pathway support
    cand = cand.sort_values(
        ['pathway_hit', 'evidence_modalities', 'composite_score'],
        ascending=[False, False, False],
    )

    # Final shortlist size
    top = cand.head(120).copy()
    return top


def main():
    enrich = pd.read_csv(BIO / 'mofa_pathway_enrichment.csv')
    rank = pd.read_csv(BIO / 'mofa_cross_modal_target_ranking.csv')

    selected = select_pathways(enrich)
    selected.to_csv(OUT / 'mofa_pathway_selected.csv', index=False)

    # Aggregated per-factor summary table
    fac_summary = (
        selected.groupby('factor', as_index=False)
        .agg(
            n_selected_terms=('Term', 'count'),
            best_adjusted_p=('Adjusted P-value', 'min'),
            best_combined_score=('Combined Score', 'max'),
        )
        .sort_values('factor')
    )
    fac_summary.to_csv(OUT / 'mofa_factor_pathway_summary.csv', index=False)

    save_pathway_figure(selected)

    final_candidates = build_candidate_table(rank, selected)
    final_candidates.to_csv(OUT / 'top_candidates_for_validation.csv', index=False)

    report = [
        '# Final candidate package',
        '',
        f'- Selected pathways: `results/mofa/final/mofa_pathway_selected.csv` ({len(selected)} rows)',
        f'- Factor pathway summary: `results/mofa/final/mofa_factor_pathway_summary.csv`',
        f'- Pathway summary figure: `results/mofa/final/mofa_pathway_summary.png`',
        f'- Validation shortlist: `results/mofa/final/top_candidates_for_validation.csv` ({len(final_candidates)} rows)',
        '',
        'Selection logic:',
        '- Pathways: adjusted p-value < 0.10; fallback to top combined-score terms per factor+gene-set if none pass.',
        '- Candidates: at least 2 modalities + top 50% composite score, prioritized by pathway-hit and cross-modal evidence.',
    ]
    (OUT / 'FINAL_CANDIDATE_REPORT.md').write_text('\n'.join(report) + '\n')

    print('Wrote final outputs to', OUT)
    print('Selected pathways:', len(selected))
    print('Validation candidates:', len(final_candidates))


if __name__ == '__main__':
    main()
