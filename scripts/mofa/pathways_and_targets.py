#!/usr/bin/env python3
# MOFA fitting, molecular integration and publication panels.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""

from __future__ import annotations

from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path

from pathlib import Path
import re
import numpy as np
import pandas as pd
import mygene


ROOT = private_path("handoff_root")
MOFA_EXPORT = ROOT / 'mofa' / 'pdac_11' / 'mofa' / 'exports'
OUT = ROOT / 'mofa' / 'pdac_11' / 'mofa' / 'biology'
OUT.mkdir(parents=True, exist_ok=True)


def clean_ensg(x: str) -> str:
    x = str(x)
    return x.split('.')[0]


def map_ensembl_to_symbol(ids: list[str]) -> dict[str, str]:
    mg = mygene.MyGeneInfo()
    unique_ids = sorted(set([clean_ensg(x) for x in ids]))
    out = {}
    try:
        res = mg.querymany(unique_ids, scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True)
        if isinstance(res, pd.DataFrame):
            for q in unique_ids:
                sym = None
                if q in res.index:
                    row = res.loc[q]
                    if isinstance(row, pd.DataFrame):
                        row = row.iloc[0]
                    sym = row.get('symbol', None)
                out[q] = str(sym) if pd.notna(sym) and sym is not None else q
        else:
            out = {q: q for q in unique_ids}
    except Exception:
        out = {q: q for q in unique_ids}
    return out


def clean_uniprot(x: str) -> str:
    s = str(x)
    if ';' in s:
        s = s.split(';')[0]
    return s.strip()


def map_uniprot_to_symbol(ids: list[str]) -> dict[str, str]:
    mg = mygene.MyGeneInfo()
    unique_ids = sorted(set([clean_uniprot(x) for x in ids]))
    out = {}
    try:
        res = mg.querymany(unique_ids, scopes='uniprot', fields='symbol', species='human', as_dataframe=True)
        if isinstance(res, pd.DataFrame):
            for q in unique_ids:
                sym = None
                if q in res.index:
                    row = res.loc[q]
                    if isinstance(row, pd.DataFrame):
                        row = row.iloc[0]
                    sym = row.get('symbol', None)
                out[q] = str(sym) if pd.notna(sym) and sym is not None else q
        else:
            out = {q: q for q in unique_ids}
    except Exception:
        out = {q: q for q in unique_ids}
    return out


def run_factor_enrichment(rna_loadings: pd.DataFrame) -> pd.DataFrame:
    try:
        import gseapy as gp
    except Exception:
        return pd.DataFrame(columns=['factor', 'gene_set', 'Term', 'Adjusted P-value', 'Combined Score', 'Genes'])

    ensg_map = map_ensembl_to_symbol(rna_loadings.index.tolist())
    mapped_index = [ensg_map.get(clean_ensg(x), clean_ensg(x)) for x in rna_loadings.index]
    rna_tmp = rna_loadings.copy()
    rna_tmp.index = mapped_index

    all_rows = []
    gene_sets = ['GO_Biological_Process_2023', 'KEGG_2021_Human', 'Reactome_2022']

    for factor in rna_tmp.columns:
        s = rna_tmp[factor].abs().sort_values(ascending=False)
        top_genes = [g for g in s.head(300).index.tolist() if re.match(r'^[A-Za-z0-9_.-]+$', str(g))]
        top_genes = [g for g in top_genes if not str(g).startswith('ENSG')]
        top_genes = list(dict.fromkeys(top_genes))
        if len(top_genes) < 10:
            continue

        for gs in gene_sets:
            try:
                enr = gp.enrichr(gene_list=top_genes, gene_sets=gs, organism='human', outdir=None, cutoff=1.0)
                if enr is None or enr.results is None or len(enr.results) == 0:
                    continue
                df = enr.results.copy()
                keep_cols = [c for c in ['Term', 'Adjusted P-value', 'Combined Score', 'Genes'] if c in df.columns]
                df = df[keep_cols].head(25)
                df.insert(0, 'gene_set', gs)
                df.insert(0, 'factor', factor)
                all_rows.append(df)
            except Exception:
                continue

    if len(all_rows) == 0:
        return pd.DataFrame(columns=['factor', 'gene_set', 'Term', 'Adjusted P-value', 'Combined Score', 'Genes'])
    return pd.concat(all_rows, ignore_index=True)


def max_abs_loading_by_feature(df: pd.DataFrame, feature_name: str) -> pd.DataFrame:
    x = df.abs().max(axis=1).to_frame('max_abs_loading').reset_index()
    x = x.rename(columns={'index': feature_name})
    return x


def main() -> None:
    wes = pd.read_csv(MOFA_EXPORT / 'mofa_loadings_wes.csv', index_col=0)
    rna = pd.read_csv(MOFA_EXPORT / 'mofa_loadings_rnaseq.csv', index_col=0)
    prot = pd.read_csv(MOFA_EXPORT / 'mofa_loadings_proteomics.csv', index_col=0)

    # 1) Factor-wise enrichment from RNA loadings
    enr = run_factor_enrichment(rna)
    enr.to_csv(OUT / 'mofa_pathway_enrichment.csv', index=False)

    # 2) Cross-modal target ranking
    wes_tbl = max_abs_loading_by_feature(wes, 'feature_key').rename(columns={'max_abs_loading': 'wes_score'})

    rna_tbl = max_abs_loading_by_feature(rna, 'feature_key')
    ens_map = map_ensembl_to_symbol(rna_tbl['feature_key'].tolist())
    rna_tbl['feature_key'] = rna_tbl['feature_key'].map(lambda x: ens_map.get(clean_ensg(x), clean_ensg(x)))
    rna_tbl = rna_tbl.groupby('feature_key', as_index=False)['max_abs_loading'].max().rename(columns={'max_abs_loading': 'rnaseq_score'})

    prot_tbl = max_abs_loading_by_feature(prot, 'feature_key')
    uni_map = map_uniprot_to_symbol(prot_tbl['feature_key'].tolist())
    prot_tbl['feature_key'] = prot_tbl['feature_key'].map(lambda x: uni_map.get(clean_uniprot(x), clean_uniprot(x)))
    prot_tbl = prot_tbl.groupby('feature_key', as_index=False)['max_abs_loading'].max().rename(columns={'max_abs_loading': 'proteomics_score'})

    rank = wes_tbl.merge(rna_tbl, on='feature_key', how='outer').merge(prot_tbl, on='feature_key', how='outer')
    for c in ['wes_score', 'rnaseq_score', 'proteomics_score']:
        rank[c] = pd.to_numeric(rank[c], errors='coerce')

    rank['n_modalities'] = rank[['wes_score', 'rnaseq_score', 'proteomics_score']].notna().sum(axis=1)

    # percentile-normalized composite score
    for c in ['wes_score', 'rnaseq_score', 'proteomics_score']:
        rank[f'{c}_pct'] = rank[c].rank(pct=True)
    rank['composite_score'] = rank[['wes_score_pct', 'rnaseq_score_pct', 'proteomics_score_pct']].fillna(0).mean(axis=1)

    rank = rank.sort_values(['n_modalities', 'composite_score'], ascending=[False, False]).reset_index(drop=True)
    rank.to_csv(OUT / 'mofa_cross_modal_target_ranking.csv', index=False)
    rank.head(200).to_csv(OUT / 'mofa_cross_modal_target_top200.csv', index=False)

    report_lines = [
        '# MOFA pathway and target outputs',
        '',
        '- Pathway enrichment table: `results/mofa/biology/mofa_pathway_enrichment.csv`',
        '- Cross-modal target ranking: `results/mofa/biology/mofa_cross_modal_target_ranking.csv`',
        '- Top 200 ranked targets: `results/mofa/biology/mofa_cross_modal_target_top200.csv`',
        '',
        'Notes:',
        '- RNA features are mapped from Ensembl IDs to symbols when possible.',
        '- Proteomics features are mapped from UniProt accessions to symbols when possible.',
    ]
    (OUT / 'MOFA_BIOLOGY_REPORT.md').write_text('\n'.join(report_lines) + '\n')

    print('Wrote biology outputs to', OUT)
    print('Enrichment rows:', len(enr))
    print('Target ranking rows:', len(rank))


if __name__ == '__main__':
    main()
