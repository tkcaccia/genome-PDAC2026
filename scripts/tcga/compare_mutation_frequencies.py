# Current mutation-frequency comparison; no historical candidate-overlap analysis.
#!/usr/bin/env python3

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

ROOT = project_root()
RUN = ROOT / 'results' / 'mofa_revised_11_20260916'
FIG = RUN / 'mofa_figures'
TAB = RUN / 'tables'
COLORS = {'wes': '#264653', 'rnaseq': '#e76f51', 'proteomics': '#2a9d8f'}
GENE_SET_COLORS = {'GO_Biological_Process_2023': '#4C78A8', 'KEGG_2021_Human': '#F58518', 'Reactome_2022': '#54A24B'}


def revised_tcga_mutation():
    afr = pd.read_csv(ROOT / 'results/wes_revised_11_20260916/outputs/primary_pdac/wes_mutation_matrix.blacklist_filtered.csv', index_col=0).fillna(0).astype(int)
    tcga = pd.read_csv(ROOT / 'external/tcga_paad/processed/tcga_mc3_gene_binary.csv', index_col=0).fillna(0)
    common = afr.index.intersection(tcga.index)
    afr = afr.loc[common]; tcga = (tcga.loc[common] != 0).astype(int)
    rows = []
    for gene in common:
        a = int(afr.loc[gene].sum()); n1 = afr.shape[1]; b = int(tcga.loc[gene].sum()); n2 = tcga.shape[1]
        p = binomtest(a, n1, b / n2).pvalue if n2 else np.nan
        rows.append({'gene': gene, 'afr_k': a, 'afr_n': n1, 'tcga_k': b, 'tcga_n': n2, 'afr_freq': a/n1, 'tcga_freq': b/n2, 'p_value': p})
    out = pd.DataFrame(rows); out['FDR'] = multipletests(out.p_value.fillna(1), method='fdr_bh')[1]; out['freq_diff'] = out.afr_freq - out.tcga_freq; out = out.sort_values(['FDR','freq_diff'], ascending=[True,False]); out.to_csv(TAB / 'revised_11_pdac_tcga_mutation_frequency.csv', index=False)
    top = out.head(30).sort_values('freq_diff')
    fig, ax = plt.subplots(figsize=(8, 9)); ax.barh(top.gene, top.freq_diff * 100, color=np.where(top.freq_diff >= 0, '#e76f51', '#457b9d')); ax.axvline(0, color='black', lw=.7); ax.set(xlabel='11-patient PDAC frequency minus TCGA PAAD frequency (percentage points)', ylabel='Gene', title='Revised PDAC versus TCGA PAAD mutation frequency'); fig.tight_layout(); fig.savefig(FIG / 'revised_pdac_tcga_mutation_frequency.png', dpi=300); plt.close(fig)



def main():
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    revised_tcga_mutation()


if __name__ == "__main__":
    main()
