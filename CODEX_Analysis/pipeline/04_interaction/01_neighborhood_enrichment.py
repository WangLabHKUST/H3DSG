# -*- coding: utf-8 -*-
"""Neighborhood enrichment: run squidpy's neighbourhood enrichment per sample and compare
the clustering coefficient between the two groups.

Input
  ../../raw_data/CODEX_adata_dict.npy      AnnData for 31 samples (archive version only)
  ../../sample_metadata.tsv                sample -> SP number, group

Output (committed)
  ../../data/neighborhood_enrichment_by_sample.csv
        one row per sample per cell type pair: count (interaction strength after degree
        normalisation), zscore, p_two_sided, greater, less
  ../../data/clustering_coefficient_group_test.csv
        the same-type self pairs (cell_type1 == cell_type2) are the clustering
        coefficient; the two groups are compared with a two-sided Mann-Whitney U test

Method
  The neighbourhood graph is r40 (40 um radius). squidpy's nhood_enrichment gives count
  and zscore; count is then normalised by "the sum of the two types' degrees minus the
  number of interactions between them", which makes it independent of the sample's cell
  density:

      interaction[i,j] = count[i,j] / (deg[i] + deg[j] - count[i,j])

  The zscore comes from squidpy's permutation test and is converted to a p value here
  under a normal approximation.

--------------------------------------------------------------------------
Differences from the original Spatial_clustering.py -- that script could not run, and
four defects are fixed here:

1. `process_sample()` was defined but never called, so `df_lst` stayed empty and the
   `pd.concat(df_lst)` that follows was bound to raise
   ValueError: No objects to concatenate. It is now actually called per sample.
2. The group comparison downstream used `sub['sample']`, but `process_sample()` never
   wrote a `sample` column. It does now.
3. `zscore_df` was merged twice, producing `zscore_x` / `zscore_y` columns, after which
   taking `'zscore'` raised KeyError. It is merged once now.
4. The `p_one_sided` column was computed by passing `alternative='ryy-one_sided'`, but
   `zscore_to_pvalue()` only recognises two-sided / greater / less and any other branch
   fell through to the end of the function and returned None -- that column was never a
   valid p value (a leftover typo). The one-sided p values are already given by the
   `greater` and `less` columns, so **the `p_one_sided` column is dropped**, and
   `zscore_to_pvalue()` now raises on an alternative it does not recognise.

Because the original script could not run, **this line has no pre-existing output to
reconcile against** -- it never produced a result before.
--------------------------------------------------------------------------
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats
import squidpy as sq

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import celltype_order, load_samples
from data_loading import ADATA_DICT, attach_radius_graph, load_adata_dict

OUT_PATH = REPO / 'data'

TH = 40          # um, the neighbourhood radius
SEED = 42        # seed of squidpy's permutation test, the same value 06_niche_density uses
# The obsp of the adata_dict in raw_data/ carries only r50_adj, not r40, so the r40 graph
# is rebuilt here from obsm['spatial']. The rebuilt graph has been verified edge for edge
# identical to the r40_connectivities already present in
# {sample}_..._modified.h5ad (diff_nnz=0); see data_loading.verify_radius_graph.


def zscore_to_pvalue(zscore, alternative='two-sided'):
    if alternative == 'two-sided':
        return 2 * (1 - stats.norm.cdf(abs(zscore)))
    if alternative == 'greater':
        return 1 - stats.norm.cdf(zscore)
    if alternative == 'less':
        return stats.norm.cdf(zscore)
    raise ValueError(
        f"alternative must be two-sided / greater / less, got {alternative!r}")


def process_sample(adata, sample_name, th=TH):
    """Return the long-form neighbourhood enrichment table for one sample."""
    adata.obs['ct_agg'] = pd.Categorical(adata.obs['ct_agg'], categories=celltype_order)

    if f'r{th}_connectivities' not in adata.obsp:
        attach_radius_graph(adata, th)

    # neighborhood enrichment to get counts/zscores
    # The seed has to be passed explicitly: zscore comes from a label permutation, so
    # without a seed both zscore and the p values derived from it change on every run.
    # count itself is deterministic and unaffected.
    # show_progress_bar=False: progress bars for 31 samples x 1000 permutations would
    # flood the log with over a hundred thousand lines.
    sq.gr.nhood_enrichment(adata, cluster_key='ct_agg', connectivity_key=f'r{th}',
                           seed=SEED, show_progress_bar=False, copy=False)
    res = adata.uns['ct_agg_nhood_enrichment']
    count = res['count']
    zscore = res['zscore']

    # normalization by summed degrees (exclude self-pairs)
    cnvt_matrix_th = adata.obsp[f'r{th}_connectivities']
    node_degree = cnvt_matrix_th.sum(axis=1).A1
    adata.obs[f'node_degree_r{th}'] = node_degree
    celltype_degree = adata.obs.groupby('ct_agg', observed=False)[f'node_degree_r{th}'].sum()
    interaction_matrix = count.copy()
    normalization_matrix = celltype_degree.values.reshape(-1, 1) + celltype_degree.values.reshape(1, -1)
    normalization_matrix = normalization_matrix - interaction_matrix
    interaction_matrix = interaction_matrix / normalization_matrix

    cats = adata.obs['ct_agg'].cat.categories
    zscore_df = pd.DataFrame(zscore, index=cats, columns=cats).stack().reset_index().fillna(0)
    zscore_df.columns = ['cell_type1', 'cell_type2', 'zscore']

    out = pd.DataFrame(interaction_matrix, index=cats, columns=cats).stack().reset_index().fillna(0)
    out.columns = ['cell_type1', 'cell_type2', 'count']
    out = out.merge(zscore_df, on=['cell_type1', 'cell_type2'], how='left')

    out['p_two_sided'] = zscore_to_pvalue(out['zscore'], alternative='two-sided')
    out['greater'] = zscore_to_pvalue(out['zscore'], alternative='greater')
    out['less'] = zscore_to_pvalue(out['zscore'], alternative='less')
    out.insert(0, 'sample', sample_name)
    return out


def clustering_coefficient_test(all_df, samples):
    """Group comparison of the same-type self pairs (two-sided Mann-Whitney U)."""
    df = all_df[all_df['cell_type1'] == all_df['cell_type2']].copy()
    g1 = set(samples.index[samples['group'] == 'Group1'])
    g2 = set(samples.index[samples['group'] == 'Group2'])

    rows = []
    for cell_type in sorted(df['cell_type1'].unique()):
        sub = df[df['cell_type1'] == cell_type]
        v1 = sub.loc[sub['sample'].isin(g1), 'count'].dropna().values
        v2 = sub.loc[sub['sample'].isin(g2), 'count'].dropna().values
        if len(v1) == 0 or len(v2) == 0:
            rows.append({'cell_type': cell_type, 'n_Group1': len(v1), 'n_Group2': len(v2),
                         'median_Group1': np.nan, 'median_Group2': np.nan,
                         'p_value': np.nan, 'enriched_in': 'NA'})
            continue
        _, p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
        m1, m2 = float(np.median(v1)), float(np.median(v2))
        rows.append({'cell_type': cell_type, 'n_Group1': len(v1), 'n_Group2': len(v2),
                     'median_Group1': m1, 'median_Group2': m2, 'p_value': p,
                     'enriched_in': 'Group1' if m1 > m2 else 'Group2' if m1 < m2 else 'Equal'})
    return pd.DataFrame(rows).sort_values('p_value').reset_index(drop=True)


def main():
    OUT_PATH.mkdir(parents=True, exist_ok=True)
    if not ADATA_DICT.exists():
        raise FileNotFoundError(
            f'{ADATA_DICT} does not exist. This step needs the raw data, which only the '
            'archive version carries.')

    print('loading CODEX_adata_dict.npy ...', flush=True)
    adata_dict = load_adata_dict()      # keyed by SP ID

    samples = load_samples()
    extra = sorted(set(adata_dict) - set(samples.index))
    short = sorted(set(samples.index) - set(adata_dict))
    if extra or short:
        raise ValueError(
            f'adata_dict does not match the analysis manifest: extra {extra}, missing {short}')
    sample_lst = [s for s in samples.index if s in adata_dict]

    df_lst = []
    for sample_name in sample_lst:
        print(f'  {sample_name} ...', flush=True)
        df_lst.append(process_sample(adata_dict[sample_name], sample_name))

    all_df = pd.concat(df_lst, ignore_index=True)
    all_df = all_df.merge(samples[['SP_ID', 'group']].rename_axis('sample').reset_index(),
                          on='sample', how='left')
    all_df.to_csv(OUT_PATH / 'neighborhood_enrichment_by_sample.csv', index=False)
    print(f'\n{len(sample_lst)} samples x {len(celltype_order)**2} cell type pairs '
          f'= {len(all_df)} rows')

    test = clustering_coefficient_test(all_df, samples)
    test.to_csv(OUT_PATH / 'clustering_coefficient_group_test.csv', index=False)
    print('\nClustering coefficient, group comparison (two-sided Mann-Whitney U):')
    print(test.round(4).to_string(index=False))

    print(f'\nsaved -> {OUT_PATH / "neighborhood_enrichment_by_sample.csv"}')
    print(f'saved -> {OUT_PATH / "clustering_coefficient_group_test.csv"}')


if __name__ == '__main__':
    main()
