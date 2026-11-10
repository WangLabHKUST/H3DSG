# -*- coding: utf-8 -*-
"""Pairwise enrichment: does the adjacency of two cell types exceed random expectation
(observed / expected)?

Input
  ../../raw_data/CODEX_adata_dict.npy      AnnData for 31 samples (archive version only)
  ../../sample_metadata.tsv                sample -> SP number, group

Output (committed)
  ../../data/pairwise_enrichment_by_sample.csv
        one row per sample per cell type pair: observed, expected, ratio, log2_ratio
  ../../data/pairwise_enrichment_pooled.csv
        observed and expected summed over the 31 samples and then divided (i.e. pooling
        weighted by edge count), plus log2(ratio) and log2(ratio + 1); the heatmap uses
        the latter

Method
  The neighbourhood graph is graph_utils.th_conn_matrix(adata, 'r50_adj', 5): starting
  from r50_adj, only each cell's 5 nearest neighbours are kept, which caps every cell's
  degree so that dense regions cannot dominate the counts. Each surviving edge (i, j) is
  counted once under the cell types of its two ends, giving the observed matrix; the
  expected value is row sum x column sum / total edges, the expected co-occurrence under
  an independence assumption.

Note that the observed matrix is symmetric (each undirected edge is counted once at
(i,j) and once at (j,i)), and so is expected, hence ratio is symmetric too; the upper
triangle of the heatmap is enough.

Differences from the original Interaction_inference.py (algorithm and numbers unchanged):
  * Statistics and plotting are separated: only tables are written here, and the heatmap
    lives in figures/interaction_enrichment_heatmap.py.
  * The sample list now comes from sample_metadata.tsv instead of a hard-coded
    g1_samples + g2_samples.
  * A few leftover notebook bare expressions were removed (`pair_df.shape`,
    `np.allclose(...)` and the like, whose return values were never used).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import celltype_order, load_samples
from data_loading import ADATA_DICT, load_adata_dict

OUT_PATH = REPO / 'data'

BASE_CONN_KEY = 'r50_adj'
N_NEIGHBORS = 5


def observed_expected(adata, sample_name):
    """Return the sample's (observed, expected) pair of 10x10 DataFrames."""
    from graph_utils import th_conn_matrix

    adata = th_conn_matrix(adata, BASE_CONN_KEY, N_NEIGHBORS)
    conn_matrix = adata.obsp[f'{BASE_CONN_KEY}_th{N_NEIGHBORS}_connectivities'].astype(int)
    conn_matrix.eliminate_zeros()

    celltype_df = adata.obs
    assert conn_matrix.shape[0] == celltype_df.shape[0], \
        f'{sample_name}: the graph and obs have different row counts'

    rows, cols = conn_matrix.nonzero()
    pair_df = pd.DataFrame({'i': rows, 'j': cols})
    pair_df['i_ct_agg'] = celltype_df.iloc[pair_df['i']]['ct_agg'].values
    pair_df['j_ct_agg'] = celltype_df.iloc[pair_df['j']]['ct_agg'].values

    ct_to_idx = {ct: idx for idx, ct in enumerate(celltype_order)}
    row_idx = pair_df['i_ct_agg'].map(ct_to_idx)
    col_idx = pair_df['j_ct_agg'].map(ct_to_idx)
    if row_idx.isna().any() or col_idx.isna().any():
        bad = sorted((set(pair_df['i_ct_agg']) | set(pair_df['j_ct_agg'])) - set(celltype_order))
        raise ValueError(f'{sample_name}: ct_agg has categories outside celltype_order: {bad}')

    observed = pd.crosstab(row_idx, col_idx).reindex(
        index=range(len(celltype_order)), columns=range(len(celltype_order)), fill_value=0)
    observed.index = celltype_order
    observed.columns = celltype_order

    n_edges = pair_df.shape[0]
    expected = pd.DataFrame(index=celltype_order, columns=celltype_order, dtype=float)
    for i, ct1 in enumerate(celltype_order):
        for j, ct2 in enumerate(celltype_order):
            expected.iloc[i, j] = observed.iloc[i, :].sum() * observed.iloc[:, j].sum() / n_edges

    return observed, expected


def to_long(mat, value_name):
    out = mat.stack().reset_index()
    out.columns = ['cell_type1', 'cell_type2', value_name]
    return out


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

    obs_all = pd.DataFrame(0.0, index=celltype_order, columns=celltype_order)
    exp_all = pd.DataFrame(0.0, index=celltype_order, columns=celltype_order)
    long_rows = []
    for sample_name in sample_lst:
        print(f'  {sample_name} ...', flush=True)
        observed, expected = observed_expected(adata_dict[sample_name], sample_name)
        obs_all += observed
        exp_all += expected

        one = to_long(observed.astype(float), 'observed').merge(
            to_long(expected, 'expected'), on=['cell_type1', 'cell_type2'])
        one['ratio'] = one['observed'] / one['expected']
        one['log2_ratio'] = np.log2(one['ratio'])
        one.insert(0, 'sample', sample_name)
        long_rows.append(one)

    by_sample = pd.concat(long_rows, ignore_index=True)
    by_sample = by_sample.merge(samples[['SP_ID', 'group']].rename_axis('sample').reset_index(),
                                on='sample', how='left')
    by_sample.to_csv(OUT_PATH / 'pairwise_enrichment_by_sample.csv', index=False)

    pooled = to_long(obs_all, 'observed').merge(to_long(exp_all, 'expected'),
                                                on=['cell_type1', 'cell_type2'])
    pooled['ratio'] = pooled['observed'] / pooled['expected']
    pooled['log2_ratio'] = np.log2(pooled['ratio'])
    ### The heatmap plots log2(ratio + 1): log2 runs to -inf as ratio approaches 0, and
    ### after adding 1, zero sits just below "exactly as expected", with negative values
    ### masked out in the heatmap (see the figures/ side)
    pooled['log2_ratio_plus1'] = np.log2(pooled['ratio'] + 1)
    pooled.to_csv(OUT_PATH / 'pairwise_enrichment_pooled.csv', index=False)

    print(f'\n{len(sample_lst)} samples x {len(celltype_order)**2} cell type pairs '
          f'= {len(by_sample)} rows')
    print('The 8 cell type pairs with the highest pooled log2(observed/expected):')
    print(pooled.sort_values('log2_ratio', ascending=False)
                .head(8)[['cell_type1', 'cell_type2', 'ratio', 'log2_ratio']]
                .round(3).to_string(index=False))
    print(f'\nsaved -> {OUT_PATH / "pairwise_enrichment_by_sample.csv"}')
    print(f'saved -> {OUT_PATH / "pairwise_enrichment_pooled.csv"}')


if __name__ == '__main__':
    main()
