# -*- coding: utf-8 -*-
"""Expression summary of every cell type on every protein channel, for the dot plot.

Input
  ../../raw_data/CODEX_adata_dict.npy      AnnData for 31 samples (archive version only)

Output (committed)
  ../../data/marker_summary_by_celltype.csv
        one row per (cell_type, marker):
          mean_expression   the mean over all cells of that type on that channel
          frac_expressing   the fraction of cells of that type with > 0 on that channel
          n_cells           the number of cells of that type

In the dot plot, dot size is frac_expressing and colour is mean_expression normalised
within each channel (scanpy's standard_scale='var': each channel is linearly stretched to
[0, 1] separately). The normalisation is left to the plotting script; this table stores
the raw, unnormalised summary values.

What X means: the X of the AnnData is the QuPath measurement multiplied by that channel's
binary positivity call (see `filtered_data = raw_data * binary_data` in
pipeline/01_build_adata/00_create_adata.py), so a negative cell is exactly 0 on that
channel and frac_expressing is equivalent to "the fraction of positive cells".

Difference from the original Show_dotplot.py: the original did sc.concat into one
full-size AnnData and then called sc.pl.dotplot, so redrawing the figure meant loading the
7.9 GB npy every time. The summary values are now computed once into a 2 KB csv and the
plotting never touches the raw data.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import celltype_order, load_samples
from data_loading import ADATA_DICT, load_adata_dict

OUT_PATH = REPO / 'data'


def main():
    OUT_PATH.mkdir(parents=True, exist_ok=True)
    if not ADATA_DICT.exists():
        raise FileNotFoundError(
            f'{ADATA_DICT} does not exist. This step needs the raw data, which only the '
            'archive version carries; on the GitHub version, run figures/ off '
            '../../data/marker_summary_by_celltype.csv.')

    print('loading CODEX_adata_dict.npy ...', flush=True)
    adata_dict = load_adata_dict()      # keyed by SP ID

    samples = load_samples()
    sample_lst = [s for s in samples.index if s in adata_dict]
    if len(sample_lst) != len(adata_dict):
        raise ValueError(f'adata_dict does not match the analysis manifest: '
                         f'{sorted(set(adata_dict) - set(sample_lst))}')

    var_names = list(adata_dict[sample_lst[0]].var_names)
    print(f'{len(sample_lst)} samples, {len(var_names)} channels', flush=True)

    ### Accumulate sample by sample, rather than concatenating everything into one big matrix
    n_ct, n_var = len(celltype_order), len(var_names)
    sum_expr = np.zeros((n_ct, n_var))
    n_pos = np.zeros((n_ct, n_var))
    n_cells = np.zeros(n_ct)
    ct_index = {c: i for i, c in enumerate(celltype_order)}

    for sample_name in sample_lst:
        print(f'  {sample_name} ...', flush=True)
        adata = adata_dict[sample_name]
        if list(adata.var_names) != var_names:
            raise ValueError(f'the channel order of {sample_name} differs from the first sample')
        X = adata.X
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        labels = adata.obs['ct_agg'].astype(str).values
        unknown = sorted(set(labels) - set(celltype_order))
        if unknown:
            raise ValueError(f'unknown ct_agg categories in {sample_name}: {unknown}')
        for ct in np.unique(labels):
            i = ct_index[ct]
            sub = X[labels == ct]
            sum_expr[i] += sub.sum(axis=0)
            n_pos[i] += (sub > 0).sum(axis=0)
            n_cells[i] += sub.shape[0]

    rows = []
    for ct, i in ct_index.items():
        if n_cells[i] == 0:
            continue
        for j, marker in enumerate(var_names):
            rows.append({'cell_type': ct, 'marker': marker,
                         'mean_expression': sum_expr[i, j] / n_cells[i],
                         'frac_expressing': n_pos[i, j] / n_cells[i],
                         'n_cells': int(n_cells[i])})
    out = pd.DataFrame(rows)
    ### Fix the row order to celltype_order x var_names, so the table is easy to diff
    out['cell_type'] = pd.Categorical(out['cell_type'], categories=celltype_order, ordered=True)
    out['marker'] = pd.Categorical(out['marker'], categories=var_names, ordered=True)
    out = out.sort_values(['cell_type', 'marker']).reset_index(drop=True)
    out.to_csv(OUT_PATH / 'marker_summary_by_celltype.csv', index=False)

    print(f'\n{int(n_cells.sum()):,} cells in total, {out["cell_type"].nunique()} types x '
          f'{len(var_names)} channels = {len(out)} rows')
    print(f'saved -> {OUT_PATH / "marker_summary_by_celltype.csv"}')


if __name__ == '__main__':
    main()
