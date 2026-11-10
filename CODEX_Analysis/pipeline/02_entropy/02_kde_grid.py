# -*- coding: utf-8 -*-
"""Entropy step 2: 2-D kernel density estimation, exporting the grid data that
figures/entropy_density_by_group.R plots.

Input:  ../../data_full/entropy_cell_level_sampled.csv   (written by 01_cell_level_stats.py)
Output: ../../data/grid_data_entropy_tumor_vs_non_malig_entropy_Group1.csv
        ../../data/grid_data_entropy_tumor_vs_non_malig_entropy_Group2.csv
        ../../data/grid_data_entropy_tumor_vs_non_malig_entropy_diff.csv

x axis  entropy_tumor      = Shannon entropy of the malignant states within the 27.5 um
                             neighbourhood
y axis  non_malig_entropy  = Shannon entropy of the non-malignant cell types in the same
                             neighbourhood
colour  density            = gaussian_kde evaluated on an nbins x nbins grid
diff                       = zi(Group1) - zi(Group2)

The three grid CSVs have 400 rows each (a 20x20 grid). They are **already post-KDE and
contain no cell-level information**, so they can be committed, and they are the only
input the GitHub version needs to redraw these three figures.

Note: the input table is already filtered and down-sampled (equal cell counts per group),
so nothing is filtered here.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import group_labels
from data_loading import DATA_FULL, read_data_full

OUT_PATH = REPO / 'data'
FIG_PATH = REPO / 'figures'

VAR1 = 'entropy_tumor'
VAR2 = 'non_malig_entropy'
PAIR_KEY = f'{VAR1}_vs_{VAR2}'
NBINS = 20
GROUPS = ['Group1', 'Group2']


def kde_on_grid(x, y, nbins=NBINS):
    """Lay an nbins x nbins grid over the range of x/y and return (xi, yi, zi)."""
    assert not np.isnan(x).any(), 'x contains NaN'
    assert not np.isnan(y).any(), 'y contains NaN'
    assert (x >= 0).all() and (y >= 0).all(), 'entropy should not be negative'

    xi, yi = np.mgrid[x.min():x.max():nbins * 1j,
                      y.min():y.max():nbins * 1j]
    kernel = gaussian_kde(np.vstack([x, y]))
    zi = kernel(np.vstack([xi.flatten(), yi.flatten()]))
    return xi, yi, zi


def main():
    OUT_PATH.mkdir(parents=True, exist_ok=True)
    FIG_PATH.mkdir(parents=True, exist_ok=True)

    src = DATA_FULL / 'entropy_cell_level_sampled.csv'
    if not src.exists():
        raise FileNotFoundError(
            f'{src} does not exist. It is written by 01_cell_level_stats.py, which needs '
            'the raw data (archive version only); if all you want is to redraw the '
            'figures, the grid CSVs in ../../data/ are already there.')

    print('[1/3] reading the downstream entry table ...', flush=True)
    df = read_data_full('entropy_cell_level_sampled.csv', index_col=0)   # sample column: SP IDs
    print(f'      cells = {len(df)}', flush=True)
    print(df['group'].value_counts().to_string(), flush=True)

    print(f'[2/3] gaussian_kde per group (nbins={NBINS}) ...', flush=True)
    xi_dict, yi_dict, zi_dict = {}, {}, {}
    for g in GROUPS:
        sub = df[df['group'] == g]
        xi, yi, zi = kde_on_grid(sub[VAR1].values, sub[VAR2].values)
        xi_dict[g], yi_dict[g], zi_dict[g] = xi, yi, zi
        print(f'      {g}: n={len(sub)}  density range=[{zi.min():.4f}, {zi.max():.4f}]',
              flush=True)

        pd.DataFrame({VAR1: xi.ravel(), VAR2: yi.ravel(), 'density': zi.ravel()}) \
          .to_csv(OUT_PATH / f'grid_data_{PAIR_KEY}_{g}.csv', index=False)

    print('[3/3] computing and exporting diff ...', flush=True)
    # Each group's grid coordinates come from its own data range, so confirm the grids
    # agree before subtracting point by point
    if not (np.allclose(xi_dict['Group1'], xi_dict['Group2']) and
            np.allclose(yi_dict['Group1'], yi_dict['Group2'])):
        print("      [WARN] the two groups' grid coordinates differ; diff uses the "
              "Group1 grid (as the original implementation did)", flush=True)

    zi_diff = zi_dict['Group1'] - zi_dict['Group2']
    xi_a, yi_a = xi_dict['Group1'], yi_dict['Group1']
    pd.DataFrame({VAR1: xi_a.ravel(), VAR2: yi_a.ravel(), 'density': zi_diff.ravel()}) \
      .to_csv(OUT_PATH / f'grid_data_{PAIR_KEY}_diff.csv', index=False)
    print(f'      diff range=[{zi_diff.min():.4f}, {zi_diff.max():.4f}]', flush=True)

    # Quick preview (the publication figures are produced by
    # figures/entropy_density_by_group.R)
    fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(16, 4.5))
    for i, g in enumerate(GROUPS):
        zi = zi_dict[g]
        cf = axes[i].pcolormesh(xi_dict[g], yi_dict[g], zi.reshape(xi_dict[g].shape),
                                shading='gouraud', cmap=plt.cm.Reds, vmin=0, vmax=20)
        axes[i].contour(xi_dict[g], yi_dict[g], zi.reshape(xi_dict[g].shape),
                        levels=20, linewidths=0.1, colors='black')
        axes[i].set_title(group_labels[g]); axes[i].set_xlabel(VAR1); axes[i].set_ylabel(VAR2)
        fig.colorbar(cf, ax=axes[i])

    max_abs = np.max(np.abs(zi_diff))
    cf = axes[2].pcolormesh(xi_a, yi_a, zi_diff.reshape(xi_a.shape),
                            shading='gouraud', cmap=plt.cm.coolwarm,
                            vmin=-max_abs / 2, vmax=max_abs / 2)
    axes[2].set_title(f'{group_labels["Group1"]} - {group_labels["Group2"]}')
    axes[2].set_xlabel(VAR1); axes[2].set_ylabel(VAR2)
    fig.colorbar(cf, ax=axes[2])
    plt.tight_layout()
    fig.savefig(FIG_PATH / f'{PAIR_KEY}_preview.png', dpi=200, bbox_inches='tight')
    plt.close()

    print(f'\ndone; grid data in {OUT_PATH}, preview figure in {FIG_PATH}', flush=True)


if __name__ == '__main__':
    main()
