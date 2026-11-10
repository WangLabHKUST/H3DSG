# -*- coding: utf-8 -*-
"""Density-adjusted clustering coefficient, PL against TL (boxplots with points).

Input:  ../data/clustering_coefficient_adjusted.csv
        ../data/clustering_coefficient_adjusted_group_test.csv
        (both written by pipeline/07_clustering_coeff/01_density_adjusted.py)
Output: clustering_coeff_boxplot_PL_vs_TL.png (written into this directory)

Reads ../data/ only; never touches the raw data.

One panel per cell type, one point per sample. The y axis is log2: count_adj spans three
orders of magnitude across cell types, and on the raw scale the ten panels cannot be read
together. The p value comes from mw_p_adjusted in the group test table (two-sided
Mann-Whitney U, sample as the unit of analysis, no multiple-testing correction).

Boxplots never show outliers (showfliers=False): every box has jittered points drawn on
top, so if the boxplot also drew its outliers the same sample would appear twice -- once
at the outlier position on the box and once at its jittered offset -- and read as two
samples.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from constants import celltype_display, celltype_order, group_colors, group_labels, group_order

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'clustering_coeff_boxplot_PL_vs_TL.png'

NCOLS = 5

# seaborn's stripplot draws its jitter offsets from numpy's global RNG; seeding it makes
# this figure byte-reproducible. Same value as JITTER_SEED in panel_helpers.py.
JITTER_SEED = 42


def main():
    src = DATA / 'clustering_coefficient_adjusted.csv'
    test_src = DATA / 'clustering_coefficient_adjusted_group_test.csv'
    for fp in (src, test_src):
        if not fp.exists():
            raise FileNotFoundError(
                f'{fp} does not exist; run '
                'pipeline/07_clustering_coeff/01_density_adjusted.py first.')

    df = pd.read_csv(src)
    test = pd.read_csv(test_src).set_index('cell_type')

    df = df[df['count_adj'] > 0].copy()
    df['log2_count_adj'] = np.log2(df['count_adj'])
    df['group_label'] = df['group'].map(group_labels)
    df['group_label'] = pd.Categorical(df['group_label'], categories=group_order, ordered=True)

    types = [c for c in celltype_order if c in set(df['cell_type'])]
    nrows = int(np.ceil(len(types) / NCOLS))
    np.random.seed(JITTER_SEED)        # see JITTER_SEED above

    fig, axes = plt.subplots(nrows, NCOLS, figsize=(2.3 * NCOLS, 3.0 * nrows), sharey=False)
    axes = np.atleast_1d(axes).ravel()

    for ax, ct in zip(axes, types):
        sub = df[df['cell_type'] == ct]
        sns.boxplot(
            data=sub, x='group_label', y='log2_count_adj', hue='group_label',
            order=group_order, hue_order=group_order, palette=group_colors,
            width=0.55, linewidth=1.0, showfliers=False, legend=False, ax=ax,
            medianprops={'color': 'black', 'linewidth': 1.6},
        )
        sns.stripplot(data=sub, x='group_label', y='log2_count_adj', order=group_order,
                      color='black', size=2.6, jitter=0.18, ax=ax)
        p = test.loc[ct, 'mw_p_adjusted'] if ct in test.index else np.nan
        ax.set_title(f'{celltype_display.get(ct, ct)}\nP = {p:.3g}', fontsize=9)
        ax.set_xlabel('')
        ax.set_ylabel('log2(clustering coeff., adj.)' if ax is axes[0] else '')
        ax.tick_params(axis='x', labelsize=9)
        ax.tick_params(axis='y', labelsize=8)
        sns.despine(ax=ax, top=True, right=True)

    for ax in axes[len(types):]:
        ax.axis('off')

    fig.patch.set_facecolor('white')
    plt.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches='tight')
    plt.close(fig)

    sig = test[test['mw_p_adjusted'] < 0.05]
    print(f'{df["sample"].nunique()} samples x {len(types)} cell types')
    print(f'p < 0.05: {list(sig.index) if len(sig) else "none"}')
    print(f'saved -> {OUT}')


if __name__ == '__main__':
    main()
