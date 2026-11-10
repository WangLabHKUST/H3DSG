# -*- coding: utf-8 -*-
"""Proportion of CD8 T cells among all cells, PL against TL (boxplot with points).

Input:  ../data/tnk_subset_proportion_by_sample.csv
Output: cd8_t_cell_proportion_boxplot_PL_vs_TL.png (written into this directory)

Reads ../data/ only; never touches the raw data. No pipeline step reproduces the input table: it
comes from the original analysis (per-sample proportions of the fine-grained cell types, taken
over all cells of the sample) and was copied into data/ with SP IDs as the sample key. CD8 T cells
are the Tcell_CD8 class of that annotation (cell_type_2 in the module); with CD4 T, Treg and NK
cells it adds up to the T_NK proportion of celltype_proportions.csv.

One point per sample (31 samples), drawn in percent of all cells of the sample.

The p value is a two-sided Mann-Whitney U test with the exact null distribution (there are no
tied values; R's wilcox.test(exact = TRUE)), sample as the unit of analysis, PL (n = 14) against
TL (n = 17).

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
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from constants import group_colors, group_labels, group_order

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'cd8_t_cell_proportion_boxplot_PL_vs_TL.png'

# seaborn's stripplot draws its jitter offsets from numpy's global RNG; seeding it makes
# this figure byte-reproducible. Same value as JITTER_SEED in panel_helpers.py.
JITTER_SEED = 42


def main():
    src = DATA / 'tnk_subset_proportion_by_sample.csv'
    if not src.exists():
        raise FileNotFoundError(f'{src} does not exist; it is copied from the original analysis, see the header.')

    df = pd.read_csv(src)
    df['percent'] = df['Tcell_CD8'] * 100
    df['group_label'] = pd.Categorical(df['group'].map(group_labels), categories=group_order, ordered=True)

    pl = df.loc[df['group_label'] == 'PL', 'percent'].to_numpy()
    tl = df.loc[df['group_label'] == 'TL', 'percent'].to_numpy()
    assert len(np.unique(np.r_[pl, tl])) == len(pl) + len(tl), 'tied values, exact test not valid'
    p = stats.mannwhitneyu(pl, tl, alternative='two-sided', method='exact').pvalue

    np.random.seed(JITTER_SEED)        # see JITTER_SEED above
    fig, ax = plt.subplots(figsize=(2.3, 3.0))
    sns.boxplot(
        data=df, x='group_label', y='percent', hue='group_label',
        order=group_order, hue_order=group_order, palette=group_colors,
        width=0.55, linewidth=1.0, showfliers=False, legend=False, ax=ax,
        medianprops={'color': 'black', 'linewidth': 1.6},
    )
    sns.stripplot(data=df, x='group_label', y='percent', order=group_order,
                  color='black', size=2.6, jitter=0.18, ax=ax)
    ax.set_title(f'$P$ = {p:.1e}', fontsize=9)
    ax.set_xlabel('')
    ax.set_ylabel('CD8 T cells / total cells (%)', fontsize=9)
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=8)
    sns.despine(ax=ax, top=True, right=True)

    fig.patch.set_facecolor('white')
    plt.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'{len(df)} samples; PL n = {len(pl)}, TL n = {len(tl)}; P = {p:.3e}')
    print(f'saved -> {OUT}')


if __name__ == '__main__':
    main()
