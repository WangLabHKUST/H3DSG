# -*- coding: utf-8 -*-
"""Share of cells in each of the four areas of the entropy map, PL against TL (boxplots with points).

Input:  ../data/entropy_region_cell_proportion_by_sample.csv
Output: entropy_region_cell_proportion_boxplot_PL_vs_TL.png (written into this directory)

Reads ../data/ only; never touches the raw data. No pipeline step reproduces the input table: it
comes from the original analysis (cells of each sample sorted into four areas of the 2-D entropy
map, area 'none' = outside them) and was copied into data/ with SP IDs as the sample key.

One panel per area (niche1..niche4 in the table; the names below are the ones printed on the
published figure), one point per sample (31 samples). prop_nicheN = cells in area N / all cells
counted for the sample (the four areas plus 'none').

The p value is a two-sided Mann-Whitney U test as R's wilcox.test() computes it: the exact
distribution when there are no tied values, the normal approximation with tie and continuity
correction when there are (zeros are tied). Sample as the unit of analysis, PL (n = 14) against
TL (n = 17), no multiple-testing correction.

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
OUT = Path(__file__).resolve().parent / 'entropy_region_cell_proportion_boxplot_PL_vs_TL.png'

# column suffix in the table -> name printed on the published figure
NICHES = {'niche1': 'Neuronal niche', 'niche2': 'Inflammatory niche',
          'niche3': 'Proliferative niche', 'niche4': 'Immune niche'}

# seaborn's stripplot draws its jitter offsets from numpy's global RNG; seeding it makes
# this figure byte-reproducible. Same value as JITTER_SEED in panel_helpers.py.
JITTER_SEED = 42


def wilcox_p(pl, tl):
    """Two-sided Mann-Whitney U p value, following R's wilcox.test() default."""
    ties = len(np.unique(np.r_[pl, tl])) < len(pl) + len(tl)
    return stats.mannwhitneyu(pl, tl, alternative='two-sided', method='asymptotic' if ties else 'exact').pvalue


def main():
    src = DATA / 'entropy_region_cell_proportion_by_sample.csv'
    if not src.exists():
        raise FileNotFoundError(f'{src} does not exist; it is copied from the original analysis, see the header.')

    df = pd.read_csv(src)
    df['group_label'] = pd.Categorical(df['group'].map(group_labels), categories=group_order, ordered=True)

    np.random.seed(JITTER_SEED)        # see JITTER_SEED above
    fig, axes = plt.subplots(1, len(NICHES), figsize=(2.3 * len(NICHES), 3.0), sharey=False)

    for ax, (key, name) in zip(axes, NICHES.items()):
        col = f'prop_{key}'
        p = wilcox_p(df.loc[df['group_label'] == 'PL', col].to_numpy(), df.loc[df['group_label'] == 'TL', col].to_numpy())
        sns.boxplot(
            data=df, x='group_label', y=col, hue='group_label',
            order=group_order, hue_order=group_order, palette=group_colors,
            width=0.55, linewidth=1.0, showfliers=False, legend=False, ax=ax,
            medianprops={'color': 'black', 'linewidth': 1.6},
        )
        sns.stripplot(data=df, x='group_label', y=col, order=group_order,
                      color='black', size=2.6, jitter=0.18, ax=ax)
        ax.set_title(f'{name}\n$P$ = {p:.1e}', fontsize=9)
        ax.set_xlabel('')
        ax.set_ylabel('Proportion', fontsize=9)
        ax.tick_params(axis='x', labelsize=9)
        ax.tick_params(axis='y', labelsize=8)
        sns.despine(ax=ax, top=True, right=True)

    fig.patch.set_facecolor('white')
    plt.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'{len(df)} samples x {len(NICHES)} areas')
    print(f'saved -> {OUT}')


if __name__ == '__main__':
    main()
