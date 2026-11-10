# -*- coding: utf-8 -*-
"""Interacting cell pairs per niche inside two areas of the entropy map, PL against TL
(boxplots with points), six panels.

Input:  ../data/entropy_region_pairs_per_niche_area1.csv
        ../data/entropy_region_pairs_per_niche_area2.csv
Output: entropy_region_pairs_per_niche_boxplot_PL_vs_TL.png (written into this directory)

Reads ../data/ only; never touches the raw data. No pipeline step reproduces the input tables:
they come from the original analysis (per sample, the number of interacting cell pairs of each
cell type pair per niche inside area 1 and area 2 of the entropy map) and were copied into data/
with SP IDs as the sample key and the manuscript cell type names. The column
'A...B' holds the pairs between cell types A and B.

Left column: area 1 (two pairs); right two columns: area 2 (four pairs). One point per sample
(31 samples). The p value is a two-sided Mann-Whitney U test as R's wilcox.test() computes it:
the exact distribution when there are no tied values, the normal approximation with tie and
continuity correction when there are (zeros are tied). Sample as the unit of analysis, PL (n = 14)
against TL (n = 17), no multiple-testing correction.

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
from constants import celltype_display, group_colors, group_labels, group_order

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'entropy_region_pairs_per_niche_boxplot_PL_vs_TL.png'

# (area, column of the table, grid row, grid column)
PANELS = [
    (1, 'OPC_OC_like...OPC_OC', 0, 0),
    (1, 'AC_MES_like...Myeloid', 1, 0),
    (2, 'Cycling...OPC_OC_like', 0, 1),
    (2, 'Myeloid...T_NK', 0, 2),
    (2, 'Cycling...OPC_OC', 1, 1),
    (2, 'AC_MES_like...T_NK', 1, 2),
]

# seaborn's stripplot draws its jitter offsets from numpy's global RNG; seeding it makes
# this figure byte-reproducible. Same value as JITTER_SEED in panel_helpers.py.
JITTER_SEED = 42


def wilcox_p(pl, tl):
    """Two-sided Mann-Whitney U p value, following R's wilcox.test() default."""
    ties = len(np.unique(np.r_[pl, tl])) < len(pl) + len(tl)
    return stats.mannwhitneyu(pl, tl, alternative='two-sided', method='asymptotic' if ties else 'exact').pvalue


def pair_title(col):
    return '+'.join(celltype_display[c] for c in col.split('...'))


def main():
    tables = {}
    for area in (1, 2):
        src = DATA / f'entropy_region_pairs_per_niche_area{area}.csv'
        if not src.exists():
            raise FileNotFoundError(f'{src} does not exist; it is copied from the original analysis, see the header.')
        df = pd.read_csv(src)
        df['group_label'] = pd.Categorical(df['group'].map(group_labels), categories=group_order, ordered=True)
        tables[area] = df

    np.random.seed(JITTER_SEED)        # see JITTER_SEED above
    fig, axes = plt.subplots(2, 3, figsize=(2.3 * 3, 2.6 * 2))

    for area, col, i, j in PANELS:
        ax, df = axes[i, j], tables[area]
        p = wilcox_p(df.loc[df['group_label'] == 'PL', col].to_numpy(), df.loc[df['group_label'] == 'TL', col].to_numpy())
        sns.boxplot(
            data=df, x='group_label', y=col, hue='group_label',
            order=group_order, hue_order=group_order, palette=group_colors,
            width=0.55, linewidth=1.0, showfliers=False, legend=False, ax=ax,
            medianprops={'color': 'black', 'linewidth': 1.6},
        )
        sns.stripplot(data=df, x='group_label', y=col, order=group_order,
                      color='black', size=2.6, jitter=0.18, ax=ax)
        ax.set_title(f'{pair_title(col)}\n$P$ = {p:.1e}', fontsize=9)
        ax.set_xlabel('')
        ax.set_ylabel('Number of pairs per niche', fontsize=9)
        ax.tick_params(axis='x', labelsize=9)
        ax.tick_params(axis='y', labelsize=8)
        sns.despine(ax=ax, top=True, right=True)

    fig.patch.set_facecolor('white')
    plt.tight_layout(rect=(0, 0.06, 1, 1))
    left, right = axes[1, 0].get_position(), axes[1, 2].get_position()
    fig.text((left.x0 + left.x1) / 2, 0.02, 'Area 1', ha='center', va='center', fontsize=10)
    fig.text((axes[1, 1].get_position().x0 + right.x1) / 2, 0.02, 'Area 2', ha='center', va='center', fontsize=10)
    fig.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'{len(tables[1])} samples x {len(PANELS)} cell type pairs')
    print(f'saved -> {OUT}')


if __name__ == '__main__':
    main()
