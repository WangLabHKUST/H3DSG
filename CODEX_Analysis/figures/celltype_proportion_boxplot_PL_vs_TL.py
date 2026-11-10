# -*- coding: utf-8 -*-
"""Cell type proportions, PL against TL (boxplots with points), six panels.

Input:  ../data/celltype_proportions.csv
        (written by pipeline/03_proportion/01_celltype_counts.py)
Output: celltype_proportion_boxplot_PL_vs_TL.png (written into this directory)

Reads ../data/ only; never touches the raw data.

One panel per cell type, one point per sample (31 samples). The proportion is the share of
that class among all cells of the sample (the 10 ct_agg classes). The left column holds the
three malignant states, the right column three non-malignant classes.

The p value is a two-sided Mann-Whitney U test with the exact null distribution (there are no
tied values), sample as the unit of analysis, PL (n = 14) against TL (n = 17), no
multiple-testing correction.

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
OUT = Path(__file__).resolve().parent / 'celltype_proportion_boxplot_PL_vs_TL.png'

# (column header, cell types in the order they are drawn, top to bottom)
COLUMNS = [
    ('Malignant', ['Cycling', 'OPC_OC_like', 'AC_MES_like']),
    ('Non-malignant', ['Myeloid', 'T_NK', 'OPC_OC']),
]

# seaborn's stripplot draws its jitter offsets from numpy's global RNG; seeding it makes
# this figure byte-reproducible. Same value as JITTER_SEED in panel_helpers.py.
JITTER_SEED = 42


def main():
    src = DATA / 'celltype_proportions.csv'
    if not src.exists():
        raise FileNotFoundError(
            f'{src} does not exist; run pipeline/03_proportion/01_celltype_counts.py first.')

    df = pd.read_csv(src)
    df['group_label'] = pd.Categorical(df['group'].map(group_labels), categories=group_order, ordered=True)

    nrows = max(len(types) for _, types in COLUMNS)
    np.random.seed(JITTER_SEED)        # see JITTER_SEED above
    fig, axes = plt.subplots(nrows, len(COLUMNS), figsize=(2.3 * len(COLUMNS), 2.1 * nrows), sharey=False,
                             squeeze=False)

    for j, (_, types) in enumerate(COLUMNS):
        for i, ct in enumerate(types):
            ax = axes[i, j]
            pl = df.loc[df['group_label'] == 'PL', ct].to_numpy()
            tl = df.loc[df['group_label'] == 'TL', ct].to_numpy()
            assert len(np.unique(np.r_[pl, tl])) == len(pl) + len(tl), f'{ct}: tied values, exact test not valid'
            p = stats.mannwhitneyu(pl, tl, alternative='two-sided', method='exact').pvalue

            sns.boxplot(
                data=df, x='group_label', y=ct, hue='group_label',
                order=group_order, hue_order=group_order, palette=group_colors,
                width=0.55, linewidth=1.0, showfliers=False, legend=False, ax=ax,
                medianprops={'color': 'black', 'linewidth': 1.6},
            )
            sns.stripplot(data=df, x='group_label', y=ct, order=group_order,
                          color='black', size=2.6, jitter=0.18, ax=ax)
            ax.set_title(f'$P$ = {p:.1e}', fontsize=9)
            ax.set_xlabel('')
            ax.set_ylabel(celltype_display.get(ct, ct), fontsize=9)
            ax.tick_params(axis='x', labelsize=9)
            ax.tick_params(axis='y', labelsize=8)
            sns.despine(ax=ax, top=True, right=True)

    fig.patch.set_facecolor('white')
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    for j, (header, _) in enumerate(COLUMNS):
        box = axes[0, j].get_position()
        fig.text((box.x0 + box.x1) / 2, 0.965, header, ha='center', va='center', fontsize=10)
    fig.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'{df["SP_ID"].nunique()} samples, {sum(len(t) for _, t in COLUMNS)} cell types')
    print(f'saved -> {OUT}')


if __name__ == '__main__':
    main()
