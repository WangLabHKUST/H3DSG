# -*- coding: utf-8 -*-
"""Ratio of T/NK cells to malignant cells, PL against TL (boxplot with points).

Input:  ../data/celltype_counts.csv
        (written by pipeline/03_proportion/01_celltype_counts.py)
Output: tcell_to_malignant_ratio_boxplot_PL_vs_TL.png (written into this directory)

Reads ../data/ only; never touches the raw data.

One point per sample (31 samples). ratio = cells of the T_NK class / cells of the three
malignant classes (constants.malig_celltype: Cycling, OPC_OC_like, AC_MES_like). T_NK is
the 10-class ct_agg label, so it includes NK cells; the CODEX cell type table behind it
(cell_type_2) separates CD4 T, CD8 T, Treg and NK, which this figure does not use.

The p value is a two-sided Mann-Whitney U test (normal approximation with continuity
correction, scipy's default for samples of this size), sample as the unit of analysis,
PL (n = 14) against TL (n = 17), no multiple-testing correction.

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
from constants import group_colors, group_labels, group_order, malig_celltype

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'tcell_to_malignant_ratio_boxplot_PL_vs_TL.png'

# seaborn's stripplot draws its jitter offsets from numpy's global RNG; seeding it makes
# this figure byte-reproducible. Same value as JITTER_SEED in panel_helpers.py.
JITTER_SEED = 42


def main():
    src = DATA / 'celltype_counts.csv'
    if not src.exists():
        raise FileNotFoundError(
            f'{src} does not exist; run pipeline/03_proportion/01_celltype_counts.py first.')

    df = pd.read_csv(src)
    df['ratio'] = df['T_NK'] / df[malig_celltype].sum(axis=1)
    df['group_label'] = pd.Categorical(df['group'].map(group_labels), categories=group_order, ordered=True)

    pl = df.loc[df['group_label'] == 'PL', 'ratio'].to_numpy()
    tl = df.loc[df['group_label'] == 'TL', 'ratio'].to_numpy()
    p = stats.mannwhitneyu(pl, tl, alternative='two-sided', method='asymptotic').pvalue

    np.random.seed(JITTER_SEED)        # see JITTER_SEED above
    fig, ax = plt.subplots(figsize=(2.3, 3.0))
    sns.boxplot(
        data=df, x='group_label', y='ratio', hue='group_label',
        order=group_order, hue_order=group_order, palette=group_colors,
        width=0.55, linewidth=1.0, showfliers=False, legend=False, ax=ax,
        medianprops={'color': 'black', 'linewidth': 1.6},
    )
    sns.stripplot(data=df, x='group_label', y='ratio', order=group_order,
                  color='black', size=2.6, jitter=0.18, ax=ax)
    ax.set_title(f'$P$ = {p:.1e}', fontsize=9)
    ax.set_xlabel('')
    ax.set_ylabel('T/NK cells / malignant cells', fontsize=9)
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
