# -*- coding: utf-8 -*-
"""Clustering coefficient before and after the density correction (violins, raw on the
left, adjusted on the right).

Input:  ../data/clustering_coefficient_adjusted.csv
        (written by pipeline/07_clustering_coeff/01_density_adjusted.py)
Output: clustering_coeff_violin_raw_vs_adj.png (written into this directory)

Reads ../data/ only; never touches the raw data.

Each violin is one cell type's distribution across the 31 samples, with the two groups
pooled -- this figure answers "did the correction change the relative ordering of the
cell types?". The group comparison is a different figure,
clustering_coeff_boxplot_PL_vs_TL.py.

The two panels cannot share a y tick step: raw runs over 0-0.7, while adj, divided by the
local cell density, is two orders of magnitude smaller (0-0.07). The steps are in
Y_TICK_STEP; set one to None to let matplotlib decide.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import seaborn as sns

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from constants import celltype_display, celltype_order, color_dict

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'clustering_coeff_violin_raw_vs_adj.png'

# key = the column that panel plots
Y_TICK_STEP = {
    'count': 0.2,
    'count_adj': 0.02,
}

PANELS = [
    ('count', 'Clustering coeff. (raw)'),
    ('count_adj', 'Clustering coeff. (density-adjusted)'),
]


def main():
    src = DATA / 'clustering_coefficient_adjusted.csv'
    if not src.exists():
        raise FileNotFoundError(
            f'{src} does not exist; run '
            'pipeline/07_clustering_coeff/01_density_adjusted.py first.')

    df = pd.read_csv(src)
    types = [c for c in celltype_order if c in set(df['cell_type'])]
    pal = {c: color_dict[c] for c in types}

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 3.8), sharey=False)
    for ax, (ycol, ylabel) in zip(axes, PANELS):
        sns.violinplot(
            data=df, x='cell_type', y=ycol, hue='cell_type',
            order=types, hue_order=types, palette=pal,
            dodge=False, legend=False, inner='box',
            linewidth=1.2, cut=0, density_norm='width', ax=ax,
        )
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_xlabel('')
        ax.set_xticks(range(len(types)))
        ax.set_xticklabels([celltype_display.get(t, t) for t in types],
                           rotation=45, ha='right')
        step = Y_TICK_STEP.get(ycol)
        if step:
            ax.yaxis.set_major_locator(mticker.MultipleLocator(step))
        sns.despine(ax=ax, top=True, right=True)

    fig.patch.set_facecolor('white')
    plt.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f'{df["sample"].nunique()} samples x {len(types)} cell types')
    print(f'saved -> {OUT}')


if __name__ == '__main__':
    main()
