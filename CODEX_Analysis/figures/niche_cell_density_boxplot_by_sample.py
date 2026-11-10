# -*- coding: utf-8 -*-
"""Local cell density of every sample (one box per sample, 100 niches each).

Input:  ../data/niche_meta.csv
        (written by pipeline/06_niche_density/01_sample_niches.py)
Output: niche_cell_density_boxplot_by_sample.png (written into this directory)

Reads ../data/ only; never touches the raw data.

Each sample contributes its 100 niches (a niche is every cell within 40 um of an index
cell, index cell included; num_cells in niche_meta.csv). The cell count is divided by the
area of the 40 um disc and drawn as cells per mm^2 (x10^3). Samples are ordered by the
median density of their niches; the box is the quartiles of the 100 niches, the whiskers
reach the most extreme niche within 1.5 x IQR, and every niche is drawn as a small point.
Cell counts are integers, so the points fall on a discrete ladder of densities.

niche_meta.csv comes from pipeline/06, which uses the niche centres sampled during the
analysis (data/niche_centers_analysis_2026-06-14.csv; see the notes in that script).

Boxplots never show outliers (showfliers=False): every box has points drawn on top, so if
the boxplot also drew its outliers the same niche would appear twice.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from constants import group_colors, group_labels, group_order

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'niche_cell_density_boxplot_by_sample.png'

RADIUS_MM = 0.040                       # niche radius, 40 um (pipeline/06 RADIUS)
DISC_AREA_MM2 = np.pi * RADIUS_MM ** 2

# seaborn's stripplot draws its jitter offsets from numpy's global RNG; seeding it makes
# this figure byte-reproducible. Same value as JITTER_SEED in panel_helpers.py.
JITTER_SEED = 42


def main():
    src = DATA / 'niche_meta.csv'
    if not src.exists():
        raise FileNotFoundError(
            f'{src} does not exist; run pipeline/06_niche_density/01_sample_niches.py first.')

    df = pd.read_csv(src)
    df['density_k'] = df['num_cells'] / DISC_AREA_MM2 / 1000        # x10^3 cells / mm^2
    df['group_label'] = df['group'].map(group_labels)

    summary = (df.groupby('SP_ID')
               .agg(group_label=('group_label', 'first'), median=('density_k', 'median'))
               .reset_index()
               .sort_values(['median', 'SP_ID']))
    order = summary['SP_ID'].tolist()
    palette = {sp: group_colors[g] for sp, g in zip(summary['SP_ID'], summary['group_label'])}

    np.random.seed(JITTER_SEED)        # see JITTER_SEED above
    fig, ax = plt.subplots(figsize=(9.0, 3.2))
    sns.boxplot(
        data=df, x='SP_ID', y='density_k', hue='SP_ID', order=order, hue_order=order,
        palette=palette, width=0.8, linewidth=0.8, showfliers=False, legend=False, ax=ax,
        medianprops={'color': 'black', 'linewidth': 1.4},
    )
    sns.stripplot(data=df, x='SP_ID', y='density_k', order=order, color='black',
                  size=0.8, jitter=0.35, alpha=0.6, ax=ax)
    ax.set_xlabel('Samples', fontsize=9)
    ax.set_ylabel('Cell density (x10$^3$ cells/mm$^2$)', fontsize=9)
    ax.tick_params(axis='x', labelsize=7, rotation=45)
    ax.tick_params(axis='y', labelsize=8)
    ax.set_ylim(bottom=0)
    ax.legend(handles=[Patch(facecolor=group_colors[g], edgecolor='black', linewidth=0.5,
                             label=f'H3-DSG-{g}') for g in group_order],
              loc='upper left', frameon=False, fontsize=8)
    sns.despine(ax=ax, top=True, right=True)

    fig.patch.set_facecolor('white')
    plt.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'{len(order)} samples x {df.groupby("SP_ID").size().iloc[0]} niches')
    print(f'saved -> {OUT}')


if __name__ == '__main__':
    main()
