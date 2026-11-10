# -*- coding: utf-8 -*-
"""Stacked bars of the per-sample cell type composition, with a gap between PL and TL.

Input:  ../data/celltype_proportions.csv
        (written by pipeline/03_proportion/01_celltype_counts.py)
Output: celltype_proportion_stacked_bar.png (written into this directory)

Reads ../data/ only; never touches the raw data.

Samples are labelled with the SP number (as in Table S1 of the paper).

The ordering follows the original script: PL first, then TL, and within a group the order
of constants.g1_samples / g2_samples.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from constants import celltype_display, celltype_order, color_dict, g1_samples, g2_samples

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'celltype_proportion_stacked_bar.png'

LABEL_COLUMN = 'SP_ID'      # the column of the table that labels the bars
GAP = 0.5                   # gap between the two groups, in bar widths


def main():
    src = DATA / 'celltype_proportions.csv'
    if not src.exists():
        raise FileNotFoundError(
            f'{src} does not exist; run pipeline/03_proportion/01_celltype_counts.py first '
            '(it needs the raw data).')

    prop = pd.read_csv(src, index_col=0)

    ### PL first, then TL; within a group, the frozen order from constants
    sample_order = [s for s in g1_samples + g2_samples if s in prop.index]
    missing = sorted(set(prop.index) - set(sample_order))
    if missing:
        raise ValueError(f'{missing} is not in g1_samples + g2_samples, so it cannot be ordered')
    prop = prop.loc[sample_order]
    n_left = sum(s in g1_samples for s in sample_order)
    n_samples = len(sample_order)

    labels = (prop[LABEL_COLUMN].tolist() if LABEL_COLUMN in prop.columns
              else list(prop.index))

    x_left = np.arange(n_left)
    x_right = np.arange(n_left, n_samples) + GAP

    fig, ax = plt.subplots(figsize=(max(14, n_samples * 0.5), 8))
    bottom_left = np.zeros(n_left)
    bottom_right = np.zeros(n_samples - n_left)
    for cell_type in celltype_order:
        values = prop[cell_type].values
        ax.bar(x_left, values[:n_left], 0.7, bottom=bottom_left,
               color=color_dict[cell_type], alpha=0.7, edgecolor='grey', linewidth=1)
        ax.bar(x_right, values[n_left:], 0.7, bottom=bottom_right,
               color=color_dict[cell_type], alpha=0.7, edgecolor='grey', linewidth=1)
        bottom_left += values[:n_left]
        bottom_right += values[n_left:]

    legend_elements = [Patch(facecolor=color_dict[ct], edgecolor='black',
                             label=celltype_display[ct], alpha=0.7)
                       for ct in celltype_order]

    ax.set_xlabel('Sample', fontsize=12)
    ax.set_ylabel('Proportion', fontsize=12)
    ax.set_xticks(np.concatenate([x_left, x_right]))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_xlim(-0.5, x_right[-1] + 0.5)
    ax.set_ylim(0, 1.05)
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for side in ('bottom', 'left'):
        ax.spines[side].set_color('black')
        ax.spines[side].set_linewidth(1)

    ### label PL / TL above the midpoint of each group
    for xs, name in ((x_left, 'PL'), (x_right, 'TL')):
        ax.text(xs.mean(), 1.02, name, ha='center', va='bottom',
                fontsize=13, fontweight='bold')

    plt.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'{n_left} PL + {n_samples - n_left} TL = {n_samples} samples')
    print(f'saved -> {OUT}')


if __name__ == '__main__':
    main()
