# -*- coding: utf-8 -*-
"""Heatmap of pairwise cell type adjacency enrichment (pooled over 31 samples).

Input:  ../data/pairwise_enrichment_pooled.csv
        (written by pipeline/04_interaction/02_pairwise_enrichment.py)
Output: interaction_enrichment_heatmap.png (written into this directory)

Reads ../data/ only; never touches the raw data.

What is plotted is log2(observed/expected + 1), with cells below 0 masked out (left
blank), matching `res_matrix_all_log2_masked` of the original Interaction_inference.py:
blank means the adjacency of that cell type pair does **not** exceed random expectation,
and the heatmap is only there to show which pairs are enriched.
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
from constants import celltype_display, celltype_order

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'interaction_enrichment_heatmap.png'

VALUE_COLUMN = 'log2_ratio_plus1'      # log2(observed/expected + 1)


def main():
    src = DATA / 'pairwise_enrichment_pooled.csv'
    if not src.exists():
        raise FileNotFoundError(
            f'{src} does not exist; run pipeline/04_interaction/02_pairwise_enrichment.py '
            'first (it needs the raw data).')

    pooled = pd.read_csv(src)
    mat = pooled.pivot(index='cell_type1', columns='cell_type2', values=VALUE_COLUMN)
    mat = mat.reindex(index=celltype_order, columns=celltype_order)

    ### mask the negatives: show only the pairs whose adjacency exceeds expectation
    masked = mat.where(mat >= 0)

    display = [celltype_display[c] for c in celltype_order]
    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(masked, cmap='Reds', ax=ax,
                xticklabels=display, yticklabels=display,
                cbar_kws={'label': 'log2(observed / expected + 1)'})
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    plt.setp(ax.get_yticklabels(), rotation=0)
    plt.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches='tight')
    plt.close()

    n_masked = int(mat.lt(0).sum().sum())
    print(f'{mat.notna().sum().sum()} cell type pairs, {n_masked} of them below '
          f'expectation (masked)')
    print(f'saved -> {OUT}')


if __name__ == '__main__':
    main()
