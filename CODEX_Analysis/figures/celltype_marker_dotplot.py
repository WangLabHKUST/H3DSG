# -*- coding: utf-8 -*-
"""Dot plot of every cell type across the protein channels.

Input:  ../data/marker_summary_by_celltype.csv
        (written by pipeline/05_marker_expression/01_marker_summary_by_celltype.py)
Output: celltype_marker_dotplot.png (written into this directory)

Reads ../data/ only; never touches the raw data.

Dot colour = that type's mean expression on that channel, **linearly stretched to [0, 1]
within each channel** (i.e. standard_scale='var' of scanpy's dotplot, as in the original
Show_dotplot.py); dot size = the fraction of cells of that type positive (> 0) on that
channel.

The original script called `sc.pl.dotplot(adata_full, ...)` directly, which needs the 31
samples concatenated into one full-size AnnData; this reads the 2 KB summary table and
draws the plot itself, with the same visual meaning.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from constants import celltype_display, celltype_order

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'celltype_marker_dotplot.png'

MAX_DOT_AREA = 180      # dot area at frac_expressing == 1


def main():
    src = DATA / 'marker_summary_by_celltype.csv'
    if not src.exists():
        raise FileNotFoundError(
            f'{src} does not exist; run '
            'pipeline/05_marker_expression/01_marker_summary_by_celltype.py first '
            '(it needs the raw data).')

    df = pd.read_csv(src)
    markers = list(dict.fromkeys(df['marker']))          # keep the channel order of the file
    cts = [c for c in celltype_order if c in set(df['cell_type'])]

    mean = df.pivot(index='cell_type', columns='marker', values='mean_expression') \
             .reindex(index=cts, columns=markers)
    frac = df.pivot(index='cell_type', columns='marker', values='frac_expressing') \
             .reindex(index=cts, columns=markers)

    ### standard_scale='var': stretch each channel to [0, 1] separately; a constant
    ### column becomes 0
    rng = mean.max(axis=0) - mean.min(axis=0)
    scaled = (mean - mean.min(axis=0)).div(rng.replace(0, np.nan), axis=1).fillna(0.0)

    xx, yy = np.meshgrid(np.arange(len(markers)), np.arange(len(cts)))
    fig, ax = plt.subplots(figsize=(0.46 * len(markers) + 4, 0.46 * len(cts) + 2.2))
    sc = ax.scatter(xx.ravel(), yy.ravel(),
                    s=frac.values.ravel() * MAX_DOT_AREA,
                    c=scaled.values.ravel(), cmap='Reds', vmin=0, vmax=1,
                    edgecolor='grey', linewidth=0.3)

    ax.set_xticks(np.arange(len(markers)))
    ax.set_xticklabels(markers, rotation=90)
    ax.set_yticks(np.arange(len(cts)))
    ax.set_yticklabels([celltype_display[c] for c in cts])
    ax.set_xlim(-0.6, len(markers) - 0.4)
    ax.set_ylim(-0.6, len(cts) - 0.4)
    ax.invert_yaxis()
    ax.grid(False)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)

    cbar = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label('Mean expression\n(scaled per marker)', fontsize=9)

    ### dot size legend
    handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor='grey',
                      markeredgecolor='grey', markersize=np.sqrt(f * MAX_DOT_AREA),
                      label=f'{int(f * 100)}%')
               for f in (0.25, 0.5, 0.75, 1.0)]
    ax.legend(handles=handles, title='Fraction positive', loc='upper left',
              bbox_to_anchor=(1.18, 1.0), fontsize=8, title_fontsize=9, frameon=False)

    plt.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'{len(cts)} cell types x {len(markers)} channels')
    print(f'saved -> {OUT}')


if __name__ == '__main__':
    main()
