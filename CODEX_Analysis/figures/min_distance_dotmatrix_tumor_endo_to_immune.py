# -*- coding: utf-8 -*-
"""Nearest-neighbour distance from tumour cells and vascular endothelium to the immune
subsets, PL against TL (dot matrix).

Input:  ../data/min_distance_group_test.csv
        (written by pipeline/08_min_distance/02_group_test.py)
Output: min_distance_dotmatrix_tumor_endo_to_immune.png (written into this directory)

Reads ../data/ only; never touches the raw data.

Read it the same way as min_distance_dotmatrix_tcells_to_all.py: red = PL closer, dot
radius proportional to -log10(p), black edge = p < 0.05.

This figure answers "from these two starting points, tumour and vasculature, which immune
cells is each closer to?". It uses map_tum (the myeloid subsets stay unfolded), so the
four myeloid populations CD68+ / CD163+ / double positive / double negative remain
distinguishable.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from constants import celltype2_display
from panel_helpers import plot_dot_matrix

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'min_distance_dotmatrix_tumor_endo_to_immune.png'

VARIANT = 'map_tum'
ROWS = ['Tumor', 'Endothelial_cell']
COLS = ['Tcell_CD8', 'Tcell_CD4', 'Treg', 'NK', 'B_cell',
        'Myeloid_CD163', 'Myeloid_CD68', 'Myeloid_M2_like', 'Myeloid_other',
        'Neutrophil', 'Dendritic_cell']

# Layout: axes width 11 x (3.9/11) = 3.9 in, axes height 2 x 0.45 = 0.9 in, exactly the
# 3.9 x 0.9 of min_distance_dotmatrix_tcells_to_all.py. That is what makes the colorbar
# and size legend the same size in both, so side by side neither looks larger.
#
# Careful: with only 2 rows, cell_h cannot be pushed lower. The colorbar and the size
# legend are anchored in ax.transAxes units, so a short axes collapses the sidebar into a
# clump with overlapping labels. 0.45 is the measured lower bound.
CELL_W, CELL_H = 3.9 / 11, 0.45
# The x labels are much longer than in that figure (CD68+CD163+ Myeloid), so reserve more
# height
FIG_H_PAD = 1.60


def main():
    src = DATA / 'min_distance_group_test.csv'
    if not src.exists():
        raise FileNotFoundError(
            f'{src} does not exist; run pipeline/08_min_distance/02_group_test.py first.')

    df = pd.read_csv(src)
    df = df[df['variant'] == VARIANT].copy()
    if df.empty:
        raise ValueError(f'{src} has no rows with variant == {VARIANT!r}')

    rows = [r for r in ROWS if r in set(df['source'])]
    cols = [c for c in COLS if c in set(df['target'])]
    missing = [c for c in ROWS + COLS if c not in set(df['source']) | set(df['target'])]
    if missing:
        raise ValueError(f'these cell types are not in the table: {missing}')

    df['source'] = df['source'].map(lambda c: celltype2_display.get(c, c))
    df['target'] = df['target'].map(lambda c: celltype2_display.get(c, c))
    rows_disp = [celltype2_display.get(c, c) for c in rows]
    cols_disp = [celltype2_display.get(c, c) for c in cols]

    plot_dot_matrix(df, rows_disp, cols_disp, out_path=OUT,
                    cell_w=CELL_W, cell_h=CELL_H, fig_h_pad=FIG_H_PAD)
    sub = df[df['source'].isin(rows_disp) & df['target'].isin(cols_disp)]
    print(f'{len(rows)} x {len(cols)} pairs, {int((sub["mw_p"] < 0.05).sum())} with p < 0.05')


if __name__ == '__main__':
    main()
