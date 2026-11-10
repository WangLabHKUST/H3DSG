# -*- coding: utf-8 -*-
"""Nearest-neighbour distance from T cell subsets to every cell type, PL against TL
(dot matrix).

Input:  ../data/min_distance_group_test.csv
        (written by pipeline/08_min_distance/02_group_test.py)
Output: min_distance_dotmatrix_tcells_to_all.png (written into this directory)

Reads ../data/ only; never touches the raw data.

How to read it
  Colour = the difference of medians, log2(TL / PL). The quantity is a distance, so
  red = TL further = PL closer = stronger interaction on the PL side; blue is the
  opposite.
  Dot radius is proportional to -log10(p); a black edge means p < 0.05 (no
  multiple-testing correction, see 02_group_test.py).

The variant is map_tum: the three malignant states are merged into a single Tumor before
the KD-tree is built. Do not instead take no_map and the minimum of its three columns to
make up a Tumor -- that is "the smallest of three per-state mean distances", which is
mathematically an upper bound of "the mean distance to the nearest tumour cell of any
kind", not the same thing.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from constants import celltype2_display, celltype2_order_map_tum
from panel_helpers import plot_dot_matrix

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'min_distance_dotmatrix_tcells_to_all.png'

VARIANT = 'map_tum'
ROWS = ['Tcell_CD4', 'Tcell_CD8', 'Treg']

# Layout: axes width 15 x 0.26 = 3.9 in, axes height 3 x 0.30 = 0.9 in.
# min_distance_dotmatrix_tumor_endo_to_immune.py uses different cell_w/cell_h to arrive at
# the same 3.9 x 0.9, which is what makes the colorbar and size legend the same size in
# both figures. Changing this means changing that one to match.
CELL_W, CELL_H = 0.26, 0.30


def main():
    src = DATA / 'min_distance_group_test.csv'
    if not src.exists():
        raise FileNotFoundError(
            f'{src} does not exist; run pipeline/08_min_distance/02_group_test.py first.')

    df = pd.read_csv(src)
    df = df[df['variant'] == VARIANT].copy()
    if df.empty:
        raise ValueError(f'{src} has no rows with variant == {VARIANT!r}')

    cols = [c for c in celltype2_order_map_tum
            if c in set(df['target'])]
    rows = [r for r in ROWS if r in set(df['source'])]
    missing = [r for r in ROWS if r not in rows]
    if missing:
        raise ValueError(f'these sources are not in the table: {missing}')

    # switch to display names (source/target and row_order/col_order must switch together)
    df['source'] = df['source'].map(lambda c: celltype2_display.get(c, c))
    df['target'] = df['target'].map(lambda c: celltype2_display.get(c, c))
    rows_disp = [celltype2_display.get(c, c) for c in rows]
    cols_disp = [celltype2_display.get(c, c) for c in cols]

    plot_dot_matrix(df, rows_disp, cols_disp, out_path=OUT,
                    cell_w=CELL_W, cell_h=CELL_H)
    sub = df[df['source'].isin(rows_disp) & df['target'].isin(cols_disp)]
    print(f'{len(rows)} x {len(cols)} pairs, {int((sub["mw_p"] < 0.05).sum())} with p < 0.05')


if __name__ == '__main__':
    main()
