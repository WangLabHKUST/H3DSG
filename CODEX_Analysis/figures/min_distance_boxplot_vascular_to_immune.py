# -*- coding: utf-8 -*-
"""Nearest-neighbour distance from vascular endothelium to the immune subsets, PL against
TL (boxplots with points).

Input:  ../data/min_distance_by_sample.csv
        ../data/min_distance_group_test.csv
        (written by the two scripts in pipeline/08_min_distance/)
Output: min_distance_boxplot_vascular_to_immune.png (written into this directory)

Reads ../data/ only; never touches the raw data.

One panel per (Vascular -> target) pair, one point per sample. The y axis is log2; for the
reason see pipeline/08_min_distance/02_group_test.py. The p value comes from the group
test table and is not corrected for multiple testing.

The matrix view of the same data is the second row of
min_distance_dotmatrix_tumor_endo_to_immune.py: that one shows every pair at a glance,
this one shows the distribution over samples.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from constants import celltype2_display, group_colors, group_labels, group_order
from panel_helpers import pair_box_grid

DATA = REPO / 'data'
OUT = Path(__file__).resolve().parent / 'min_distance_boxplot_vascular_to_immune.png'

VARIANT = 'map_tum'
SOURCE = 'Endothelial_cell'
TARGETS = ['Tcell_CD8', 'Tcell_CD4', 'Treg', 'NK', 'B_cell',
           'Myeloid_CD163', 'Myeloid_CD68', 'Myeloid_M2_like', 'Myeloid_other',
           'Neutrophil', 'Dendritic_cell']


def main():
    src = DATA / 'min_distance_by_sample.csv'
    test_src = DATA / 'min_distance_group_test.csv'
    for fp in (src, test_src):
        if not fp.exists():
            raise FileNotFoundError(
                f'{fp} does not exist; run the scripts in pipeline/08_min_distance/ first.')

    df = pd.read_csv(src)
    df = df[(df['variant'] == VARIANT) & (df['avg_min_distance'] > 0)].copy()
    df['log2_d'] = np.log2(df['avg_min_distance'])
    df['group_label'] = df['group'].map(group_labels)

    test = pd.read_csv(test_src)
    test = test[test['variant'] == VARIANT]

    targets = [t for t in TARGETS if t in set(df['target'])]
    missing = [t for t in TARGETS if t not in targets]
    if missing:
        raise ValueError(f'these targets are not in the table: {missing}')

    pairs = [(SOURCE, t, f'{celltype2_display.get(SOURCE, SOURCE)} to '
                         f'{celltype2_display.get(t, t)}') for t in targets]

    pair_box_grid(df, pairs, out_path=OUT, test_df=test, ncols=4,
                  group_order=group_order, group_colors=group_colors)
    print(f'{len(pairs)} pairs, {df["sample"].nunique()} samples')


if __name__ == '__main__':
    main()
