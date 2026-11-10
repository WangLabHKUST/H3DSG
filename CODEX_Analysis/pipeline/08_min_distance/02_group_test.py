# -*- coding: utf-8 -*-
"""Group comparison of nearest-neighbour distances: test PL against TL for every
(source, target) pair.

Input
  ../../data/min_distance_by_sample.csv    output of 01_nearest_neighbour_distance.py
  ../../sample_metadata.tsv

Output (committed; read by figures/min_distance_*.py)
  ../../data/min_distance_group_test.csv
        variant, source, target, n_group1, n_group2,
        median_log2_group1, median_log2_group2, log2_fc, mw_p

Method
  log2 first: distances span two or three orders of magnitude across samples and cell
  types, so differences on the raw scale are not comparable. Non-positive values and
  NaNs (a class with no cells in that sample) are dropped outright.
  Two-sided Mann-Whitney U test with the sample as the unit of analysis.

  log2_fc = median(log2 d, TL) - median(log2 d, PL), i.e. the log2 ratio of the medians.
  Because the quantity is a distance, log2_fc > 0 means the TL side is further apart,
  that is, PL is closer. Figures colour this with RdBu_r, so red = PL closer = stronger
  interaction in PL.

  No multiple-testing correction is applied. The no_map variant alone has 17x17 = 289
  pairs, and filtering at p < 0.05 directly would yield many false positives. The table
  carries raw p values and figures/ only maps them to dot size; if a conclusion is to go
  into the text, correct them yourself (BH or similar) or restrict attention to
  pre-specified pairs.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import group_labels

DATA = REPO / 'data'
SRC = DATA / 'min_distance_by_sample.csv'
OUT = DATA / 'min_distance_group_test.csv'

G1, G2 = 'Group1', 'Group2'


def main():
    if not SRC.exists():
        raise FileNotFoundError(
            f'{SRC} does not exist; run '
            'pipeline/08_min_distance/01_nearest_neighbour_distance.py first')

    df = pd.read_csv(SRC)
    df = df[df['avg_min_distance'] > 0].copy()
    df['log2_d'] = np.log2(df['avg_min_distance'])

    rows = []
    for (variant, src, tgt), sub in df.groupby(['variant', 'source', 'target'], observed=True):
        a = sub.loc[sub['group'] == G1, 'log2_d'].dropna()
        b = sub.loc[sub['group'] == G2, 'log2_d'].dropna()
        if len(a) >= 1 and len(b) >= 1:
            p = float(stats.mannwhitneyu(a, b, alternative='two-sided').pvalue)
            fc = float(b.median() - a.median())
        else:
            p, fc = np.nan, np.nan
        rows.append({
            'variant': variant, 'source': src, 'target': tgt,
            'group1_label': group_labels[G1], 'group2_label': group_labels[G2],
            'n_group1': int(len(a)), 'n_group2': int(len(b)),
            'median_log2_group1': float(a.median()) if len(a) else np.nan,
            'median_log2_group2': float(b.median()) if len(b) else np.nan,
            'log2_fc': fc, 'mw_p': p,
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)

    print(f'{len(out)} pairs')
    # Named aggregation on a DataFrameGroupBy requires a (column, function) tuple for
    # every entry; n='size' as used on a SeriesGroupBy is not accepted here
    print(out.groupby('variant').agg(
        n=('mw_p', 'size'),
        testable=('mw_p', lambda s: int(s.notna().sum())),
        p_lt_05=('mw_p', lambda s: int((s < 0.05).sum()))).to_string())
    print(f'\nsaved -> {OUT}')


if __name__ == '__main__':
    main()
