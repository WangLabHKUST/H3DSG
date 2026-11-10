# -*- coding: utf-8 -*-
"""Density-adjusted clustering coefficient: divide same-type self-aggregation by the
sample's local cell density.

Why the correction is needed
  The clustering coefficient is the diagonal of the neighbourhood enrichment matrix
  (cell_type1 == cell_type2) and measures how strongly a cell type sits next to its own
  kind. It is already normalised once by degree (see 04_interaction/01), but cell density
  differs fourfold between samples (density_factor runs from 10 to 40), and in a dense
  sample any two cells are more likely to fall within 40 um of each other, which lifts
  the whole diagonal. Dividing by each sample's local cell density removes that term.

Input
  ../../data/neighborhood_enrichment_by_sample.csv   output of 04_interaction/01
  ../../data/sample_density_factor.csv               output of 06_niche_density/01
  ../../sample_metadata.tsv

Output (committed; read by figures/clustering_coeff_*.py)
  ../../data/clustering_coefficient_adjusted.csv
        one row per sample per cell type: count (raw), density_factor,
        count_adj (corrected)
  ../../data/clustering_coefficient_adjusted_group_test.csv
        one row per cell type: n and median for each group and the p value of a
        two-sided Mann-Whitney U test, raw and adjusted side by side so that it is easy
        to see whether the correction changes the conclusion

Method
  count_adj = count / density_factor

  density_factor is the mean cell count over that sample's 100 niches (a niche being a
  40 um neighbourhood, index cell included), computed by 06_niche_density/01.

  The group comparison is a two-sided Mann-Whitney U test on log2(count_adj), with the
  sample as the unit of analysis. The log2 is taken because count_adj spans three orders
  of magnitude across samples and cell types, so differences on the raw scale are not
  comparable.

Differences from the original
202604_rebuttal_related/20260506_cluster_coeff_num_cells_adjust_violin.py
  * The original mixed statistics and plotting in one file; only the statistics are kept
    here, and the plotting lives in figures/.
  * The original read the development-era fig3e_*_interaction_df_ct_agg_40_prop.csv
    directly; here the same table is read as recomputed by 04_interaction/01, so a script
    inside this repository is responsible for how it is defined.
  * Where density_factor was missing or zero, the original only warned and set count_adj
    to NaN; here it raises -- all 31 samples should have a density_factor, and a missing
    one means the previous step did not finish.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import celltype_order, load_samples, group_labels

DATA = REPO / 'data'
ENRICHMENT = DATA / 'neighborhood_enrichment_by_sample.csv'
DENSITY = DATA / 'sample_density_factor.csv'


def main():
    for fp, who in ((ENRICHMENT, 'pipeline/04_interaction/01_neighborhood_enrichment.py'),
                    (DENSITY, 'pipeline/06_niche_density/01_sample_niches.py')):
        if not fp.exists():
            raise FileNotFoundError(f'{fp} does not exist; run {who} first')

    samples = load_samples()
    enr = pd.read_csv(ENRICHMENT)
    dens = pd.read_csv(DENSITY)

    # clustering coefficient = the diagonal of the neighbourhood enrichment matrix
    cc = enr[enr['cell_type1'] == enr['cell_type2']].copy()
    cc = cc.rename(columns={'cell_type1': 'cell_type'}).drop(columns=['cell_type2'])
    cc = cc[['sample', 'cell_type', 'count', 'zscore']]

    missing = sorted(set(cc['sample']) - set(dens['sample']))
    if missing:
        raise ValueError(f'these samples have no density_factor: {missing}')

    cc = cc.merge(dens[['sample', 'density_factor']], on='sample', how='left')
    bad = cc['density_factor'].isna() | (cc['density_factor'] <= 0)
    if bad.any():
        raise ValueError(f'{int(bad.sum())} rows have a missing or non-positive density_factor')

    cc['count_adj'] = cc['count'] / cc['density_factor']
    for col in ('SP_ID', 'group', 'group_as_run'):
        cc.insert(1, col, samples.loc[cc['sample'], col].to_numpy())

    # Within a cell type and group the samples keep the order they have in
    # sample_density_factor.csv (the order of the original sample keys, see step 06), not
    # the alphabetical order of the SP IDs, so that the rows stay where the committed
    # table has them.
    cc['cell_type'] = pd.Categorical(cc['cell_type'], categories=celltype_order, ordered=True)
    cc['_sample_rank'] = cc['sample'].map({s: i for i, s in enumerate(dens['sample'])})
    cc = (cc.sort_values(['cell_type', 'group', '_sample_rank'])
            .drop(columns='_sample_rank').reset_index(drop=True))
    cc.to_csv(DATA / 'clustering_coefficient_adjusted.csv', index=False)

    # ---- group test: raw and adjusted side by side ----
    g1, g2 = 'Group1', 'Group2'
    rows = []
    for ct, sub in cc.groupby('cell_type', observed=True):
        rec = {'cell_type': ct,
               'group1_label': group_labels[g1], 'group2_label': group_labels[g2]}
        for metric, col in (('raw', 'count'), ('adjusted', 'count_adj')):
            # load_samples() maps the group column to Group1/Group2 while group_as_run
            # stays PL/TL, so filtering by group must use the group column
            v = np.log2(sub[col].where(sub[col] > 0))
            a = v[sub['group'] == g1].dropna()
            b = v[sub['group'] == g2].dropna()
            if len(a) >= 1 and len(b) >= 1:
                p = float(stats.mannwhitneyu(a, b, alternative='two-sided').pvalue)
            else:
                p = np.nan
            rec[f'n_group1_{metric}'] = int(len(a))
            rec[f'n_group2_{metric}'] = int(len(b))
            rec[f'median_log2_group1_{metric}'] = float(a.median()) if len(a) else np.nan
            rec[f'median_log2_group2_{metric}'] = float(b.median()) if len(b) else np.nan
            rec[f'mw_p_{metric}'] = p
        rows.append(rec)

    test = pd.DataFrame(rows)
    test.to_csv(DATA / 'clustering_coefficient_adjusted_group_test.csv', index=False)

    print(f'{cc["sample"].nunique()} samples x {cc["cell_type"].nunique()} cell types'
          f' = {len(cc)} rows')
    print(f'density_factor ranges {dens["density_factor"].min():.1f} - '
          f'{dens["density_factor"].max():.1f}')
    print('\nGroup test (rows with p < 0.05):')
    sig = test[(test['mw_p_raw'] < 0.05) | (test['mw_p_adjusted'] < 0.05)]
    cols = ['cell_type', 'mw_p_raw', 'mw_p_adjusted']
    print(sig[cols].to_string(index=False) if len(sig) else '  none')
    print(f'\nsaved -> {DATA / "clustering_coefficient_adjusted.csv"}')
    print(f'saved -> {DATA / "clustering_coefficient_adjusted_group_test.csv"}')


if __name__ == '__main__':
    main()
