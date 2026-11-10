# -*- coding: utf-8 -*-
"""Cell type composition: cell counts and proportions of the 10 classes, per sample.

Input
  ../../raw_data/CODEX_adata_dict.npy      AnnData for 31 samples (archive version only)
  ../../sample_metadata.tsv                sample -> SP number, group

Output (committed; read by figures/celltype_proportion_stacked_bar.py)
  ../../data/celltype_counts.csv           cell count per sample per class, plus total
  ../../data/celltype_proportions.csv      as above, divided by total

Both tables are keyed on CODEX_ID and carry the SP_ID, group and group_as_run columns, so
that a plotting script never has to go back to the raw data and can label with the SP
number directly.

Differences from the original Proportion.py (the numbers are unchanged):
  * The original computed over `sample_lst = list(adata_dict.keys())` but plotted only
    g1+g2, so as soon as the two differed in length the csv and the figure were not the
    same set of samples. Here the analysis manifest in sample_metadata.tsv is the single
    reference, and the keys of adata_dict are asserted to match it.
  * The original mixed statistics and plotting in one file, so redrawing a figure meant
    loading the 7.9 GB npy again. The statistics now stay in pipeline/ and the plotting
    reads only the two CSVs above.
"""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import celltype_order, load_samples
from data_loading import ADATA_DICT, load_adata_dict

OUT_PATH = REPO / 'data'


def main():
    OUT_PATH.mkdir(parents=True, exist_ok=True)
    if not ADATA_DICT.exists():
        raise FileNotFoundError(
            f'{ADATA_DICT} does not exist. This step needs the raw data, which only the '
            'archive version carries; on the GitHub version, run figures/ off '
            '../../data/celltype_{counts,proportions}.csv.')

    print('loading CODEX_adata_dict.npy ...', flush=True)
    adata_dict = load_adata_dict()      # keyed by SP ID

    samples = load_samples()
    extra = sorted(set(adata_dict) - set(samples.index))
    short = sorted(set(samples.index) - set(adata_dict))
    if extra or short:
        raise ValueError(f'adata_dict does not match the analysis manifest in '
                         f'sample_metadata.tsv: extra {extra}, missing {short}')
    sample_lst = [s for s in samples.index if s in adata_dict]   # in metadata order

    counts = pd.DataFrame(index=pd.Index(sample_lst, name='CODEX_ID'), columns=celltype_order, dtype=int)
    for sample_name in sample_lst:
        vc = adata_dict[sample_name].obs['ct_agg'].value_counts()
        unknown = sorted(set(vc.index.astype(str)) - set(celltype_order))
        if unknown:
            raise ValueError(f'unknown ct_agg categories in {sample_name}: {unknown}')
        counts.loc[sample_name] = [int(vc.get(c, 0)) for c in celltype_order]

    counts['total'] = counts[celltype_order].sum(axis=1)
    prop = counts[celltype_order].div(counts['total'], axis=0).astype(float)

    for df in (counts, prop):
        for col in ('SP_ID', 'group', 'group_as_run'):
            df.insert(0, col, samples.loc[df.index, col])

    counts.to_csv(OUT_PATH / 'celltype_counts.csv')
    prop.to_csv(OUT_PATH / 'celltype_proportions.csv')

    print(f'{len(counts)} samples, {counts["total"].sum():,} cells in total')
    print(counts.groupby('group_as_run')['total'].agg(['size', 'sum']).to_string())
    print(f'\nsaved -> {OUT_PATH / "celltype_counts.csv"}')
    print(f'saved -> {OUT_PATH / "celltype_proportions.csv"}')


if __name__ == '__main__':
    main()
