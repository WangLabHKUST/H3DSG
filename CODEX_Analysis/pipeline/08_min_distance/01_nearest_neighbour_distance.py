# -*- coding: utf-8 -*-
"""Nearest-neighbour distance: one cell type x cell type matrix of mean nearest-neighbour
distances per sample.

Complementary to neighbourhood enrichment (04_interaction): enrichment asks "is there
one within 40 um", a binary notion of adjacency, whereas this asks "how far away is the
nearest one", a continuous quantity that can tell "neither is adjacent, but one is 50 um
away and the other 500 um" apart.

Input
  ../../raw_data/CODEX_adata_dict.npy     AnnData for 31 samples (archive version only)
  ../../raw_data/cell_type_2.parquet      the 17-class fine-grained annotation
  ../../sample_metadata.tsv

Output (committed; read by figures/min_distance_*.py)
  ../../data/min_distance_by_sample.csv
        sample, SP_ID, group, group_as_run, variant, source, target, avg_min_distance

Method
  For an ordered pair (A, B): take the Euclidean distance from each A cell to the nearest
  B cell, then average over all A cells. When A == B the cell itself is excluded and the
  second nearest is taken. Coordinates come from obsm['spatial'], in um. If either class
  has no cells in that sample the entry is NaN; so is the A == B entry when the class has
  only one cell.

The three variants
  no_map            the 17 classes as they are
  map_tum           the three malignant states merged into a single Tumor
                    (Cycling / OPC_OC_like / AC_MES_like)
  map_tum_map_mye   as above, and additionally the four myeloid subsets merged into a
                    single Myeloid (Myeloid_CD68 / CD163 / M2_like / other;
                    Dendritic_cell is not merged in)

  Note: merging relabels the cells before the KD-tree is built, so under map_tum
  "Tumor -> X" is the mean over every tumour cell of its distance to the nearest X of any
  kind, not the minimum of three per-state averages. The latter (what the development-era
  02/04 dot matrices did) is mathematically an upper bound of the former, not an
  equivalent. Only the former is produced here.

  Naming trap: in manuscript naming the malignant oligodendroglial/astrocytic states are
  OPC_OC_like / AC_MES_like, while OPC_OC means the normal oligodendrocyte precursor,
  which must **not** be merged into Tumor. The merge uses constants.malig_celltype;
  never write the class names out by hand.

Differences from the original codex/01_ccc_cal_min_dist.py (the numbers are unchanged)
  * The original rebuilt a KD-tree for every (A, B), i.e. 17x17 = 289 trees per sample.
    Here one tree is built per target class and reused across all sources, bringing the
    count down to 17.
  * The original ran one variant at a time, switched on the command line, with the
    outputs scattered over three sets of npy files. Here all three run in one pass into
    a single long table.
  * For A == B the original queried with k=2 and then dropped the self match row by row,
    with a Python-level for loop. Here k=2 is used as well, but the self match is removed
    in a vectorised way.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import load_samples, malig_celltype
from data_loading import load_adata_dict, load_cell_type_2, attach_cell_type_2

DATA = REPO / 'data'

MYELOID_SUBTYPES = ['Myeloid_CD68', 'Myeloid_CD163', 'Myeloid_M2_like', 'Myeloid_other']

VARIANTS = {
    'no_map': {},
    'map_tum': {ct: 'Tumor' for ct in malig_celltype},
    'map_tum_map_mye': {**{ct: 'Tumor' for ct in malig_celltype},
                        **{ct: 'Myeloid' for ct in MYELOID_SUBTYPES}},
}


def avg_min_distance_matrix(coords, labels, order):
    """Return DataFrame[order x order] of mean nearest-neighbour distance from the source
    in each row to the target in each column."""
    idx = {ct: np.flatnonzero(labels == ct) for ct in order}
    trees = {ct: (cKDTree(coords[i]) if len(i) else None) for ct, i in idx.items()}

    out = pd.DataFrame(np.nan, index=order, columns=order, dtype=float)
    for tgt in order:
        tree = trees[tgt]
        if tree is None:
            continue
        n_tgt = len(idx[tgt])
        for src in order:
            src_i = idx[src]
            if len(src_i) == 0:
                continue
            if src != tgt:
                d, _ = tree.query(coords[src_i], k=1, workers=-1)
                d = np.asarray(d, dtype=float)
            else:
                if n_tgt < 2:
                    continue            # only one cell, so there is no "other one of its kind"
                d, _ = tree.query(coords[src_i], k=2, workers=-1)
                d = np.asarray(d, dtype=float)[:, 1]   # column 0 is the cell itself, distance 0
            d = d[np.isfinite(d)]
            if d.size:
                out.loc[src, tgt] = float(d.mean())
    return out


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    samples = load_samples()
    print(f'{len(samples)} samples, variants: {list(VARIANTS)}')

    print('loading raw_data/CODEX_adata_dict.npy ...', flush=True)
    adata_dict = load_adata_dict(samples.index)
    ct2 = load_cell_type_2()

    frames = []
    for i, sample_name in enumerate(samples.index, 1):
        adata = adata_dict[sample_name]
        attach_cell_type_2(adata, ct2[sample_name])
        coords = np.asarray(adata.obsm['spatial'], dtype=np.float64)
        base = adata.obs['cell_type_2'].astype(str).to_numpy()
        base_order = list(adata.obs['cell_type_2'].cat.categories)

        for variant, mapping in VARIANTS.items():
            labels = np.array([mapping.get(c, c) for c in base])
            # keep the original class order; a merged class takes the position of its
            # first constituent
            order, seen = [], set()
            for c in base_order:
                c2 = mapping.get(c, c)
                if c2 not in seen:
                    seen.add(c2)
                    order.append(c2)

            mat = avg_min_distance_matrix(coords, labels, order)
            # melt rather than stack: newer pandas no longer accepts the dropna argument
            # of stack, whereas melt keeps every cell (NaN included) by definition and is
            # not version-sensitive
            long = mat.reset_index(names='source').melt(
                id_vars='source', var_name='target', value_name='avg_min_distance')
            long.insert(0, 'variant', variant)
            long.insert(0, 'sample', sample_name)
            frames.append(long)

        print(f'  [{i:2d}/{len(samples)}] {sample_name}: n={adata.n_obs:,}', flush=True)
        del adata.obs['cell_type_2']

    out = pd.concat(frames, ignore_index=True)
    for col in ('SP_ID', 'group', 'group_as_run'):
        out.insert(1, col, samples.loc[out['sample'], col].to_numpy())

    out.to_csv(DATA / 'min_distance_by_sample.csv', index=False)

    print(f'\n{len(out):,} rows')
    print(out.groupby('variant')['avg_min_distance'].agg(
        n='size', non_na=lambda s: int(s.notna().sum()),
        min='min', median='median', max='max').to_string())
    print(f'\nsaved -> {DATA / "min_distance_by_sample.csv"}')


if __name__ == '__main__':
    main()
