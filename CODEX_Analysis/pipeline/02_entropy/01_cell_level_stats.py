# -*- coding: utf-8 -*-
"""Entropy step 1: compute the cell-level statistics the entropy analysis needs from
CODEX_adata_dict.npy.

Input
  ../../raw_data/CODEX_adata_dict.npy      AnnData for 31 samples (archive only, 7.9 GB)
  ../../sample_metadata.tsv                sample -> group
  ../../data_full/sampled_cell_ids.csv     the cell IDs of the existing down-sample
                                           (not resampled, see below)

Output (written to ../../data_full/; too large to commit)
  Full-size layer (7.09 M cells, parquet, so that a filter change later does not mean
  starting from the raw data again)
    1. entropy_cell_level.parquet     the cell-level main table
    2. neighbor_counts_r27_5.parquet  per-cell-type counts within the 27.5 um
                                      neighbourhood (index cell included)
  Downstream entry layer (1.32 M cells, csv, readable from both Python and R)
    3. entropy_cell_level_sampled.csv the cell-level main table after filtering and
                                      down-sampling  <= this is what downstream reads

The downstream computation and plotting scripts read only file 3 and never touch an
h5ad / adata_dict again.

--------------------------------------------------------------------------
Design rationale (all verified empirically; see
../../notes/Entropy_data_provenance.md):

* Spatial graph: the adata_dict holds only r50_adj, not r27_5, so the graph is rebuilt
  here from obsm['spatial']. The graph rebuilt with
  sklearn radius_neighbors_graph(radius=27.5) is edge for edge identical to the
  r27_5_connectivities in the original _modified.h5ad (diff_nnz=0). Coordinates are in
  um, and 27.5 is the diameter of a Visium spot.

* include_self=True corresponds to the original implementation's mask_diagonal=False
  (adj.setdiag(1)), i.e. the index cell counts towards its own neighbourhood, so
  total_num >= 1.

* Cell type naming: the ct_agg in the adata_dict has already been renamed and differs
  from the old scripts. old OPC_OC -> new OPC_OC_like; old AC_MES -> new AC_MES_like;
  old Oligo_pro -> new OPC_OC. This is a one-to-one rename and the malignant /
  non-malignant sets themselves are unchanged, so recomputing under the new naming is
  numerically equivalent to the old chain (measured entropy max|diff| = 2.2e-16).

* Down-sampling: the original script used .sample() without random_state and is not
  reproducible, but the cell IDs that were kept had been written to disk, and this
  script reuses data_full/sampled_cell_ids.csv directly. That list holds 1,320,620 IDs,
  all of which hit the obs_names of the adata_dict, so it can be used for a direct .loc.
  The file originally lived at
    20250929/20251210_dt/try_adjust_tumor/interacting_pairs_data/
    20260106_show_niche/all_samples_full_27_5/
    20260126_max_group1_index_filtered_non_malig_equal.csv
  and has been copied into data_full/; this script no longer reads that path.

* Load compatibility: CODEX_adata_dict.npy was pickled by an older version of anndata,
  so under anndata 0.10.x a plain np.load raises KeyError: '_adata_ref'.
  data_loading.load_adata_dict() applies the patch and returns the samples keyed by SP ID.

* Sample keys: data_full/ is keyed by the sample keys raw_data/ uses, not by SP IDs.
  Everything in this script works on SP IDs; data_loading.write_data_full() turns the
  sample column back into those keys when it writes, and read_data_full() turns them into
  SP IDs when it reads, so the files in data_full/ are the same whichever way they were
  made.

* Groups come from the group_as_run column of sample_metadata.tsv, i.e. the grouping that
  actually produced the existing outputs. The original implementation read
  common_data/samples.csv; the two are equivalent for these 31 samples.
--------------------------------------------------------------------------
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats
import scipy.sparse as sp
from sklearn.neighbors import radius_neighbors_graph
from tqdm import tqdm

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import (celltype_order, load_samples, malig_celltype,
                       non_malig_celltype, old2new_celltype)
from data_loading import ADATA_DICT, DATA_FULL, load_adata_dict, write_data_full

# ======================== paths and parameters ========================
OUT_PATH = DATA_FULL
SAMPLED_IDS = OUT_PATH / 'sampled_cell_ids.csv'

RADIUS = 27.5              # um, the diameter of a Visium spot
MALIG_CELLTYPE = malig_celltype
NON_MALIG_CELLTYPE = non_malig_celltype

# Filter thresholds, matching the original 20260129_entropy_cal_clean_v3_visium_spot_size.py
TUMOR_NUM_LIMIT = 2        # tumor_num >= 2
NON_MALIG_NUM_LIMIT = 2    # non_malig_num >= 2
TOTAL_NUM_LIMIT = 20       # total_num <= 20

# new naming -> old naming, only for comparing against historical results
NEW2OLD_CT = {v: k for k, v in old2new_celltype.items()}


def normalize_celltype(adata, k='ct_agg'):
    """Bring ct_agg onto the canonical manuscript naming.

    CODEX_adata_dict.npy is already in the new naming, so this usually does nothing; but
    if the input is switched to the old _modified.h5ad / cell_type_csv (which carry the
    development-era short names), this converts them, so that everything downstream has
    only one vocabulary to know about.

    Note that Series.replace cannot be applied item by item: in the old naming
    OPC_OC -> OPC_OC_like while Oligo_pro -> OPC_OC, so the wrong order would rewrite a
    value that has just been set. A single map() is used here.
    """
    labels = adata.obs[k].astype(str)
    present = set(labels.unique())

    is_old = bool(present & set(old2new_celltype))          # OPC_OC/AC_MES/Oligo_pro present
    # Oligo_pro exists only in the old naming and is the most reliable discriminator
    definitely_old = 'Oligo_pro' in present

    if definitely_old:
        labels = labels.map(lambda x: old2new_celltype.get(x, x))
    elif is_old and not definitely_old:
        # only OPC_OC / AC_MES and no Oligo_pro: possible under either naming, so leave
        # it to the caller to confirm
        unknown = present - set(celltype_order)
        if unknown:
            raise ValueError(
                f'the ct_agg naming cannot be determined; unknown categories: {sorted(unknown)}')

    unknown = set(labels.unique()) - set(celltype_order)
    if unknown:
        raise ValueError(f'ct_agg has categories outside celltype_order: {sorted(unknown)}')

    adata.obs[k] = labels
    return adata


# ======================== computation ========================
def entropy_calculation(df, cols=None):
    """Row-wise Shannon entropy. scipy returns NaN for an all-zero row."""
    if cols is None:
        cols = df.columns
    return scipy.stats.entropy(df[cols].values, axis=1)


def aggregate_neighbors_counts(adata, radius=RADIUS, k='ct_agg', celltype_order=celltype_order):
    """Count the cells of each type within radius of each cell (index cell included).

    Equivalent to
    aggregate_neighbors_counts(graph_key='r27_5_connectivities', mask_diagonal=False)
    in the original 20260129_process_all.py.
    """
    coords = np.asarray(adata.obsm['spatial'])
    adj = radius_neighbors_graph(coords, radius=radius, mode='connectivity', include_self=True)
    adj = (adj > 0).astype(np.float32)

    labels = adata.obs[k].astype('category')
    cats = labels.cat.categories
    onehot = sp.csr_matrix(
        (np.ones(len(labels), dtype=np.float32), (np.arange(len(labels)), labels.cat.codes.values)),
        shape=(len(labels), len(cats)))

    counts = (adj @ onehot).toarray()
    df = pd.DataFrame(counts, index=adata.obs_names, columns=cats)

    # fill in cell types absent from this sample, so all samples share the same columns
    for c in celltype_order:
        if c not in df.columns:
            df[c] = 0
    return df[celltype_order].astype(int)


def process_one_sample(adata, sample_name):
    """Return (neighbour count df, cell-level statistics df)."""
    adata = normalize_celltype(adata)
    nc = aggregate_neighbors_counts(adata)

    meta = pd.DataFrame(index=nc.index)
    meta.index.name = 'Object_ID'
    meta['sample'] = sample_name
    meta['ct_agg'] = adata.obs['ct_agg'].astype(str).values

    # carry over the columns the original script brought along, where present
    for col in ['r', 'l', 'r_cell', 'r50_adj_degree']:
        if col in adata.obs.columns:
            meta[col] = adata.obs[col].values

    meta['entropy_tumor'] = entropy_calculation(nc, cols=MALIG_CELLTYPE)
    meta['non_malig_entropy'] = entropy_calculation(nc, cols=NON_MALIG_CELLTYPE)
    meta['entropy_all'] = entropy_calculation(nc)

    meta['tumor_num'] = nc[MALIG_CELLTYPE].sum(axis=1)
    meta['non_malig_num'] = nc[NON_MALIG_CELLTYPE].sum(axis=1)
    meta['total_num'] = nc.sum(axis=1)

    # entropy between the two major classes, malignant and non-malignant
    major = pd.DataFrame({'malig': meta['tumor_num'], 'non_malig': meta['non_malig_num']})
    meta['entropy_major'] = entropy_calculation(major, cols=['malig', 'non_malig'])

    # the dominant malignant state in the neighbourhood; 'None' when there is no
    # malignant cell
    meta['dominated_tumor_state'] = nc[MALIG_CELLTYPE].idxmax(axis=1)
    meta.loc[meta['tumor_num'] == 0, 'dominated_tumor_state'] = 'None'

    return nc, meta


# ======================== main ========================
def main():
    OUT_PATH.mkdir(parents=True, exist_ok=True)
    if not ADATA_DICT.exists():
        raise FileNotFoundError(
            f'{ADATA_DICT} does not exist. This step needs the raw data, which only the '
            'archive version carries; on the GitHub version, run figures/ off the grid '
            'data already exported into ../../data/.')

    print('[1/4] loading CODEX_adata_dict.npy ...', flush=True)
    adata_dict = load_adata_dict()      # keyed by SP ID, in the order of the stored dict
    sample_lst = list(adata_dict.keys())
    print(f'      samples = {len(sample_lst)}', flush=True)

    sample_info = load_samples()
    unknown = sorted(set(sample_lst) - set(sample_info.index))
    if unknown:
        raise ValueError(
            f'samples in adata_dict are not in the analysis manifest of '
            f'sample_metadata.tsv: {unknown}')

    print(f'[2/4] rebuilding the {RADIUS}um neighbourhood graph per sample and computing '
          f'statistics ...', flush=True)
    nc_lst, meta_lst = [], []
    for sample_name in tqdm(sample_lst):
        nc, meta = process_one_sample(adata_dict[sample_name], sample_name)
        nc_lst.append(nc)
        meta_lst.append(meta)

    neighbor_counts_all = pd.concat(nc_lst, axis=0)
    meta_all = pd.concat(meta_lst, axis=0)
    meta_all['group'] = sample_info.loc[meta_all['sample'], 'group'].values

    assert meta_all.index.is_unique, 'duplicate Object_ID'
    assert (neighbor_counts_all.index == meta_all.index).all(), 'the two tables are not in the same row order'

    print(f'      total cells = {len(meta_all)}', flush=True)
    print(meta_all['group'].value_counts().to_string(), flush=True)

    print('[3/4] writing the full-size layer (parquet) ...', flush=True)
    # data_full/ is keyed by the original sample keys; write_data_full translates the
    # sample column back, so a rerun leaves these files as they were
    write_data_full(meta_all, 'entropy_cell_level.parquet')
    write_data_full(neighbor_counts_all, 'neighbor_counts_r27_5.parquet')

    print('[4/4] filtering and reusing the existing down-sampled cell IDs (no resampling) ...',
          flush=True)
    sel = pd.read_csv(SAMPLED_IDS, index_col=0)['index'].astype(str)

    missing = set(sel) - set(meta_all.index.astype(str))
    print(f'      sampled IDs = {len(sel)}, not found = {len(missing)}', flush=True)
    assert not missing, 'some cell IDs cannot be indexed in this run'

    # filter on the same conditions as the original script first, then check that every
    # sampled ID falls inside the filtered set
    mask = ((meta_all['tumor_num'] >= TUMOR_NUM_LIMIT) &
            (meta_all['non_malig_num'] >= NON_MALIG_NUM_LIMIT) &
            (meta_all['total_num'] <= TOTAL_NUM_LIMIT))
    meta_filtered = meta_all[mask]
    print(f'      passing the filter = {len(meta_filtered)}', flush=True)
    print(meta_filtered['group'].value_counts().to_string(), flush=True)

    not_in_filtered = set(sel) - set(meta_filtered.index.astype(str))
    if not_in_filtered:
        print(f'      [WARN] {len(not_in_filtered)} sampled IDs are not in the filtered '
              f'set, which means the filter conditions differ from the original run',
              flush=True)

    meta_sampled = meta_all.loc[sel.tolist()]
    write_data_full(meta_sampled, 'entropy_cell_level_sampled.csv')
    print(f'      rows in the downstream entry table = {len(meta_sampled)}', flush=True)
    print(meta_sampled['group'].value_counts().to_string(), flush=True)

    print(f'\ndone, output directory: {OUT_PATH}', flush=True)


if __name__ == '__main__':
    main()
