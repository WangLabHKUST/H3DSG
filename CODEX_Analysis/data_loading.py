# -*- coding: utf-8 -*-
"""Unified access layer for raw_data/ and data_full/.

Why this is a module of its own
  raw_data/CODEX_adata_dict.npy has three properties that every pipeline script would
  otherwise have to handle again: it was pickled by an old version of anndata, its obsp
  holds only r50_adj, and its obs holds only ct_agg. The steps used to carry one copy of
  the compatibility patch each, which did not scale, so the three concerns live here.
  constants.py holds constants only; the logic for reading data lives here.

What the dict has and what it is missing
  has      X (masked counts for 23 proteins), obsm['spatial'] (um), obs['ct_agg']
           (10 classes, already in manuscript naming),
           obsp['r50_adj_{connectivities,distances}']
  missing  the r40 / r27_5 / r50 / de graphs in obsp  -> rebuild from the coordinates
           with attach_radius_graph()
           cell_type_2 (17 classes) in obs            -> attach from the parquet with
           attach_cell_type_2()

  They are rebuilt rather than shipped because a graph is fully determined by
  obsm['spatial'], so there is no reason to spend another 1 GB on it. The rebuilt graph
  has been verified edge for edge against the r40_connectivities in
  {sample}_masked_counts_no_other_modified.h5ad (diff_nnz=0 on two samples);
  02_entropy/01 ran the same check for r27_5 with the same result.
  verify_radius_graph() keeps that check available.

Sample keys
  raw_data/ and data_full/ are keyed by the sample keys the data was first stored under,
  not by the SP IDs used everywhere else. The key layer below (SAMPLE_KEY, to_sp_id,
  to_raw_key, legacy_sample_order, read_data_full, write_data_full, ...) is the only
  place where the two meet: every loader here translates the keys to SP IDs on the
  way in, and every writer that stores into raw_data/ or data_full/ translates back, so
  code outside this module sees SP IDs only and those two directories stay as they were.
  The translation table is raw_data/sample_key.tsv (columns raw_key, SP_ID); it is
  written by pipeline/00_sample_metadata/01_build_sample_metadata.py and exists in the
  archive version only. Translating is idempotent: a value that is already an SP ID is
  left alone, so a table written with SP IDs loads the same way.
"""
import weakref
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.neighbors import radius_neighbors_graph

REPO = Path(__file__).resolve().parent
RAW_DATA = REPO / 'raw_data'
DATA_FULL = REPO / 'data_full'
ADATA_DICT = RAW_DATA / 'CODEX_adata_dict.npy'
CELL_TYPE_2 = RAW_DATA / 'cell_type_2.parquet'
TISSUE_BOUNDARY = RAW_DATA / 'tissue_boundary'
SAMPLE_KEY = RAW_DATA / 'sample_key.tsv'

_MISSING_RAW = (
    '{path} does not exist. This step needs the raw data, which only the archive '
    'version carries in raw_data/; on the GitHub version, run figures/ directly off '
    'the precomputed tables in data/.'
)
_MISSING_KEY = (
    '{path} does not exist. raw_data/ and data_full/ are keyed by the original sample '
    'keys, and this table (columns raw_key, SP_ID) is what turns them into SP IDs. It is '
    'written by pipeline/00_sample_metadata/01_build_sample_metadata.py and, like '
    'raw_data/ itself, only the archive version carries it.'
)


# ---------------------------------------------------------------------------
# Key layer: original sample key <-> SP ID
# ---------------------------------------------------------------------------
@lru_cache(maxsize=1)
def _key_pairs():
    """Return ((raw_key, SP_ID), ...) in the row order of raw_data/sample_key.tsv."""
    if not SAMPLE_KEY.exists():
        raise FileNotFoundError(_MISSING_KEY.format(path=SAMPLE_KEY))
    key = pd.read_csv(SAMPLE_KEY, sep='\t', dtype=str)
    if list(key.columns) != ['raw_key', 'SP_ID']:
        raise ValueError(f'{SAMPLE_KEY}: expected the columns raw_key, SP_ID, '
                         f'got {list(key.columns)}')
    if key.isna().any().any() or not key['raw_key'].is_unique or not key['SP_ID'].is_unique:
        raise ValueError(f'{SAMPLE_KEY}: raw_key and SP_ID must be unique and complete')
    return tuple(zip(key['raw_key'], key['SP_ID']))


def raw_to_sp():
    """{original key: SP ID}."""
    return dict(_key_pairs())


def sp_to_raw():
    """{SP ID: original key}."""
    return {sp_id: raw for raw, sp_id in _key_pairs()}


def to_raw_key(sp_id):
    """The key a sample has in raw_data/ and data_full/."""
    try:
        return sp_to_raw()[sp_id]
    except KeyError:
        raise KeyError(f'{sp_id!r} is not an SP ID in {SAMPLE_KEY}') from None


def to_sp_id(key):
    """SP ID of a sample given by its original key; an SP ID is returned unchanged."""
    lookup = raw_to_sp()
    if key in lookup:
        return lookup[key]
    if key in sp_to_raw():
        return key
    raise KeyError(f'{key!r} is neither an original key nor an SP ID in {SAMPLE_KEY}')


def legacy_sample_order(sp_ids):
    """Sort SP IDs the way the original keys sorted, i.e. the row order of sample_key.tsv.

    The tables in data/ were first written with the samples ordered by their original key
    (a groupby or sort on the sample column). SP IDs sort differently, so wherever a table
    is ordered by sample, use this to keep the rows where they were.
    """
    rank = {sp_id: i for i, (_, sp_id) in enumerate(_key_pairs())}
    return sorted(sp_ids, key=rank.__getitem__)


def _translate(values, lookup, what):
    out = values.map(lookup)
    bad = out.isna() & values.notna()
    if bad.any():
        raise ValueError(f'{int(bad.sum())} sample value(s) not in {SAMPLE_KEY} '
                         f'(e.g. {sorted(set(values[bad].astype(str)))[:3]}); {what}')
    return out


def translate_sample_column(df, col='sample'):
    """Return df with the original keys in column col replaced by SP IDs (copy)."""
    lookup = {**{sp_id: sp_id for sp_id in sp_to_raw()}, **raw_to_sp()}
    df = df.copy(deep=False)
    df[col] = _translate(df[col], lookup, 'cannot translate to SP IDs')
    return df


def to_raw_keys(values):
    """Series of SP IDs -> Series of original keys, for tables stored in raw_data/ or
    data_full/ (they are keyed by the original key). Original keys pass through."""
    lookup = {**{raw: raw for raw in raw_to_sp()}, **sp_to_raw()}
    return _translate(values, lookup, 'cannot translate to original keys')


def tissue_boundary_path(sp_id):
    """raw_data/tissue_boundary/<original key>.geojson for a sample given by SP ID."""
    return TISSUE_BOUNDARY / f'{to_raw_key(sp_id)}.geojson'


def read_data_full(name, **kwargs):
    """Read data_full/<name> (csv or parquet); a 'sample' column is returned as SP IDs."""
    path = DATA_FULL / name
    if not path.exists():
        raise FileNotFoundError(_MISSING_RAW.format(path=path))
    if path.suffix == '.csv':
        df = pd.read_csv(path, **kwargs)
    else:
        df = pd.read_parquet(path, **kwargs)
    if 'sample' in df.columns:
        df = translate_sample_column(df, 'sample')
    return df


def write_data_full(df, name):
    """Write df to data_full/<name> (csv or parquet) with the 'sample' column, if any, in
    the original keys, which is how data_full/ is keyed."""
    DATA_FULL.mkdir(parents=True, exist_ok=True)
    if 'sample' in df.columns:
        df = df.copy(deep=False)
        df['sample'] = to_raw_keys(df['sample'])
    path = DATA_FULL / name
    if path.suffix == '.csv':
        df.to_csv(path)
    else:
        df.to_parquet(path)
    return path


def patch_anndata_compat():
    """Make objects pickled by an old version of anndata loadable.

    The old AnnDataFileManager stored _adata (a strong reference) whereas 0.10.x expects
    _adata_ref (a weakref), and the old AnnData had no _n_obs / _n_vars cache. Without
    this patch, np.load raises KeyError: '_adata_ref' inside
    AnnDataFileManager.__setstate__.
    """
    import anndata
    from anndata._core.file_backing import AnnDataFileManager

    def _fm_setstate(self, state):
        self.__dict__ = dict(state)
        ref = state.get('_adata_ref', state.get('_adata', None))
        self.__dict__['_adata_ref'] = weakref.ref(ref) if ref is not None else (lambda: None)

    AnnDataFileManager.__setstate__ = _fm_setstate

    def _ad_setstate(self, state):
        self.__dict__.update(state)
        if '_n_obs' not in self.__dict__:
            obs = self.__dict__.get('_obs', None)
            self.__dict__['_n_obs'] = len(obs) if obs is not None else 0
        if '_n_vars' not in self.__dict__:
            var = self.__dict__.get('_var', None)
            self.__dict__['_n_vars'] = len(var) if var is not None else 0

    anndata.AnnData.__setstate__ = _ad_setstate


def load_adata_dict(samples=None):
    """Load raw_data/CODEX_adata_dict.npy and return {SP ID: AnnData}.

    The dict is stored under the original sample keys; they are translated to SP IDs here
    (the order of the dict is kept), and so is a 'sample' column in obs, if there is one.

    When samples is given (usually constants.load_samples().index), the keys of the dict
    are asserted to match it exactly and the result is returned in the order of samples,
    so that nothing downstream depends on the dict's own key order.
    """
    if not ADATA_DICT.exists():
        raise FileNotFoundError(_MISSING_RAW.format(path=ADATA_DICT))
    _key_pairs()                      # fail on a missing key table before the long load

    patch_anndata_compat()
    d = np.load(ADATA_DICT, allow_pickle=True).item()
    d = {to_sp_id(k): v for k, v in d.items()}
    for adata in d.values():
        if 'sample' in adata.obs.columns:
            adata.obs['sample'] = translate_sample_column(adata.obs, 'sample')['sample']

    if samples is not None:
        samples = list(samples)
        extra = sorted(set(d) - set(samples))
        short = sorted(set(samples) - set(d))
        if extra or short:
            raise ValueError(
                f'adata_dict does not match the sample manifest: extra {extra}, missing {short}')
        d = {s: d[s] for s in samples}
    return d


def attach_radius_graph(adata, radius, key=None, include_self=False):
    """Rebuild a radius neighbourhood graph from obsm['spatial'] and write it into obsp.

    radius is in um, the unit of obsm['spatial']. key defaults to 'r{radius}' following
    squidpy's convention, with the decimal point replaced by an underscore
    (27.5 -> 'r27_5'), so the keys written are '{key}_connectivities' and
    '{key}_distances'.

    include_self=False matches squidpy's spatial_neighbors(radius=...): no self loops.
    Pass True only where the index cell has to count towards its own neighbourhood, such
    as when counting the cells in a niche.
    """
    if key is None:
        key = 'r' + str(radius).replace('.', '_')
    coords = np.asarray(adata.obsm['spatial'], dtype=np.float64)

    conn = radius_neighbors_graph(coords, radius=radius, mode='connectivity',
                                  include_self=include_self)
    dist = radius_neighbors_graph(coords, radius=radius, mode='distance',
                                  include_self=include_self)
    adata.obsp[f'{key}_connectivities'] = conn.tocsr()
    adata.obsp[f'{key}_distances'] = dist.tocsr()
    return adata


def verify_radius_graph(adata, radius, h5ad_path, key=None):
    """Compare the rebuilt graph edge for edge against the one in {sample}_..._modified.h5ad.

    Usable only on the archive machine, where the h5ad files exist; this is the evidence
    that attach_radius_graph is correct. Returns (diff_nnz, ref_nnz, rebuilt_nnz);
    diff_nnz=0 means the two are identical.
    """
    import h5py
    if key is None:
        key = 'r' + str(radius).replace('.', '_')

    with h5py.File(h5ad_path, 'r') as f:
        g = f[f'obsp/{key}_connectivities']
        ref = sp.csr_matrix((g['data'][:], g['indices'][:], g['indptr'][:]),
                            shape=tuple(g.attrs['shape']))
    rebuilt = adata.obsp[f'{key}_connectivities']
    diff = ((ref > 0) != (rebuilt > 0)).nnz
    return diff, ref.nnz, rebuilt.nnz


def load_cell_type_2():
    """Read raw_data/cell_type_2.parquet and return {SP ID: Series} indexed by cell_id.

    This column is the 17-class fine-grained annotation, exported from the external
    adata_obs_dict.npy by pipeline/01_build_adata/04_export_cell_type_2.py and already
    converted to manuscript naming. The sample column is stored under the original keys
    and translated to SP IDs here.
    """
    if not CELL_TYPE_2.exists():
        raise FileNotFoundError(_MISSING_RAW.format(path=CELL_TYPE_2))
    _key_pairs()

    df = translate_sample_column(pd.read_parquet(CELL_TYPE_2), 'sample')
    return {s: g.set_index('cell_id')['cell_type_2'] for s, g in df.groupby('sample', observed=True)}


def attach_cell_type_2(adata, series):
    """Attach one sample's Series from load_cell_type_2() to adata.obs['cell_type_2'].

    Aligned explicitly on obs_names rather than relying on order; a mismatch raises.
    """
    s = series.reindex(adata.obs_names)
    if s.isna().any():
        n = int(s.isna().sum())
        raise ValueError(
            f'{n} cells are absent from the cell_type_2 table; obs_names differ from '
            'those at export time')
    adata.obs['cell_type_2'] = pd.Categorical(
        s.values, categories=list(series.cat.categories), ordered=True)
    return adata
