# -*- coding: utf-8 -*-
"""Export cell_type_2 (the 17-class fine-grained annotation) for the downstream
nearest-neighbour distance analysis.

Background
  The obs of raw_data/CODEX_adata_dict.npy carries only ct_agg (10 classes), not
  cell_type_2. cell_type_2 is produced by 00_update_adata_obs.py in this directory and
  lands in the external h5ad/adata_obs_dict.npy. This script takes that column for the
  31 analysis samples, converts it to the manuscript naming, and stores it under
  raw_data/ so that nothing downstream has to touch an external path.

  Like the other three scripts in this directory, this is an **archival record**: it
  reads external files under h5ad/ and can only run on the original analysis machine. A
  collaborator receives its output, raw_data/cell_type_2.parquet.
  (00_create_adata.py cannot even run there any more -- the {sample}_full_raw.h5ad files
  it depends on no longer exist; see section 1.0 of
  notes/Entropy_data_provenance.md.)

Input (not part of the repository; expected under raw_data/)
  ../../raw_data/h5ad/adata_obs_dict.npy
                                   obs slices for 33 samples;
                                   columns = ct_agg / cell_type / cell_type_2
                                   (410 MB, so not shipped with the archive version;
                                   place or link it there to run this script)

Output (archive version only; raw_data/ is in .gitignore)
  ../../raw_data/cell_type_2.parquet
        one row per cell: sample, cell_id, cell_type_2
        cell_type_2 already converted to manuscript naming via
        constants.old2new_celltype

Cell ordering
  Verified: the index of adata_obs_dict is identical to the obs_names of
  CODEX_adata_dict for all 31 samples, one by one (the cell counts match too, 7,092,596
  in total). The cell_id column is kept anyway so that downstream code can join
  explicitly instead of relying on order.

Sample keys
  adata_obs_dict.npy and the parquet written here are keyed by the original sample keys
  (raw_data/ is). The script works on SP IDs and translates with data_loading.to_raw_key()
  for the lookup and for the sample column it writes; data_loading.load_cell_type_2()
  translates back to SP IDs on reading.

Naming
  adata_obs_dict stores the development-era naming, whereas the ct_agg in the dict has
  already been renamed -- the two do not use the same vocabulary. This script applies
  constants.old2new_celltype uniformly:

      OPC_OC    -> OPC_OC_like     (malignant state)
      AC_MES    -> AC_MES_like     (malignant state)
      Oligo_pro -> OPC_OC          (normal oligodendrocyte precursor)

  Note this is a three-way exchange, not a simple rename: the old OPC_OC and the new
  OPC_OC are not the same population. It has to be a single atomic .map(); a chained
  .replace() would merge Oligo_pro into OPC_OC_like by mistake.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import old2new_celltype, load_samples
from data_loading import ADATA_DICT, CELL_TYPE_2, load_adata_dict, to_raw_key

H5AD_PATH = REPO / 'raw_data' / 'h5ad'
OBS_DICT = H5AD_PATH / 'adata_obs_dict.npy'
OUT_PATH = CELL_TYPE_2

# Display order of the 17 classes: malignant states -> myeloid -> granulocytes ->
# lymphoid -> stromal/normal
CELLTYPE2_ORDER = [
    'Cycling', 'OPC_OC_like', 'AC_MES_like',
    'Dendritic_cell', 'Myeloid_CD68', 'Myeloid_CD163', 'Myeloid_M2_like', 'Myeloid_other',
    'Neutrophil',
    'Tcell_CD4', 'Tcell_CD8', 'Treg', 'NK', 'B_cell',
    'Endothelial_cell', 'Astrocyte', 'OPC_OC',
]


def main():
    if not OBS_DICT.exists():
        raise FileNotFoundError(
            f'{OBS_DICT} does not exist. This step is an archival record and can only '
            'run on the original analysis machine; collaborators should use the '
            'already generated raw_data/cell_type_2.parquet.')

    samples = load_samples()
    print(f'{len(samples)} analysis samples')

    print('loading adata_obs_dict.npy ...', flush=True)
    obs_dict = np.load(OBS_DICT, allow_pickle=True).item()

    missing = sorted(s for s in samples.index if to_raw_key(s) not in obs_dict)
    if missing:
        raise ValueError(f'adata_obs_dict is missing analysis samples: {missing}')

    # Reconcile against the obs_names of the dict (only the archive machine has the
    # dict; if it is absent, skip the check but say so)
    if ADATA_DICT.exists():
        print('loading CODEX_adata_dict.npy to reconcile the cell order ...', flush=True)
        adata_dict = load_adata_dict()
        for s in samples.index:
            a = np.asarray(adata_dict[s].obs_names)
            c = np.asarray(obs_dict[to_raw_key(s)].index)
            if len(a) != len(c) or not (a == c).all():
                raise ValueError(
                    f'{s}: the cell order of adata_dict and adata_obs_dict differs')
        print('  the cell order agrees for all 31 samples')
        del adata_dict
    else:
        print(f'  WARNING: {ADATA_DICT} is absent, skipping the cell order check')

    frames = []
    for s in samples.index:
        obs = obs_dict[to_raw_key(s)]
        ct2 = obs['cell_type_2'].astype(str)

        unknown = sorted(set(ct2.unique()) - set(old2new_celltype) - set(CELLTYPE2_ORDER))
        if unknown:
            raise ValueError(f'{s}: unknown cell_type_2 categories {unknown}')

        # single atomic mapping: a three-way exchange cannot be done by chained replace
        ct2 = ct2.map(lambda c: old2new_celltype.get(c, c))

        frames.append(pd.DataFrame({
            'sample': to_raw_key(s),
            'cell_id': np.asarray(obs.index),
            'cell_type_2': ct2.values,
        }))

    out = pd.concat(frames, ignore_index=True)
    out['cell_type_2'] = pd.Categorical(out['cell_type_2'], categories=CELLTYPE2_ORDER, ordered=True)
    if out['cell_type_2'].isna().any():
        bad = sorted(set(pd.concat(frames)['cell_type_2']) - set(CELLTYPE2_ORDER))
        raise ValueError(f'categories still outside CELLTYPE2_ORDER after renaming: {bad}')

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUT_PATH, index=False)

    print(f'\n{len(out):,} cells, {out["sample"].nunique()} samples')
    print(out['cell_type_2'].value_counts().reindex(CELLTYPE2_ORDER).to_string())
    print(f'\nsaved -> {OUT_PATH}  ({OUT_PATH.stat().st_size / 1e6:.1f} MB)')


if __name__ == '__main__':
    main()
