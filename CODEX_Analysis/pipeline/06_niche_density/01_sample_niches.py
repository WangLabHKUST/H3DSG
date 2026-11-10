# -*- coding: utf-8 -*-
"""Niches and local cell density: 100 non-adjacent niches per sample, the cells in each
neighbourhood counted.

The density_factor this step produces is the denominator of the density correction in
the next step (07_clustering_coeff).

Input
  ../../raw_data/CODEX_adata_dict.npy        AnnData for 31 samples (archive version only)
  ../../raw_data/tissue_boundary/*.geojson   tissue boundaries, used to check that the niche
                                             centres lie inside the tissue
  ../../data/niche_centers_analysis_2026-06-14.csv
                                             the niche centres sampled during the analysis
                                             (see below); SP_ID, cell_id
  ../../sample_metadata.tsv                  sample -> SP number, group

Output (committed)
  ../../data/niche_meta.csv                  one row per niche, 31 x 100 = 3100 rows
        sample, SP_ID, group, group_as_run, cell_id, index_ct_agg,
        num_cells (cells in the neighbourhood, including the index cell),
        num_tumor_cells
  ../../data/sample_density_factor.csv       one row per sample
        n_niches, density_factor (= mean of num_cells), num_cells_median, num_cells_sd

Method
  A niche is every cell within 40 um of an index cell, including the index cell itself.
  Index cells were chosen under three constraints:

    1. the cell must lie inside the tissue boundary eroded inwards by 40 um, so that the
       whole neighbourhood falls within the tissue -- otherwise niches at the tissue
       edge have an artificially low cell count because part of the disc is cut off;
    2. no two index cells are adjacent in the r40 graph (i.e. they are more than 40 um
       apart), which limits the overlap between niches;
    3. the index cell has r40 degree > 2, which excludes isolated cells.

  At most 100 centres were taken per sample, drawn at random from a candidate pool of ten
  times that many cells satisfying constraint 1; the first cell of the pool was accepted
  without the degree check.

  The niche centres used in the paper are the ones sampled during the analysis, on
  2026-06-14 (the run behind Fig. 4d, Fig. 4e and Extended Data Fig. 5c). That draw had no
  random seed, so a new draw would give different niches; the centres are therefore shipped
  as data/niche_centers_analysis_2026-06-14.csv and this step reads them instead of drawing
  again. It checks the three constraints on them and reports any centre that does not
  meet one.

  The neighbourhood graph is not shipped but rebuilt from obsm['spatial']; see the notes
  on data_loading.attach_radius_graph. num_cells is taken directly as "r40 degree + 1",
  which is equivalent to the original implementation's building of a sub-AnnData and
  summing value_counts (get_adata_sub returns exactly the index cell plus its r40
  neighbours) but much faster.

Differences from the original niche/00_sample_niche_w_boundary.py
  * The original drew the niches on every run; here the niches of the analysis run are
    read from data/ (see above), so the same input always gives the same result.
  * The original did adata[idx].copy() for every niche and saved it as npy, producing
    tens of GB of intermediates for what amounts to a cell count. No sub-objects are
    built here.
  * The original also computed pair counts (th_conn_matrix + get_pairs_in_niche), which
    belong to a different analysis and have nothing to do with the density correction;
    they are not computed here.
"""
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
from constants import load_samples, malig_celltype
from data_loading import (load_adata_dict, attach_radius_graph, legacy_sample_order,
                          tissue_boundary_path)

OUT_PATH = REPO / 'data'
CENTERS = REPO / 'data' / 'niche_centers_analysis_2026-06-14.csv'

RADIUS = 40            # um, the niche radius and also the minimum spacing between centres
N_NICHES = 100         # niches per sample
BOUNDARY_MARGIN = 40   # um, how far the tissue boundary is eroded inwards; None = no erosion
MIN_DEGREE = 2         # lower bound on the index cell's r40 degree (strictly greater)

# geojson coordinates are full-resolution pixels; obsm['spatial'] is in um. 0.5081 um/px
# is the calibration of this batch of CODEX images, carried over from the constant in
# the original niche/00_sample_niche_w_boundary.py.
GEOJSON_SCALE = 0.5081


def interior_mask(adata, sample_name, margin=BOUNDARY_MARGIN):
    """Return a boolean array marking cells inside the boundary eroded by margin."""
    if margin is None:
        return np.ones(adata.n_obs, dtype=bool)

    fp = tissue_boundary_path(sample_name)
    if not fp.exists():
        raise FileNotFoundError(f'{sample_name} has no tissue boundary file: {fp}')

    boundary = gpd.read_file(fp)
    # [a, b, d, e, xoff, yoff]: x' = a*x + b*y + xoff; pure scaling, so b=d=xoff=yoff=0
    boundary['geometry'] = boundary.geometry.affine_transform(
        [GEOJSON_SCALE, 0, 0, GEOJSON_SCALE, 0, 0])

    inner = boundary.geometry.buffer(-margin)
    inner = inner.union_all() if hasattr(inner, 'union_all') else inner.unary_union
    if inner.is_empty:
        raise ValueError(f'{sample_name}: tissue too small, empty after eroding {margin} um')

    coords = np.asarray(adata.obsm['spatial'])
    pts = gpd.GeoSeries(gpd.points_from_xy(coords[:, 0], coords[:, 1]))
    return pts.within(inner).to_numpy()


def check_centers(conn, mask, centers):
    """Count the centres that miss a constraint: outside the eroded boundary, adjacent to
    another centre, degree <= MIN_DEGREE (the first centre of a sample is exempt)."""
    sub = conn[centers][:, centers]
    return {
        'outside': int((~mask[centers]).sum()),
        'adjacent': int((sub.sum(axis=1) > 0).sum()),
        'low_degree': int(sum(conn[i].nnz <= MIN_DEGREE for i in centers[1:])),
    }


def main():
    OUT_PATH.mkdir(parents=True, exist_ok=True)
    samples = load_samples()
    centers_all = pd.read_csv(CENTERS, dtype=str)
    assert set(centers_all['SP_ID']) == set(samples.index), 'centres and samples differ'
    print(f'{len(samples)} samples, {N_NICHES} niches each (centres from {CENTERS.name}), '
          f'radius {RADIUS} um, boundary eroded by {BOUNDARY_MARGIN} um')

    print('loading raw_data/CODEX_adata_dict.npy ...', flush=True)
    adata_dict = load_adata_dict(samples.index)

    rows, missed = [], []
    for i, sample_name in enumerate(samples.index, 1):
        adata = adata_dict[sample_name]
        attach_radius_graph(adata, RADIUS)
        conn = adata.obsp[f'r{RADIUS}_connectivities'].tocsr()

        names = centers_all.loc[centers_all['SP_ID'] == sample_name, 'cell_id']
        pos = pd.Series(np.arange(adata.n_obs), index=adata.obs_names)
        absent = names[~names.isin(pos.index)]
        if len(absent):
            raise KeyError(f'{sample_name}: {len(absent)} niche centres are not cells of the sample')
        centers = pos[names].to_numpy()
        assert len(centers) == N_NICHES, (sample_name, len(centers))

        mask = interior_mask(adata, sample_name)
        miss = check_centers(conn, mask, centers)
        missed.append({'sample': sample_name, **miss})

        ct = adata.obs['ct_agg'].astype(str).to_numpy()
        is_malig = np.isin(ct, malig_celltype)
        obs_names = np.asarray(adata.obs_names)

        for idx in centers:
            members = conn[idx].indices          # r40 neighbours
            rows.append({
                'sample': sample_name,
                'cell_id': obs_names[idx],
                'index_ct_agg': ct[idx],
                'num_cells': len(members) + 1,              # includes the index cell
                'num_tumor_cells': int(is_malig[members].sum() + is_malig[idx]),
            })

        print(f'  [{i:2d}/{len(samples)}] {sample_name}: n={adata.n_obs:,} '
              f'-> {len(centers)} niches; centres missing a constraint: {miss}', flush=True)

        # the graphs are memory-hungry; drop them as soon as they are used
        del adata.obsp[f'r{RADIUS}_connectivities']
        del adata.obsp[f'r{RADIUS}_distances']

    niche = pd.DataFrame(rows)
    for col in ('SP_ID', 'group', 'group_as_run'):
        niche.insert(1, col, samples.loc[niche['sample'], col].to_numpy())

    g = niche.groupby('sample')['num_cells']
    factor = pd.DataFrame({
        'n_niches': g.size(),
        'density_factor': g.mean(),
        'num_cells_median': g.median(),
        'num_cells_sd': g.std(),
    })
    factor.index.name = 'sample'
    # One row per sample, in the order of the original sample keys (not sorted by SP ID):
    # that is the row order the committed table has, and step 07 takes its sample order
    # from this file.
    factor = factor.loc[legacy_sample_order(factor.index)].reset_index()
    for col in ('SP_ID', 'group', 'group_as_run'):
        factor.insert(1, col, samples.loc[factor['sample'], col].to_numpy())

    niche.to_csv(OUT_PATH / 'niche_meta.csv', index=False)
    factor.to_csv(OUT_PATH / 'sample_density_factor.csv', index=False)

    missed = pd.DataFrame(missed).set_index('sample')
    print(f'\ncentres missing a constraint, all samples: {missed.sum().to_dict()}')
    print(f'{len(niche)} niches, density_factor ranges '
          f'{factor["density_factor"].min():.1f} - {factor["density_factor"].max():.1f}')
    print(factor.groupby('group_as_run')['density_factor'].agg(['size', 'mean']).to_string())
    print(f'\nsaved -> {OUT_PATH / "niche_meta.csv"}')
    print(f'saved -> {OUT_PATH / "sample_density_factor.csv"}')


if __name__ == '__main__':
    main()
