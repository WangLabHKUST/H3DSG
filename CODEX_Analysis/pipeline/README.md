# Pipeline

Run in numbered order. Each script states its inputs and outputs in its own file header,
and each one runs standalone — no script is `import`ed or `source`d by another, so
reading one script never means following a call chain into other files.

| Step | Script | Input | Output |
|---|---|---|---|
| 00 | `00_sample_metadata/01_build_sample_metadata.py` | `../raw_data/common_data/`: the paper's `Table S1_0921.xlsx`, `SP_ID_pathology_ID_mapping.xlsx`, `samples.csv` (archive version only) | `../sample_metadata.tsv`, 35 rows; `../raw_data/sample_key.tsv`, the table that turns the sample keys of `raw_data/` and `data_full/` into SP IDs |
| 01 | `01_build_adata/` | QuPath measurement tables + cell type tables + the existing graph structures | `../raw_data/CODEX_adata_dict.npy` — **archived verbatim, not runnable as it stands**; see that directory's README |
| 02 | `02_entropy/01_cell_level_stats.py` | `../raw_data/CODEX_adata_dict.npy` | `../data_full/entropy_cell_level.parquet`, `neighbor_counts_r27_5.parquet`, `entropy_cell_level_sampled.csv` |
| 02 | `02_entropy/02_kde_grid.py` | `../data_full/entropy_cell_level_sampled.csv` | `../data/grid_data_*_{Group1,Group2,diff}.csv`, 400 rows each |
| 03 | `03_proportion/01_celltype_counts.py` | `../raw_data/CODEX_adata_dict.npy` | `../data/celltype_counts.csv`, `celltype_proportions.csv` |
| 04 | `04_interaction/01_neighborhood_enrichment.py` | `../raw_data/CODEX_adata_dict.npy` | `../data/neighborhood_enrichment_by_sample.csv`, `clustering_coefficient_group_test.csv` |
| 04 | `04_interaction/02_pairwise_enrichment.py` | `../raw_data/CODEX_adata_dict.npy` | `../data/pairwise_enrichment_by_sample.csv`, `pairwise_enrichment_pooled.csv` |
| 05 | `05_marker_expression/01_marker_summary_by_celltype.py` | `../raw_data/CODEX_adata_dict.npy` | `../data/marker_summary_by_celltype.csv` |
| 06 | `06_niche_density/01_sample_niches.py` | `../raw_data/CODEX_adata_dict.npy`, `../raw_data/tissue_boundary/*.geojson` (one per sample), `../data/niche_centers_analysis_2026-06-14.csv` (niche centres sampled during the analysis) | `../data/niche_meta.csv` (3100 rows), `sample_density_factor.csv` (31 rows) |
| 07 | `07_clustering_coeff/01_density_adjusted.py` | the outputs of steps 04 and 06 (**does not read the `.npy`**) | `../data/clustering_coefficient_adjusted.csv`, `clustering_coefficient_adjusted_group_test.csv` |
| 08 | `08_min_distance/01_nearest_neighbour_distance.py` | `../raw_data/CODEX_adata_dict.npy`, `../raw_data/cell_type_2.parquet` | `../data/min_distance_by_sample.csv` |
| 08 | `08_min_distance/02_group_test.py` | the output of the row above (**does not read the `.npy`**) | `../data/min_distance_group_test.csv` |

Steps 02 through 06 and `08/01` all read the 7.9 GB `.npy` (archive version only), so
each of those scripts checks for it up front and, if it is missing, fails with a pointer
to the ready-made small tables in `../data/`. Steps `07` and `08/02` never read the
`.npy`; they read what the earlier steps left in `../data/`.

`raw_data/` and `data_full/` are keyed by the sample keys they were first written with,
not by SP IDs. Every script that reads or writes them does so through `../data_loading.py`
(`load_adata_dict()`, `load_cell_type_2()`, `read_data_full()`, `write_data_full()`,
`tissue_boundary_path()`), which translates with `../raw_data/sample_key.tsv`; the scripts
themselves see SP IDs only. See "Sample keys" in the top-level README. Two details of
`06_niche_density/01` follow from this: the seed of each sample is a hash of its key in
`raw_data/` (so the committed niches are reproduced), and the rows of
`sample_density_factor.csv` are written in the order of those keys; step 07 keeps that
order.

There are only three ordering dependencies: inside step 02 (the output of
`01_cell_level_stats.py` is the input of `02_kde_grid.py`), step 07 needs both 04 and 06
to have finished, and `08/02` needs `08/01`. Steps 03, 04, 05, 06 and `08/01` are
mutually independent and depend only on the `.npy` and `../sample_metadata.tsv`.

## Four neighbourhood graphs — do not mix them up

Each AnnData in the `.npy` was built with four graphs in `obsp` (by
`01_build_adata/00_modify_adata.py`). Different steps use different ones:

| Graph | How it is built | Used by |
|---|---|---|
| `de` | Delaunay triangulation | unused in this repository |
| `r40` | 40 μm radius | `04_interaction/01`, `06_niche_density/01` — **not present in the `.npy` either**, see below |
| `r50` / `r50_adj` | 50 μm radius | `04_interaction/02_pairwise_enrichment.py` (then truncated to the 5 nearest neighbours) |
| `r27_5` | 27.5 μm radius | entropy — but **not present in the `.npy`**, see below |

**Both `r40` and `r27_5` are rebuilt** rather than read from `obsp`: the `obsp` of
`CODEX_adata_dict.npy` holds `r50_adj` only. The table above describes the four graphs
that `01_build_adata/00_modify_adata.py` wrote into `_modified.h5ad`; the archived `.npy`
kept only `r50_adj` when it was exported.

`r40` is rebuilt from `obsm['spatial']` by `data_loading.attach_radius_graph(adata, 40)`
using `sklearn.radius_neighbors_graph(radius=40, include_self=False)`
(`include_self=False` matches squidpy's `spatial_neighbors(radius=...)`, i.e. no self
loops). It has been verified **edge for edge identical** to the `r40_connectivities` in
`_modified.h5ad` (diff_nnz = 0 on three samples, the largest of which has 2,835,714
edges). The check itself is kept in `data_loading.verify_radius_graph()`.

**The 27.5 μm entropy graph is likewise rebuilt**, not read from `obsp`:
`02_entropy/01_cell_level_stats.py` rebuilds it from `obsm['spatial']` with
`sklearn.radius_neighbors_graph(radius=27.5, include_self=True)`, and this too is
**edge for edge identical** to the `r27_5_connectivities` in the original
`_modified.h5ad` (diff_nnz = 0). `include_self=True` corresponds to the `adj.setdiag(1)`
of the original `mask_diagonal=False` implementation, i.e. the index cell counts towards
its own neighbourhood, so `total_num >= 1`. The number 27.5 μm is the diameter of a
Visium spot.

## The down-sampling on the entropy line is not reproducible

The original script called `.sample()` without `random_state`, so **the same set of cells
cannot be drawn again**. The cell IDs that were kept at the time were written to disk,
however, and are now frozen in `../data_full/sampled_cell_ids.csv` (1,320,620 IDs,
verified to hit `obs_names` of the `.npy` with 100% coverage).
`01_cell_level_stats.py` reuses that list and does **not** resample. Choosing a different
sample means accepting that the results change — in particular the colour breaks of the
difference map in `figures/entropy_density_by_group.R` (now generated dynamically from
the data range; the original script hard-coded `-10.05` / `1.47`).

The filter thresholds match the original
`20260129_entropy_cal_clean_v3_visium_spot_size.py`: `tumor_num >= 2`,
`non_malig_num >= 2`, `total_num <= 20`. Of 7.09 M cells, 2.13 M pass the filter, and the
ID list above down-samples those to 1.32 M (660,310 per group, deliberately equal).

## Which group column to use

Every script takes its samples and groups from `constants.load_samples()`, which defaults
to the `group_as_run` column of `sample_metadata.tsv` — the grouping that actually
produced `../data/` and `../figures/`. It **agrees completely** with the authoritative
`group` column (taken from Table S1 of the paper), so the two are currently equivalent
and the existing results match the paper. See "Samples and groups" in the top-level
README, which also records why the `Group` column of
`SP_ID_pathology_ID_mapping.xlsx` must not be used to assign groups.
