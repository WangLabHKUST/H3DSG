# CODEX (multiplex immunofluorescence) spatial analysis

Spatial analysis of 31 CODEX samples of spinal cord glioma: cell type composition,
Shannon entropy density maps over 27.5 μm neighbourhoods, 40 μm neighbourhood
enrichment and clustering coefficient, observed/expected enrichment of pairwise cell
type adjacency, clustering coefficient corrected for local cell density, and mean
nearest-neighbour distances between cell types.

## How this directory is organised

Three layers; each layer reads only one kind of input.

| Layer | Reads | Writes | Runtime |
|---|---|---|---|
| `pipeline/01_build_adata/` | QuPath measurement tables + cell type tables | `raw_data/CODEX_adata_dict.npy` | hours, one sample at a time |
| `pipeline/02` – `pipeline/08` | `raw_data/CODEX_adata_dict.npy` (7.9 GB) and friends | the full-size intermediates in `data_full/`, the small tables in `data/` | minutes to hours, memory-hungry |
| `figures/` | `data/` only | figures | seconds |

**One rule runs through the whole repository, and it is the only one you have to
remember:**

> **No script under `figures/` ever opens `raw_data/`.** Every figure is drawn from a
> small table in `data/`, and every table in `data/` is written by `pipeline/`.

So redrawing a figure needs nothing but this repository. Verified:
`figures/entropy_density_by_group.R` rerun from the three 400-row grid CSVs in `data/`
produces 5 PNGs that are **byte-identical** to the ones already in the repository.

## Two versions

| Version | Contents | Size | Purpose |
|---|---|---|---|
| **Archive** | this entire directory, including `raw_data/` and `data_full/` | ~9.9 GB | hand the whole folder to a collaborator; every step from raw input to figure is present |
| **GitHub** | the same code + the small tables in `data/` + the PNGs in `figures/`, without those two directories | ~5 MB | public repository |

There is only one copy of the code. What the GitHub version leaves out is decided by
`.gitignore` alone (there is no second list): `raw_data/`, `data_full/`, `results/`,
`logs/`, the h5ad / npy / rds files and the non-PNG figure formats. Nothing committed may
exceed GitHub's 100 MB per-file limit.

## Data tiers

| Tier | Contents | Size | Where |
|---|---|---|---|
| 0 | `CODEX_adata_dict.npy`: AnnData for 31 samples (X = channel intensities after positivity calling, `obsm` spatial coordinates, the `r50_adj` neighbourhood graph in `obsp`, the 10-class `ct_agg` in `obs`) | 7.9 GB | `raw_data/`, archive only |
| 0 | `cell_type_2.parquet`: the 17-class fine-grained annotation for 7.09 M cells (a subdivision of `ct_agg`, used for nearest-neighbour distances) | 257 MB | `raw_data/`, archive only |
| 0 | `tissue_boundary/*.geojson`: 37 tissue boundaries, used to keep niche centres inside the tissue | 5.2 MB | `raw_data/`, archive only |
| 0 | `sample_key.tsv`: the table that translates the sample keys used inside `raw_data/` and `data_full/` to SP IDs (see "Sample keys" below); `common_data/`: the small input tables of `pipeline/00_sample_metadata/` | < 1 MB | `raw_data/`, archive only |
| 1 | full-size cell-level intermediates: `entropy_cell_level.parquet` (7.09 M rows), `neighbor_counts_r27_5.parquet`, the downstream entry point `entropy_cell_level_sampled.csv` (1.32 M rows), `sampled_cell_ids.csv` | ~1 GB | `data_full/`, archive only |
| 2 | the small tables the figures are drawn from | ~2.5 MB | `data/`, committed |

Tier 1 is kept so that changing a filter or a down-sampling rule later does not mean
touching the 7.9 GB raw data again. Tier 2 is the only input `figures/` has.

## Quick start

Redraw the entropy density maps, with nothing to download:

```bash
cd figures
Rscript entropy_density_by_group.R
```

Rerun the steps that need raw data (archive version):

```bash
python pipeline/00_sample_metadata/01_build_sample_metadata.py   # seconds
python pipeline/02_entropy/01_cell_level_stats.py                # hours, needs a lot of RAM
python pipeline/02_entropy/02_kde_grid.py
python pipeline/03_proportion/01_celltype_counts.py
python pipeline/04_interaction/01_neighborhood_enrichment.py
python pipeline/04_interaction/02_pairwise_enrichment.py
python pipeline/05_marker_expression/01_marker_summary_by_celltype.py
python pipeline/06_niche_density/01_sample_niches.py
python pipeline/08_min_distance/01_nearest_neighbour_distance.py       # minutes

# These two do not read the 7.9 GB raw data; they read what the steps above wrote
# into data/, and take seconds.
python pipeline/08_min_distance/02_group_test.py
python pipeline/07_clustering_coeff/01_density_adjusted.py             # needs 04 and 06 first
```

## Randomness and reproducibility

Four places use random numbers, and **none of the original implementations set a
seed**, so a new run would not reproduce the original outputs bit for bit. Two of them
therefore do not draw again but read what was drawn during the analysis; the other two are
now seeded (`SEED = 42`). Either way, the same input always gives the same result:

| Where | What is random | What it affects |
|---|---|---|
| `06_niche_density/01` | drawing niche centres | `density_factor`, and through it `count_adj` in step 07 — this one does **not** draw again; it reads the niche centres sampled during the analysis (2026-06-14, the run behind Fig. 4d, 4e and Extended Data Fig. 5c), shipped as `data/niche_centers_analysis_2026-06-14.csv` |
| `04_interaction/01` | label permutation in squidpy's neighbourhood enrichment | `zscore` and the p values derived from it; `count` itself is deterministic and unaffected |
| `02_entropy/01` | cell down-sampling | the entropy line — this one does **not** resample; it reuses the cell IDs that were written to disk at the time, see `pipeline/README.md` |
| `figures/panel_helpers.py`, `figures/clustering_coeff_boxplot_PL_vs_TL.py` | the jitter offsets of seaborn's `stripplot` | only where the points land inside a box; no statistic depends on it. Verified: the three boxplot figures now come out byte-identical on a rerun, whereas before seeding every run produced a different PNG |

In `06_niche_density/01` the niche centres file holds, for each of the 31 samples, the
100 index cells of the analysis run (`SP_ID`, `cell_id` = the cell's name in the AnnData).
The script recomputes every niche from them (the cells within 40 um of the index cell) and
checks that the centres meet the three sampling constraints described in its header. How
much a different draw would change the result: the standard error of each sample's
`density_factor` is about 4.3% of its mean (7.6% at worst), roughly ±0.06 on
`log2(count_adj)`, against a median between-sample SD within a group of 1.105, i.e. about
**0.3% of the variance**; no conclusion of the group test depends on the draw.

Each script states its own inputs and outputs in its file header;
`pipeline/README.md` has the summary table.

## Samples and groups

`sample_metadata.tsv` is the single sample manifest: 35 rows, of which **31 enter the
analysis**. Code always reads it through `constants.load_samples()`; never hard-code a
sample list in a script.

### Sample names

Every sample is named by the Sample ID of Table S1 of the paper (`SP11`, ...), everywhere
in the published files. Three columns of `sample_metadata.tsv` carry it:

| Column | Example | Used by |
|---|---|---|
| `CODEX_ID` | `SP11` | **all file names, variable names and intermediate-table keys**, and the keys of `adata_dict` as seen by the code |
| `Original_ID` | `SP11` | kept for compatibility with older scripts; the hospital record number is not published |
| `SP_ID` | `SP11` | the Sample ID in Table S1 of the paper — **use this in figures and in anything shown externally** |

The three columns are identical; they are kept so that a script written against any of
them keeps working. The record-number to SP-number correspondence of all 35 samples has
been checked against Table S1 in both directions: 35/35 agree.

### Sample keys (archive version only)

`raw_data/` and `data_full/` were written before the samples were renamed and are keyed by
the original sample keys (the keys of the dict in `CODEX_adata_dict.npy`, the `sample`
column of the parquet files, the names of the `tissue_boundary/*.geojson` files). They have
not been rewritten. The table `raw_data/sample_key.tsv` (columns `raw_key`, `SP_ID`; it is
written by `pipeline/00_sample_metadata/01_build_sample_metadata.py` and not published)
translates between the two; `data_loading.py` is the only code that reads it, and every
script that needs a translation calls the functions below:

| Function | What it does |
|---|---|
| `load_adata_dict()`, `load_cell_type_2()`, `read_data_full()` | read `raw_data/` / `data_full/` and return SP IDs |
| `write_data_full()`, `to_raw_key()`, `tissue_boundary_path()` | go the other way where a file name or a stored column has to carry the original key |
| `legacy_sample_order()` | sorts SP IDs in the order of the original keys, the row order of `sample_key.tsv` |

Everything else sees SP IDs only. If `raw_data/sample_key.tsv` is missing, these
functions stop with an error that says so. Two places depend on the original key or its
order, and both would change the results if they were switched to SP IDs without care:

- **the niche seed** of `06_niche_density/01` is a hash of the sample's original key, so
  the niches drawn are the same as in the committed tables;
- **the row order** of `data/sample_density_factor.csv` and of
  `data/clustering_coefficient_adjusted.csv` follows the order of the original keys, not
  the alphabetical order of the SP IDs (a sort on the sample column used to give that
  order; `06` now writes it explicitly through `legacy_sample_order()` and `07` follows
  `sample_density_factor.csv`).

### PL / TL grouping

**The authoritative source is `Group/Predict Group` in Table S1 of the paper** (training
labels for the discovery cohort; for the validation cohort, the call of the SCG
classifier, written as e.g. `TL-SCG classifier` — strip the suffix to get PL/TL).
`sample_metadata.tsv` carries two columns with the same content but different
provenance:

| Column | From | Distribution over the 31 samples |
|---|---|---|
| `group` | **Table S1**, authoritative | PL 14 / TL 17 |
| `group_as_run` | `constants.g1_samples` / `g2_samples`, i.e. the grouping that actually produced `data/` and `figures/` | PL 14 / TL 17 |

**The two columns agree completely**, so the results in `data/` and `figures/` match the
paper and do not need to be rerun. The generating script asserts this: if Table S1 ever
changes a group while the existing outputs still use the old one, the script fails
loudly instead of passing silently.

```python
load_samples()                              # group_as_run, the default
load_samples(group_column='group')          # Table S1; currently equivalent
```

> **Warning: do not use the `Group` column of
> `common_data/SP_ID_pathology_ID_mapping.xlsx` to assign groups.** It agrees with
> Table S1 for 34 of the 35 CODEX samples, but records **SP41** as `PL`, whereas Table S1
> has `TL-SCG classifier`, the classifier gives `Prob = 0.5173 > 0.5`, and `samples.csv`
> on the CODEX side also puts it in Group2 (TL). **SP41 is TL**; that one cell of the
> spreadsheet is stale. `pipeline/00_sample_metadata/` takes only the ID columns from it,
> and uses its `Group` column solely to compare against Table S1 and print a warning.
> The latest Table S1 is the authority for every group.

### The 4 samples that do not enter the analysis

| `SP_ID` | Reason for exclusion |
|---|---|
| `SP02` | low_quality |
| `SP50` | low_quality |
| `SP42` | recurrent |
| `SP43` | recurrent |

The exclusion happens two steps upstream (`sample_lst_33` and `sample_lst_31` in
`codes/cleaned_202604/codex/utils.py`); these samples **are simply not present** in
`CODEX_adata_dict.npy` (verified: the full-size layer holds exactly 31 samples and
7,092,596 cells).

## Cell type naming

The 10 classes in `constants.celltype_order` are the names used in the submitted
manuscript. During development (`codes/all`, `codes/cleaned_202604`) a different set of
short names was in use; the difference is a `_like` suffix on the malignant states:

```
old OPC_OC     (malignant state) -> new OPC_OC_like
old AC_MES     (malignant state) -> new AC_MES_like
old Oligo_pro  (normal cell)     -> new OPC_OC
```

**Careful: `OPC_OC` does not refer to the same population in the two sets of names**, so
when porting code from an old script, first establish which set it uses.
`CODEX_adata_dict.npy` already uses the new names; `normalize_celltype()` in
`pipeline/02_entropy/01_` converts the old names over, using the presence of
`Oligo_pro` as the test.

The canonical display names used in figures are in `constants.celltype_display`
(`OPC/OC-like`, `AC/MES-like`, `T/NK` and so on).

## Runtime environment

Python 3.12 (anaconda base) + R 4.5.1. Tested versions are recorded in `env/`.

Two environment-related traps:

- **Loading `CODEX_adata_dict.npy` requires a monkey patch.** It was pickled by an
  older version of anndata; under anndata 0.10.x a plain `np.load` raises
  `KeyError: '_adata_ref'`. `data_loading.load_adata_dict()` applies the patch, and every
  script that needs the dict goes through it. There is no way around this on a newer
  environment either; see section 3.3 of `notes/Entropy_data_provenance.md`.
- **All text in figures is English.** The Windows `pdf()` device cannot embed CJK
  glyphs and renders them as garbage.

## Known issues

1. **The three scripts in `pipeline/01_build_adata/` are archived verbatim and cannot be
   run as they stand.** They are notebook-style, handle one sample per run (with
   `sample_name = 'SP08'` written into the file), and their inputs (QuPath tables, h5ad
   files, cell type tables) are expected under `raw_data/`, where the archive version does
   not ship them because of their size. They are kept so that the
   "raw input to `.npy`" leg is documented; the `.npy` already exists, so normally there
   is nothing to rerun. See `pipeline/01_build_adata/README.md`.
2. **The original `Spatial_clustering.py` had four defects that stopped it running**
   (`process_sample()` was never called, `zscore_df` was merged twice, a typo
   `alternative='ryy-one_sided'`, and a missing `sample` column). These are fixed in
   `pipeline/04_interaction/01_neighborhood_enrichment.py` and each one is recorded in
   that file's header. **This line therefore has no pre-existing output to reconcile
   against** — it never produced a result before.
3. **No multiple-testing correction is applied** in the clustering-coefficient or
   nearest-neighbour-distance group tests. The tables carry raw p values and the figures
   map them to dot size; correct them yourself (BH or similar), or restrict attention to
   pre-specified pairs, before drawing conclusions in text.

## Provenance

This directory was assembled from `Multiplex-IHC/spatial_analysis/codes/submit/` (the
version in active use) rather than from `GitHub_ref/CODEX_Analysis/`, which is an
earlier partial snapshot missing `Preprocessing.py`, `Entropy.py`, `data/`, `figures/`
and the provenance notes. The upstream adata-building scripts come from
`codes/cleaned_202604/codex/`.

`notes/Entropy_data_provenance.md` (29 KB) is a layer-by-layer provenance record of the
entropy line, including the reproducibility checks run at each step. Worth reading
before changing anything on that line.
