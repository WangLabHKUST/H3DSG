# data — the small tables the figures are drawn from

The scripts under `figures/` read only this directory, so every table here is committed
and redrawing a figure requires no download. Every table is written by `pipeline/`, apart
from the four listed under "Tables copied from the original analysis"; not one of them was
edited by hand, so the others can be regenerated and diffed.

The bulky intermediates live in `../data_full/` (~1 GB, not committed) and the raw data
in `../raw_data/` (7.9 GB, not committed); both exist only in the archive version.

## Table index

Grouped by the `pipeline/` step that produces them.

### 02 entropy

| File | Rows | Columns | Written by |
|---|---|---|---|
| `grid_data_entropy_tumor_vs_non_malig_entropy_Group1.csv` | 400 | `entropy_tumor`, `non_malig_entropy`, `density` | `pipeline/02_entropy/02_kde_grid.py` |
| `grid_data_entropy_tumor_vs_non_malig_entropy_Group2.csv` | 400 | as above | as above |
| `grid_data_entropy_tumor_vs_non_malig_entropy_diff.csv` | 400 | as above, `density` = Group1 − Group2 | as above |

These three hold a 2-D Gaussian kernel density estimate on a 20×20 grid. They are
**already post-KDE and contain no cell-level information**, so they are safe to publish.
`figures/entropy_density_by_group.R` reads only these, and has been verified to
reproduce the 5 PNGs in the repository byte for byte.

### 03 cell composition

| File | Rows | Written by | Consumed by |
|---|---|---|---|
| `celltype_counts.csv` | 31 | `pipeline/03_proportion/01_celltype_counts.py` | `figures/celltype_proportion_stacked_bar.py` |
| `celltype_proportions.csv` | 31 | as above | as above |

### 04 neighbourhood enrichment and pairwise enrichment

| File | Rows | Written by | Consumed by |
|---|---|---|---|
| `neighborhood_enrichment_by_sample.csv` | 3100 | `pipeline/04_interaction/01_neighborhood_enrichment.py` | `pipeline/07` |
| `clustering_coefficient_group_test.csv` | 10 | as above | — |
| `pairwise_enrichment_by_sample.csv` | 3100 | `pipeline/04_interaction/02_pairwise_enrichment.py` | — |
| `pairwise_enrichment_pooled.csv` | 100 | as above | `figures/interaction_enrichment_heatmap.py` |

3100 = 31 samples × 100 ordered cell type pairs (10 × 10).

### 05 marker expression

| File | Rows | Written by | Consumed by |
|---|---|---|---|
| `marker_summary_by_celltype.csv` | 230 | `pipeline/05_marker_expression/01_marker_summary_by_celltype.py` | `figures/celltype_marker_dotplot.py` |

230 = 10 cell types × 23 markers, pooled over all samples.

### 06 niche density

| File | Rows | Columns | Written by |
|---|---|---|---|
| `niche_centers_analysis_2026-06-14.csv` | 3100 | `SP_ID`, `cell_id` | input, not written by `pipeline/` (see below) |
| `niche_meta.csv` | 3100 | `sample`, `group_as_run`, `group`, `SP_ID`, `cell_id`, `index_ct_agg`, `num_cells`, `num_tumor_cells` | `pipeline/06_niche_density/01_sample_niches.py` |
| `sample_density_factor.csv` | 31 | `n_niches`, `density_factor`, `num_cells_median`, `num_cells_sd` | as above |

`niche_centers_analysis_2026-06-14.csv` holds the niche centres sampled during the analysis
(2026-06-14; the run behind Fig. 4d, 4e and Extended Data Fig. 5c): 100 index cells per
sample, `cell_id` being the cell's name in the AnnData. That draw had no random seed, so the
centres are kept as data and step 06 reads them instead of drawing again.

`density_factor` is the mean cell count over that sample's 100 niches, and serves as the
denominator of the density correction in step 07.

### 07 density-adjusted clustering coefficient

| File | Rows | Written by | Consumed by |
|---|---|---|---|
| `clustering_coefficient_adjusted.csv` | 310 | `pipeline/07_clustering_coeff/01_density_adjusted.py` | `figures/clustering_coeff_violin_raw_vs_adj.py`, `figures/clustering_coeff_boxplot_PL_vs_TL.py` |
| `clustering_coefficient_adjusted_group_test.csv` | 10 | as above | as above |

310 = 31 samples × 10 cell types.

### 08 nearest-neighbour distance

| File | Rows | Written by | Consumed by |
|---|---|---|---|
| `min_distance_by_sample.csv` | 20,398 | `pipeline/08_min_distance/01_nearest_neighbour_distance.py` | `figures/min_distance_boxplot_*.py` |
| `min_distance_group_test.csv` | 658 | `pipeline/08_min_distance/02_group_test.py` | `figures/min_distance_dotmatrix_*.py` |

20,398 = 31 samples × 658 pairs, and 658 is the sum of the squared class counts of the
three variants (no_map 17² + map_tum 15² + map_tum_map_mye 12²). The `variant` column
says which cell type merging scheme a row uses; its values correspond one to one with
`VARIANTS` in `pipeline/08_min_distance/01`.

**Neither table is corrected for multiple testing**; both carry raw p values. The figures
only map them to dot size or print them in a panel title. Correct them yourself before
drawing a conclusion in text, or restrict attention to pre-specified pairs.

### Tables copied from the original analysis (no pipeline step)

These four tables are the exception to the rule above: no step in `pipeline/` writes them. They are copies of tables
from the original analysis, with SP IDs. `group` follows Table S1, and `group_as_run` is the group used in that run.

| File | Rows x columns | Columns | Consumed by |
|---|---|---|---|
| `entropy_region_cell_proportion_by_sample.csv` | 31 x 13 | `SP_ID`, `group_as_run`, `group`, `cells_{niche1..niche4,none}`, `prop_{niche1..niche4,none}` | `figures/entropy_region_cell_proportion_boxplot_PL_vs_TL.py` |
| `entropy_region_pairs_per_niche_area{1,2}.csv` | 31 x 68 | `SP_ID`, `group_as_run`, `group`, and 65 columns of cells or pairs per niche; `A...B` is a pair of A and B | `figures/entropy_region_pairs_per_niche_boxplot_PL_vs_TL.py` |
| `tnk_subset_proportion_by_sample.csv` | 31 x 9 | `SP_ID`, `group_as_run`, `group`, `Tcell_CD4`, `Tcell_CD8`, `Treg_CD4`, `Treg_CD8`, `Treg_other`, `NK`; share of all cells of the sample, an empty `NK` means the class is absent | `figures/cd8_t_cell_proportion_boxplot_PL_vs_TL.py` |

## Join keys

- `sample` / `CODEX_ID` join to `../sample_metadata.tsv`. The two names mean the same
  thing (`sample` in long tables, the row index `CODEX_ID` in wide ones), with values of
  the form `SP11`.
- When a figure needs to label samples, use `SP_ID` (the number used in Table S1 of the
  paper); the tables above already carry that column.
- Two tables are ordered by sample: `sample_density_factor.csv` lists the samples in the
  order of the original sample keys (see "Sample keys" in the top-level README), and
  `clustering_coefficient_adjusted.csv` follows it within each cell type and group. That
  is not the alphabetical order of the SP IDs; the pipeline writes this order explicitly,
  so a rerun gives the same rows in the same places. The other tables list the samples in
  the order of `../sample_metadata.tsv`.
- `group` takes the values `Group1` / `Group2`, and only becomes `PL` / `TL` at figure
  time via `constants.group_labels`. Note that the `group_as_run` column in the same
  table holds the literal **`PL` / `TL`** (it is the original column from
  `sample_metadata.tsv`), whereas `group` has already been mapped to `Group1` /
  `Group2` by `constants.load_samples()`. Filter on `group`.
- The nearest-neighbour-distance tables use `source` / `target` for the ordered pair,
  with values from the 17 classes of `cell_type_2` (or the merged `Tumor` / `Myeloid`).
  That is not the same vocabulary as the `cell_type` column (10 `ct_agg` classes) of the
  tables above; see `constants.celltype2_order`.

## Naming rule

**A file name says what the table contains, never which figure it feeds.** Figure
numbers move between revisions, and a file called `fig3b_*.csv` is wrong the moment a
panel is reordered.
