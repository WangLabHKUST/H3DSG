# Provenance of the entropy figures (density_r.Group*.auto.png / density_r_diff.png)

> Investigated 2026-09-28. **This document is a record only; no existing file was changed
> while writing it.**
> Every conclusion was verified by actually reading the file or loading the data, not
> inferred. Verified points are marked "✅ verified".

---

## 0. What the three figures actually are

| Figure | Generating script | Direct input |
|---|---|---|
| `density_r.Group1.auto.png`<br>`density_r.Group2.auto.png` | `codes/submit/Entropy.R` (= `codes/all/20260129_show_contour_group.R`) | `grid_data_entropy_tumor_vs_non_malig_entropy_{Group1,Group2}.csv` |
| `density_r_diff.png` | `codes/all/20260129_show_contour_diff.R` | `grid_data_entropy_tumor_vs_non_malig_entropy_diff.csv` |

What is plotted: for the neighbourhood of **radius 27.5 μm** (the diameter of a Visium
spot) around each cell,
x axis = the Shannon entropy of the **malignant cell states** in that neighbourhood
(`entropy_tumor`),
y axis = the Shannon entropy of the **non-malignant cell types** in the same
neighbourhood (`non_malig_entropy`),
colour = a 2-D **Gaussian kernel density estimate (gaussian_kde)** evaluated on a 20×20
grid (`density`).
The diff figure is `zi(Group1) - zi(Group2)`.

The grid csv itself has only **400 rows** (a 20×20 grid) plus a header, with the columns
`entropy_tumor, non_malig_entropy, density`.
**It is already post-KDE and contains no cell-level information.**

---

## 1. The full data chain (bottom up)

### 1.0 Further upstream: where the h5ad files come from

> Source: the scripts under `codes/cleaned_202604/codex/` (brought in from another
> session on 2026-09-28 and checked against the code).

```
QuPath measurements  all_samples_v1/csv/{sample}.csv ──┐
cell type tables     20251210_dt/cell_type_csv/{sample}.csv ──┤
existing graphs      h5ad/{sample}_full_raw.h5ad (.obsp) ──┘
        │
        ▼  codex/00_create_adata.py   <- where ad.AnnData(...) is actually called
   {sample}_masked_counts.h5ad
        │
        ▼  codex/00_modify_adata.py   <- filtering + building the spatial graphs
   {sample}_masked_counts_no_other.h5ad          (drops ct_agg in {Other, Tumor_other})
   {sample}_masked_counts_no_other_modified.h5ad (builds the de / r40 / r50 / r27_5 graphs)
        │
        ▼  codex/00_update_adata_obs.py
   h5ad/adata_obs_dict.npy   (stores only the three columns ct_agg / cell_type / cell_type_2)
```

**What each part of the AnnData means** (useful when writing Methods):

| Part | Definition |
|---|---|
| `X` | **masked counts** = the intensity of 23 proteins x the corresponding binary positivity column; a negative cell is set to zero for that protein |
| `obs['r']` | nuclear perimeter / 2π |
| `obs['l']` | the mean of the major and minor axes |
| `obs['r_cell']` | the mean of the two above |
| `obsm['spatial']` | `Centroid_X_µm`, `Centroid_Y_µm` -> **the coordinate unit is indeed µm**, consistent with the measurement in section 3.1 |
| `obsp` | carried over from `_full_raw.h5ad` as is, with four graphs appended by `00_modify_adata.py` |

Note: `_full_raw.h5ad` is **not** produced inside `cleaned_202604`; it comes from the
earlier `codes/all/20251024_write_adata.py:137` (written from an even earlier
`adata_dict.npy`, 35 samples). So strictly speaking the ultimate origin is in
`codes/all/`.

**Why samples were excluded** (`codex/utils.py:11-12`; until this was found, the
exclusion was known but not its reason):

```python
sample_lst_33 = [f for f in sample_lst_full if f not in ['SP02', 'SP50']]  ### low-quality
sample_lst_31 = [f for f in sample_lst_33 if f not in ['SP42', 'SP43']]    ### recurrent
```

- `SP02`, `SP50` -- **failed quality control**
- `SP42`, `SP43` -- **recurrent samples**

-> `sample_lst_31` agrees **exactly** with the entropy chain and with the 31 samples of
`CODEX_adata_dict.npy` ✅

### 1.1 The entropy chain itself

```
[L0] original h5ad: 20250929/h5ad/{sample}_masked_counts_no_other_modified.h5ad   (35 sample files)
       │  obsp already holds six spatial graphs: r27_5 / r40 / r50 / r50_adj / spatial / de
       │  obs['ct_agg'] = the old naming (AC_MES, OPC_OC, Oligo_pro, ...)
       ▼
[L1] codes/all/20260129_process_all.py
       │  aggregate_neighbors_counts(graph_key='r27_5_connectivities',
       │                             groupby_column='ct_agg', mask_diagonal=False)
       │  note mask_diagonal=False -> adj.setdiag(1), so the index cell counts towards
       │  its own neighbourhood
       ▼
     {sample}_neighbor_counts.csv / {sample}_neighbor_prop.csv   (31 samples x 2)
       directory: .../20260106_show_niche/all_samples_full_27_5/
       ▼
[L2] codes/all/20260129_entropy_cal_clean_v3_visium_spot_size.py  <- the core script
       │  (1) reads obs_dict_modified.npy to recover r/l/r_cell/r50_*/ct_agg
       │  (2) computes entropy_tumor / entropy_all / tumor_num / total_num /
       │      dominated_tumor_state per sample
       │  (3) merges -> all_samples_meta_data.csv (674 MB, 3,624,857 rows)
       │      and computes non_malig_entropy / non_malig_num / entropy_major
       │  (4) filter: non_malig_num>=2 & tumor_num>=2 & total_num<=20
       │  (5) DOWN-SAMPLING (see section 2)
       │  (6) gaussian_kde, nbins=20 -> xi/yi/zi
       ▼
     20260129_data_dict_filtered_normal_equal.npy   (only 19,909 bytes: just the 20x20 grid)
       ▼
[L3] exported at the end of the same script -> grid_data_*_{Group1,Group2,diff}.csv  (400 rows each)
       ▼
[L4] Entropy.R / 20260129_show_contour_diff.R -> PNG
```

### Input files at each layer (all under `all_samples_full_27_5/`)

| File | Size | Notes |
|---|---|---|
| `all_samples_meta_data.csv` | 674 MB | **the cell-level main table**, 16 columns, 3,624,857 rows |
| `all_samples_meta_data_filtered.csv` | 674 MB | its md5 is **identical to the file above**: a copy of the same content, not a filtered result (the code that wrote it is commented out in the script) |
| `all_samples_neighbor_counts_filtered.csv` | 211 MB | neighbourhood counts for the 10 cell types |
| `all_samples_neighbor_prop_filtered.csv` | 398 MB | the above, normalised row-wise |
| `20260126_max_group1_index_filtered_non_malig_equal.csv` | 59.6 MB | **the cell IDs kept after down-sampling**, 1,320,620 rows |
| `20260127_data_dict_filtered_normal_equal.npy` | 19,909 B | md5 **identical** to the `20260129_` version |
| `20260129_data_dict_filtered_normal_equal.npy` | 19,909 B | the KDE grid result |

The 16 columns of `all_samples_meta_data.csv`:
`Object_ID, r, l, r_cell, r50_adj_degree, r50_degree, ct_agg, sample, entropy_tumor, entropy_all, tumor_num, total_num, dominated_tumor_state, group, entropy_major, non_malig_entropy`

---

## 2. The sampling step: not reproducible, but the cell IDs were saved

Source (`20260129_entropy_cal_clean_v3_visium_spot_size.py`):

```python
max_num = min(data_df_filtered['group'].value_counts())
g2_sampled_idx = data_df_filtered[data_df_filtered['group'] == 'Group2'] \
                     .sample(max_num, replace=False).index.tolist()   # <- no random_state
selected_idx = data_df_filtered[data_df_filtered['group'] == 'Group1'].index.tolist() + g2_sampled_idx
pd.DataFrame(selected_idx, columns=['index']).to_csv(
    out_path / '20260126_max_group1_index_filtered_non_malig_equal.csv')
```

- **`.sample()` was given no `random_state`, and there is no global `np.random.seed` or
  `random.seed` either, so the result is not reproducible.** ✅ verified (read the whole
  script; no seed is set anywhere)
- **But the cell IDs that were kept were written to disk**, in exactly
  `20260126_max_group1_index_filtered_non_malig_equal.csv`.
  Format: two columns, `,index`, where `index` holds `Object_ID` as a UUID, e.g.
  `d2efe245-d546-48af-909a-969adf5efe5e`.

### Count reconciliation ✅ verified (awk over the 674 MB main table)

| | All cells | Passing the filter |
|---|---|---|
| Group1 | 1,224,858 | **660,310** |
| Group2 | 2,399,999 | 1,473,928 |

`max_num = min(660310, 1473928) = 660310`
-> `selected = all 660,310 of Group1 + 660,310 sampled from Group2 = 1,320,620`
-> **exactly the 1,320,620 rows of the index file.**

**Conclusion: read that index file and subset with `.loc[selected_idx]`; do not
resample.**
(Note: the `non_malig_num` used in the filter is not a column of the main table; the
script computed it on the fly from neighbor_counts. For the reconciliation it was
approximated as `non_malig_num = total_num - tumor_num`, and the numbers matched exactly,
which confirms that the 10 classes sum to total.)

---

## 3. The key trap: the cell type naming of `CODEX_adata_dict.npy` differs from the old chain

`D:\MyDownload\CODEX_adata_dict.npy` (7.9 GB) and
`20250929/h5ad/CODEX_adata_dict.npy` are the same thing (the submit scripts read the
latter).

✅ Verified by loading it:

- **31 samples** (`SP02, SP50, SP42, SP43` already excluded), matching the
  sample list of the entropy chain
- `obs` has only 5 columns: `r, l, r_cell, ct_agg, r50_adj_degree`
- `obsm`: `spatial` (coordinates, in μm)
- **`obsp` holds only `r50_adj_connectivities` / `r50_adj_distances`**
- `uns`: only `r50_adj_neighbors`

### 3.1 There is no `r27_5_connectivities`

`CODEX_adata_dict.npy` was built from **`{sample}_masked_counts_no_other.h5ad`** (without
the `_modified` suffix), whose obsp holds only the `r50_adj` graph, whereas the entropy
chain used **`_modified.h5ad`**, which does carry `r27_5`. ✅ verified by comparing the
obsp of the two

-> **the preprocessing script has to rebuild the 27.5 μm radius graph itself**, which is
possible because `obsm['spatial']` is present:
```python
sq.gr.spatial_neighbors(adata, coord_type='generic', radius=27.5, ...)
```
✅ Verified in `_modified.h5ad`: `r27_5_distances.data.max() == 27.500`, mean degree 5.05
(excluding self), coordinate ranges x:1.14-10985.00, y:15003.80-21654.30 -> **the
coordinate unit is indeed μm and 27.5 can be used directly**.

### 3.2 The `ct_agg` labels were renamed

A rename was applied when `CODEX_adata_dict.npy` was built (see the comments at
`codes/all/20260223_inspect_opc_h3.py:33-36`):

| Old naming (`_modified.h5ad` / the entropy chain) | New naming (`CODEX_adata_dict.npy` / `constants.py`) |
|---|---|
| `OPC_OC` | -> `OPC_OC_like` |
| `AC_MES` | -> `AC_MES_like` |
| `Oligo_pro` | -> `OPC_OC` |
| the other 7 classes | unchanged |

✅ Verified: the ct_agg of `_modified.h5ad` is
`['AC_MES','Astrocyte','B_cell','Cycling','Endothelial_cell','Myeloid','Neutrophil','OPC_OC','Oligo_pro','T_NK']`;
the ct_agg of `CODEX_adata_dict.npy` is
`['AC_MES_like','OPC_OC_like','OPC_OC','Cycling','Myeloid','T_NK','Endothelial_cell','Neutrophil','Astrocyte','B_cell']`

**The definition of malignant / non-malignant therefore has to be rewritten to match:**

```python
# old (as written in the entropy scripts)
malig_celltype = ['Cycling', 'OPC_OC', 'AC_MES']
# new (required when working from CODEX_adata_dict.npy)
malig_celltype = ['Cycling', 'OPC_OC_like', 'AC_MES_like']
```

**The good news: this is a one-to-one rename and the sets themselves did not change**
(both the old `Oligo_pro` and the new `OPC_OC` are non-malignant), so recomputing under
the new naming is **exactly equivalent** to the old chain.

#### Handled: the naming is settled once, in preprocessing

The submitted version uses the canonical naming throughout, and the development-era short
names are converted exactly once at the entry point (`Preprocessing.py`), after which
only one set of names exists anywhere (intermediate files, downstream scripts, figures).
The relevant definitions are collected in `codes/submit/constants.py`:

| Constant | Purpose |
|---|---|
| `old2new_celltype` | development short name -> canonical name, used by `Preprocessing.normalize_celltype()` |
| `celltype_display` | canonical name -> figure display name (`OPC_OC_like` -> `OPC/OC-like`, etc.) |
| `group_labels` | `{'Group1': 'PL', 'Group2': 'TL'}` |
| `group_colors` | `{'PL': '#990000', 'TL': '#2b8cbe'}` |
| `malig_celltype` / `non_malig_celltype` | the malignant / non-malignant split, so no script has to redefine it |

`normalize_celltype()` has an **ordering trap**, avoided in the implementation and
covered by a unit test: in the old naming `OPC_OC -> OPC_OC_like` while
`Oligo_pro -> OPC_OC`, so an item-by-item `Series.replace` may rewrite a value that has
just become `OPC_OC_like`. A single `map` is used instead, and the presence of
`Oligo_pro` decides which naming the input uses. The function is idempotent on new-naming
input and raises on any category outside `celltype_order`.

#### The two namings that coexist in the project (always establish which one you have)

| | Old naming | New naming |
|---|---|---|
| Defined in | `ct_agg_order` in `codes/cleaned_202604/codex/utils.py` | `celltype_order` in `codes/submit/constants.py` |
| Used by | `cleaned_202604/codex/*`, the old entropy scripts, `_modified.h5ad`, `{sample}_neighbor_counts.csv` | `CODEX_adata_dict.npy`, `codes/submit/*` (including the new scripts) |
| The three malignant states | `Cycling, OPC_OC, AC_MES` | `Cycling, OPC_OC_like, AC_MES_like` |
| Oligodendrocyte precursor | `Oligo_pro` | `OPC_OC` -- **the same name as the old malignant state but a completely different meaning** |

> The name `OPC_OC` does **not** refer to the same population in the two systems, and
> this is the easiest thing to get wrong. Always establish which naming the data at hand
> uses before copying code between scripts.

#### The other annotation: the 17-class `cell_type_2`

`codex/00_update_adata_obs.py` takes the fine-grained `cell_type` from `cell_type_csv`,
substitutes it only for cells with `ct_agg` in `{T_NK, Myeloid}`, and merges `Treg_*` into
`Treg`, giving `cell_type_2` (17 classes; the order is `celltype_order2` in
`codex/utils.py`):

```
Cycling, OPC_OC, AC_MES,
Dendritic_cell, Myeloid_CD68, Myeloid_CD163, Myeloid_M2_like, Myeloid_other,
Neutrophil, Tcell_CD4, Tcell_CD8, Treg, NK,
B_cell, Endothelial_cell, Astrocyte, Oligo_pro
```

It lives in `h5ad/adata_obs_dict.npy` (410 MB), and the downstream usage is
`adata.obs = adata_obs_dict[sample_name]` -- because the obs inside the h5ad is stale and
this npy is the authority on cell type.
**Note that it uses the old naming**, and **the entropy analysis does not use it** (entropy
works on the 10-class `ct_agg`). If entropy by fine-grained type is ever wanted, the data
is here.

#### Timeline (confirming that a later edit does not affect these results)

`cell_type_csv/*.csv` was tidied up in 2026-09 (the header was unified to
`Object_ID, cell_type, ct_agg, suffix` and the redundant columns `ct_agg_prev`,
`ct_agg_0102_bak`, `ct_agg_0209`, `ct_agg_w_umap_pred`, `Object_ID.1` and so on were
removed; `SP50` and `SP02` were left alone because they lack the `suffix` column).
`CODEX_adata_dict.npy` was generated on **2026-03-03**, before that tidy-up, and the
tidy-up **did not touch the `ct_agg` column itself**. Together with section 3.5(d) below,
where the `ct_agg` of adata_dict was verified cell by cell against the existing
`meta_data` after mapping, **these results are unaffected**.

### 3.3 Loading that npy needs a monkey patch (a fresh environment does not help)

✅ Verified: a plain `np.load(..., allow_pickle=True).item()` raises
`KeyError: '_adata_ref'` (at `anndata/_core/file_backing.py:45`).
The reason: the file was pickled by a **much older** anndata version, and
`AnnDataFileManager.__setstate__` in the current anndata 0.10.9 insists that the state
contain `_adata_ref`.

✅ Verified: the state keys of `AnnDataFileManager` in the pickle are
`['_adata', '_file', '_filemode', '_filename']` -- a strong `_adata` reference, with no
`_adata_ref`.

**Tested on freshly created environments (neither version works):**

| Environment | anndata | Deserialises | Object usable | Error |
|---|---|---|---|---|
| `codex_old` | 0.9.2 | yes | no | `AttributeError: 'AnnData' object has no attribute '_n_obs'` |
| `codex_old8` | 0.8.0 | yes | no | as above |
| `base` + patch | 0.10.9 | yes | **fully usable** | -- |

Conclusion: 0.8/0.9 get past `AnnDataFileManager` (because they use `_adata`), but the
`AnnData` of those versions depends on the `_n_obs` cache attribute, which the pickle does
not contain -> **the original file was written by anndata 0.7.x or earlier**. Loading it
natively would mean going back to 0.7.x and dragging the whole scanpy/pandas/numpy stack
back with it, which is not worth it. **Use the base environment plus the patch below**
(verified to load all 31 samples with exact numbers):

```python
import weakref, anndata
from anndata._core.file_backing import AnnDataFileManager

def _fm_setstate(self, state):                      # accept either _adata or _adata_ref
    self.__dict__ = dict(state)
    ref = state.get('_adata_ref', state.get('_adata', None))
    self.__dict__['_adata_ref'] = weakref.ref(ref) if ref is not None else (lambda: None)
AnnDataFileManager.__setstate__ = _fm_setstate

def _ad_setstate(self, state):                      # supply the missing _n_obs / _n_vars
    self.__dict__.update(state)
    if '_n_obs' not in self.__dict__:
        obs = self.__dict__.get('_obs', None)
        self.__dict__['_n_obs'] = len(obs) if obs is not None else 0
    if '_n_vars' not in self.__dict__:
        var = self.__dict__.get('_var', None)
        self.__dict__['_n_vars'] = len(var) if var is not None else 0
anndata.AnnData.__setstate__ = _ad_setstate
```
> This is wrapped in `_patch_anndata_compat()` in `codes/submit/Preprocessing.py`.
> The two new environments `codex_old` / `codex_old8` can be deleted, though keeping them
> does no harm.
> Note: the existing `Spatial_clustering.py` / `Proportion.py` /
> `Interaction_inference.py` / `Show_dotplot.py` under `codes/submit/` all `np.load` this
> file directly, so **they should hit the same error in the current environment** (not
> verified by running each of them).
> Also, `Spatial_clustering.py` requires `adata.obsp['r40_distances']` to exist, but the
> adata_dict was measured to hold only `r50_adj`, so that script reaches its
> `raise RuntimeError`. Both points are worth confirming.

Current environment versions: anndata 0.10.9 / scanpy 1.12 / numpy 1.26.4 / pandas 3.0.5
/ scipy 1.13.1 / squidpy 1.2.2

---

## 3.5 Reproducibility checks (measured 2026-09-28, all passed)

### (a) The Object_IDs match exactly, so a direct subset works ✅

| Check | Result |
|---|---|
| Total cells in adata_dict | 7,092,596 |
| Are the `obs_names` of adata_dict globally unique | yes (no duplicates across samples) |
| Sampled cell IDs found in adata_dict | **1,320,620 / 1,320,620 = 100.0000%** |
| IDs from `meta_data.csv` found in adata_dict | **3,624,857 / 3,624,857 = 100%** |

**Conclusion: `.loc[sampled_ids]` can be used directly, with no ID alignment needed.**

### (b) Why meta_data has 3.62 M rows while adata_dict has 7.09 M cells

✅ Verified: `max(total_num) = 20` in `all_samples_meta_data.csv`, and **100% of its rows
satisfy `total_num<=20`** -> **that file is in fact already the filtered data** (the file
name has no `_filtered` but the content does), which is corroborated by its md5 being
identical to `all_samples_meta_data_filtered.csv`.

Meanwhile the row counts of `{sample}_neighbor_counts.csv` agree **exactly** with
adata_dict (SP11: 25619 / SP31: 497190) -> the cell set used to build the graphs is
the same as adata_dict's, and no cell was lost.

**So the complete filter is equivalent to:
`tumor_num>=2 & non_malig_num>=2 & total_num<=20`** (the notna condition on `entropy` is
implied by the `>=2` conditions and does not need to be written separately).

### (c) The sklearn-rebuilt 27.5 μm graph is edge for edge identical to the original ✅

```python
radius_neighbors_graph(coords, radius=27.5, mode='connectivity', include_self=False)
```
Measured on SP11: `sklearn_nnz = orig_nnz = 129,490`, and **`diff_nnz = 0`, i.e.
identical**.
-> Preprocessing **does not depend on squidpy**; sklearn is enough, and
`include_self=True` corresponds exactly to the original `setdiag(1)`.

### (d) End-to-end numerical reproduction ✅ (the most important check)

Rebuilding the graph from adata_dict and recomputing, compared cell by cell against the
existing `meta_data.csv`:

| Sample | Cells compared | tumor_num | total_num | entropy_tumor | non_malig_entropy | ct_agg (after mapping) |
|---|---|---|---|---|---|---|
| SP11 | 14,567 | 0 mismatches | 0 | max\|diff\| = **2.22e-16** | **2.22e-16** | 0 |
| SP05 | 25,129 | 0 mismatches | 0 | max\|diff\| = **2.22e-16** | **2.22e-16** | 0 |

The post-filter cell counts also match exactly: SP11 recomputed 7,111 vs existing
7,111; SP05 recomputed 12,591 vs existing 12,591.

**2.22e-16 is float64 machine precision, i.e. identical. The whole chain reproduces
exactly.**

---

## 4. The two-step reorganisation (implemented as described)

### Step 1: `codes/submit/Preprocessing.py` ✅ written

Input: `20250929/h5ad/CODEX_adata_dict.npy`
Output directory: `codes/submit/data/`

| File | Contents |
|---|---|
| `entropy_cell_level.parquet` | the full cell-level main table (~7.09 M rows), so a filter change later can start from here |
| `neighbor_counts_r27_5.parquet` | the full 10-class neighbourhood counts |
| **`entropy_cell_level_sampled.csv`** | **the filtered and down-sampled table (~1.32 M rows); the single downstream entry point** |
| `sampled_cell_ids.csv` | the sampled cell IDs copied from the existing results; **no resampling** |

Points already handled in the script:
- `_patch_anndata_compat()` for the load compatibility (section 3.3)
- `radius_neighbors_graph(radius=27.5, include_self=True)` to rebuild the graph, where
  `include_self=True` corresponds to the original `setdiag(1)` under
  `mask_diagonal=False` (sections 3.1 / 3.5c)
- `MALIG_CELLTYPE = ['Cycling','OPC_OC_like','AC_MES_like']`, the new naming (section
  3.2), with a `NEW2OLD_CT` mapping kept for reference
- the filter thresholds raised to constants
  `TUMOR_NUM_LIMIT=2 / NON_MALIG_NUM_LIMIT=2 / TOTAL_NUM_LIMIT=20`
- assertions that Object_ID is unique, that the two tables are in the same row order, and
  that every sampled ID can be indexed

> The full-size layer uses parquet rather than csv: 7.09 M rows x 16 columns is about
> 1.3 GB as csv, and parquet is much better on both size and read/write time. The
> downstream entry layer stays csv so that R can read it directly.

### Step 2: `Entropy.py` + `Entropy.R` -- compute and visualise from the intermediate files

Both take `data/entropy_cell_level_sampled.csv` as their single entry point:
1. `Entropy.py`: group, then `gaussian_kde` with `nbins=20`, then export
   `grid_data_*_{Group1,Group2,diff}.csv`
2. `Entropy.R`: switched to an absolute path (or `here::here()`), reads the grid CSVs and
   draws the three figures, with the hard-coded colour breaks of the diff figure replaced
   by values computed from the data (see item 2 of section 5)

### Step 3: end-to-end rerun ✅ (actually run 2026-09-28)

`Preprocessing.py` -> `Entropy.py` -> `Entropy.R` all ran through, and were reconciled
item by item against the original outputs:

**The cell counts match exactly** (against the numbers counted independently from the old
chain)

| | Rerun | Old chain |
|---|---|---|
| All cells | 7,092,596 | -- |
| Passing the filter | 2,134,238 (G1 660,310 / G2 1,473,928) | 660,310 / 1,473,928 ✅ |
| Downstream entry table | 1,320,620 (660,310 per group) | 1,320,620 ✅ |
| Sampled IDs not found | 0 | -- |

**KDE grid values**: grid coordinates `max|diff| = 0` (identical), density
`max|diff| ~= 1e-14` (float64 machine precision).
The Group1 density peak `14.6187` matches the first row of the original grid csv,
`14.618703823238665`; the diff range `[-10.0534, 1.4671]` is exactly the `-10.05` and
`1.47` hard-coded in the original R script.

**Pixel-level figure comparison**

The figure naming has been changed to PL / TL (see "Figure naming" below), so the file
names changed with it:

| New figure | Original | Differing pixels, whole image | **Differing pixels inside the plot area** |
|---|---|---|---|
| `density_r.PL.auto.png` | `density_r.Group1.auto.png` | 35,440 (1.26%) | **0** |
| `density_r.PL.max.png` | `density_r.Group1.max.png` | 35,440 (1.26%) | **0** |
| `density_r.TL.auto.png` | `density_r.Group2.auto.png` | 35,386 (1.26%) | **0** |
| `density_r.TL.max.png` | `density_r.Group2.max.png` | 35,386 (1.26%) | **0** |
| `density_r_diff.png` | same name | -- | colour breaks are now dynamic, see below |

> **The plot area of all four figures (rows 100-1450, columns 100-1800) differs by zero
> pixels**, i.e. the heatmap and the contours are unchanged down to the pixel. The 1.26%
> difference over the whole image falls entirely on three pieces of text, each checked:
> columns 27-76 are the y axis label (17,423 px), rows 26-86 the title (5,131 px), and
> rows 1483-1532 the x axis label (12,886 px). The changes are the title
> `Group1/Group2 -> PL/TL` and the axis labels switching to the canonical names; the plot
> area's dimensions are unaffected by text length.
>
> In an earlier round of comparison, before the switch to PL/TL, the two Group1 figures
> differed from the originals by **zero pixels over the whole image**, and the 853 px of
> difference in Group2 fell entirely in rows 27-73, columns 910-941 -- a single 31 px wide
> character in the title. The cause is that the original script,
> `20260129_show_contour_group.R:68`, hard-coded `title = "Group1"`, so the Group2 figure
> was also titled "Group1". **That is a fix, not a deviation.**

> **The diff figure's** pixel differences come from the change in colour mapping (the
> original hard-coded the breaks, the new script generates them from the relative
> proportions of the data extremes). The visual structure is the same, and the title was
> corrected from the erroneous "Group1" to "PL - TL".

### Figure naming

The data layer (intermediate files, script variables, the `group` column of the CSVs)
always uses `Group1` / `Group2`; **figures (titles and file names) always use PL / TL**,
with the mapping collected in `group_labels` of `constants.py`:

| Data layer | Figures | n | Colour |
|---|---|---|---|
| `Group1` | **PL** | 14 | `#990000` |
| `Group2` | **TL** | 17 | `#2b8cbe` |

R cannot read `constants.py`, so `Entropy.R` carries an explicit copy of the same
`group_labels`; change both together.

The axes also switched to canonical display names (`axis_labels` in `Entropy.R`):

| Column | Axis label |
|---|---|
| `entropy_tumor` | Malignant-state entropy |
| `non_malig_entropy` | Non-malignant cell-type entropy |

---

## 5. Other pre-existing issues worth noting

1. **The `base_path` of `Entropy.R` is a bare relative path**,
   `'all_samples_full_27_5/'`, whose full path is
   `D:\local_storage\Work\SingleCell\Multiplex-IHC\spatial_analysis\20250929\20251210_dt\try_adjust_tumor\interacting_pairs_data\20260106_show_niche\all_samples_full_27_5\`
   so you have to `setwd()` to its parent `20260106_show_niche` for it to run.

2. **The colour breaks in `20260129_show_contour_diff.R` are hard-coded**:
   `values = scales::rescale(c(-10.05, -7, -4.5, -2.5, -0.8, 0, 0.4, 0.75, 1.1, 1.47))`
   and `breaks = c(-10, -7.5, -5, -2.5, 0, 1.47)`.
   Those numbers are tied to the density range of this particular sampling, and a
   resample would shift them out of place.
   -> The new `Entropy.R` generates them from the relative proportions of `min/max`
   (keeping the visual structure of the original).
   Note that the `colours` and `values` of `scale_fill_gradientn` **must be the same
   length** (10 here); unequal lengths shift the colour mapping and tint the whole figure.

3. **The two groups' KDE grid coordinates are in fact different**: the original
   implementation laid each group's grid over its own `x.min()/x.max()`, subtracted by
   array position, and used Group1's coordinates. `Entropy.py` preserves that behaviour
   (matching the original result) but prints a WARN. For strict alignment, the two groups
   should share one grid.

4. **`x = np.where(x > 60, 60, x)`** (in the L2 script, before the KDE) is meaningless for
   entropy -- entropy has an upper bound of about 1.1. The line was copied from elsewhere
   and is dead code; it does not affect the result.

5. **`if os.getenv("data_dict") is None: data_dict = {}`** does not do what it looks like
   (`os.getenv` reads an environment variable and is always None here), so in practice it
   resets every time; it does not affect the result.

6. In `Entropy.R`, `df = df[df$non_malig_entropy < 1.2,]` and
   `coord_cartesian(ylim = c(0, 1.2))` are plot clipping, not an analysis filter.

7. **None of the scripts involved here contains a boxplot**, so no outlier-related change
   was needed.

---

## 6. Methods reference

Two Nature-style Methods paragraphs from the same body of work are already final, and the
entropy section can follow their style (radii written as `40 μm`, with a space between
number and unit; `Mann–Whitney` with an en dash):

**Spatial coherence analysis.** A connectivity graph was constructed by connecting cells
located within 40 μm of each other. Based on this graph, the spatial coherence of cell type A
was defined as the number of edges connecting two cells of type A divided by the total number
of edges involving A. To account for differences in cell density across samples, we also
calculated a density-adjusted spatial coherence. For each sample, 100 index cells were randomly
selected with a minimum spacing of 40 μm, and the number of cells within 40 μm of each index
cell was counted; the mean of these counts was taken as the local cell density of the sample.
Density-adjusted spatial coherence was defined as spatial coherence divided by local cell
density. Differences between groups were assessed using a two-sided Mann–Whitney U test.

**Nearest-neighbour distance analysis.** For each sample, the average nearest-neighbour
distance from cell type A to cell type B was calculated as the mean Euclidean distance between
each A cell and its nearest B cell (excluding the cell itself when A = B), using a k-d tree built
on cell centroid coordinates. ... Distances were log2-transformed and compared between PL and TL
samples for each source–target pair using a two-sided Mann–Whitney U test, with individual
samples as the unit of analysis.

> **PL = Group1, TL = Group2** ✅ confirmed: several scripts state
> `GROUP_LABELS = {"Group1": "PL", "Group2": "TL"}`
> (e.g. `cleaned_202604/202604_rebuttal_related/20260504_t_to_tumor_cont.py:26` and
> `ccc_summary_bbx_grp_comp_map_tum_only.py:18`), and the old entropy script
> `20260129_..._visium_spot_size.py:147` agrees with its
> `ax[1,0].set_ylabel('PL'); ax[2,0].set_ylabel('TL')`.
> Group colours: `{"PL": "#990000", "TL": "#2b8cbe"}` (PL red, TL blue).
> This document and the new scripts use Group1 / Group2 throughout; **figures shown
> externally need PL / TL**.

Points to state when writing the entropy section (matching the actual implementation):
1. Neighbourhood definition: radius 27.5 μm (the diameter of a Visium spot), **including
   the index cell itself**
2. Two entropies: the Shannon entropy across the three malignant states, and the Shannon
   entropy across the seven non-malignant classes
3. Inclusion criteria: at least 2 malignant and at least 2 non-malignant cells in the
   neighbourhood, and at most 20 cells in total
4. The two groups have unequal cell counts, so the larger group (Group2) was
   down-sampled without replacement to match Group1
5. The 2-D distribution is a Gaussian kernel density estimate evaluated on a 20 x 20 grid
