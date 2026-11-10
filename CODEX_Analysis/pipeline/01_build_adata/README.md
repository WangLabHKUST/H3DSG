# 01_build_adata — from raw input to `CODEX_adata_dict.npy`

This layer records where the `.npy` came from.

```
QuPath measurements  raw_data/qupath_csv/{sample}.csv ──┐
cell type tables     raw_data/cell_type_csv/{sample}.csv ──┤
existing graphs      raw_data/h5ad/{sample}_full_raw.h5ad (.obsp) ──┘
        │
        ▼  00_create_adata.py        <- where ad.AnnData(...) is actually called
   {sample}_masked_counts.h5ad
        │
        ▼  00_modify_adata.py        <- filtering + building the spatial graphs
   {sample}_masked_counts_no_other.h5ad           drops ct_agg in {Other, Tumor_other}
   {sample}_masked_counts_no_other_modified.h5ad  builds the de / r40 / r50 / r27_5 graphs
        │
        ▼  00_update_adata_obs.py    <- adds cell_type and the 17-class cell_type_2
   raw_data/h5ad/adata_obs_dict.npy
        │
        ▼ (the final gathering step, not included in this directory)
   raw_data/CODEX_adata_dict.npy     <- the file in raw_data/ of this repository

({sample} in the file names above is the sample key the files were created with, which is
what `data_loading.to_raw_key()` returns for an SP ID; the scripts work on SP IDs and
translate where they build a file name.)
```

## The first three scripts are archived verbatim and cannot be run

They are notebook-style working scripts, kept as they were so that this leg is
documented — not so that it can be rerun:

- **One sample per run**: `sample_name = 'SP08'` is written straight into the file
  (under a comment reading `################# input here ##############`); changing
  sample means editing that line and running again. `00_update_adata_obs.py` is the only
  one with a loop.
- **The inputs are expected under `raw_data/`**, at the locations `path_constants.py`
  names: `qupath_csv/` (22 GB), `h5ad/` (357 GB, including the 410 MB
  `adata_obs_dict.npy`), `cell_type_csv/` (2.8 GB), `cell_loc/` (0.5 GB),
  `detailed_types/` (0.7 GB) and `common_data/`. Only `common_data/` (the two small tables
  `pipeline/00_sample_metadata/` reads) is shipped there; the others are too large to
  copy and have to be placed or linked there to run these scripts.
- **They depend on two modules this repository does not contain**:
  `00_create_adata.py` imports `all_proteins` and `binary_features` from
  `codes/cleaned_202604/codex/utils.py`, and `00_modify_adata.py` imports
  `save_graph_obsp` and `clean_uns_dict` from `graph_interaction_plotting_utils.py`.
- `00_create_adata.py` also imports `adata_obs_dict_path` from `path_constants`, and the
  copy of `path_constants.py` here **does not define that name** — a consequence of
  copying verbatim, which was deliberate.
- In `00_modify_adata.py`, `adata = adata[~adata.obs['ct_agg'].isin(...)]` yields a view
  and is written out with `write_h5ad` immediately afterwards, which newer versions of
  anndata warn about.

**The `.npy` already exists (in `../../raw_data/` of the archive version), so normally
there is nothing to rerun in this layer.**

Beyond that, `00_create_adata.py` can **no longer run even on the original analysis
machine**: the `{sample}_full_raw.h5ad` files it reads no longer exist (verified: zero
of them), having come from the earlier `codes/all/20251024_write_adata.py`.

## `04_export_cell_type_2.py` — the one runnable script here

Unlike the three above, this one was written for this repository and does run on the
original analysis machine. Its job is narrow: take the `cell_type_2` column (the
17-class fine-grained annotation) for the 31 samples out of
`raw_data/h5ad/adata_obs_dict.npy`, convert it to the manuscript naming, and store it as
`../../raw_data/cell_type_2.parquet` (257 MB).

Why it is needed: the `obs` of `CODEX_adata_dict.npy` carries only `ct_agg` (10
classes), whereas `pipeline/08_min_distance/` works on the 17 classes. A collaborator
receives its output and does not need to run it.

The naming conversion uses `constants.old2new_celltype`, and is **a three-way exchange,
not a simple rename**:

```
old OPC_OC     (malignant state) -> new OPC_OC_like
old AC_MES     (malignant state) -> new AC_MES_like
old Oligo_pro  (normal cell)     -> new OPC_OC          <- takes over the name OPC_OC
```

It therefore has to be a single atomic `.map()`. A chained `.replace()` would first turn
`Oligo_pro` into `OPC_OC` and then let the first rule turn that into `OPC_OC_like`,
merging normal oligodendrocyte precursors into a malignant state. The script asserts the
cell ordering (the `obs_names` of all 31 samples are compared against the `.npy` one by
one).

## Two things copied over verbatim that will trip you up while reading

**Channel positivity calling.** In `00_create_adata.py`:

```python
raw_data = raw_df[all_proteins].values.astype(float)      # measured intensity, 23 channels
binary_data = raw_df[binary_features].values.astype(int)  # QuPath's positive/negative call
filtered_data = raw_data * binary_data                    # the X of the AnnData
```

So **X is intensity × positivity call**, and a negative cell is exactly 0 in that
channel. `frac_expressing` in `pipeline/05_marker_expression/` (the fraction with X > 0)
is therefore equivalent to the fraction of positive cells. `binary_features` is
`all_proteins` minus HLA_A/HLA_DR, each with an underscore suffix appended, plus
`MCH_I_APC` and `MCH_II_APC`.

**Samples are excluded in two stages.** In
`codes/cleaned_202604/codex/utils.py`:

```python
sample_lst_full = sample_info.index.tolist()                                    # 35
sample_lst_33 = [f for f in sample_lst_full if f not in ['SP02','SP50']]  # low-quality
sample_lst_31 = [f for f in sample_lst_33 if f not in ['SP42','SP43']]    # recurrent
```

`00_update_adata_obs.py` repeats `['SP02', 'SP50']` on its own. 31 samples end up
in the `.npy` (verified). Those four samples and the reason each was dropped are
recorded in the `in_analysis` and `exclusion_reason` columns of
`../../sample_metadata.tsv`.

## The 17-class subdivision

`00_update_adata_obs.py` also writes `cell_type_2`, which unfolds the `T_NK` and
`Myeloid` classes of `ct_agg` into Dendritic_cell / Myeloid_CD68 / Myeloid_CD163 /
Myeloid_M2_like / Myeloid_other / Tcell_CD4 / Tcell_CD8 / Treg / NK and so on, 17 classes
in total (`celltype_order2` in the script). Steps 02 to 07 of this repository all work on
the 10-class `ct_agg`; only `pipeline/08_min_distance/` uses the 17 classes, through
`../../raw_data/cell_type_2.parquet` as exported by `04_export_cell_type_2.py`.

Note that the naming inside `adata_obs_dict.npy` is the development-era naming
(`OPC_OC` meaning the malignant state, `Oligo_pro` meaning normal oligodendrocyte
precursors), which is not the same vocabulary as `constants.celltype_order`; see
"Cell type naming" in the top-level README. The parquet export has already been
converted to the manuscript naming.
