# Figure-ready tables

Everything a figure needs, and nothing else. These tables are committed to the repository
so that a panel can be redrawn without downloading anything; 12 MB in total. None of them
is hand-edited.

All but two are written by `../pipeline/99_export_figure_data.R` from the three released
Seurat objects, and can be regenerated with it and diffed. The two others,
`sample_level/k27_alt_freq_by_sample.tsv` and `panel_data/cxcl14_dotplot_data.csv`, were
exported by the analysis code that drew their panels and are not yet produced by
`99_export_figure_data.R`; each carries one column set from an earlier version of the
analysis that the released objects do not hold (see the notes below), so they cannot be
regenerated from the objects alone.

Three kinds of table, which is the answer to "some panels need the object, some need
per-sample numbers":

- **`cell_level/`** — one row per cell. Replaces opening a Seurat object: it holds the
  annotation columns and the UMAP coordinates, which is everything a UMAP, a stacked bar
  or a per-cell histogram needs. 132,343 rows; gzipped, 2.5 MB for the annotation and
  2.5 MB for the coordinates.
- **`sample_level/`** — one row per sample, or per sample and category. This is what the
  comparison plots use: a boxplot of a state fraction in PL against TL has 19 points, one
  per sample, and reads 19 rows. Tens of kilobytes.
- **`panel_data/`** — tables for panels that need something neither of the above holds,
  in practice expression values. A dot plot of one gene across cell types needs that one
  gene, not the expression matrix, so the export writes out that gene's values alone.

File names here say what the table holds, never which figure uses it: panel numbers change
between revisions, and the panel index in `../figures/README.md` is the one place a figure
number appears.

Expression values are the only reason a panel would otherwise have to open an object.
Whenever a new panel needs a gene, add it to the gene list at the top of
`99_export_figure_data.R` rather than reaching for the `.rds` from the figure script.

## Keys

`cell_id` joins every table in `cell_level/`. It is the column name of the released
object, which carries a per-sample prefix; barcode formats differ between samples, so
`cell_id` is a string to join on and never to parse. `sid` joins to
`../sample_metadata.tsv`; `renamed_id` is the name to print in a figure.

`patient_id` holds the patient's Table S1 `Sample ID` (`SP..`). The UMAP panels of the major cell types and of the K27M reads are drawn on `umap_corr_1` and `umap_corr_2`.

## Files

| File | Rows | Columns |
|---|---|---|
| `cell_level/cell_annotation.tsv.gz` | 132,343 | `cell_id`, `sid`, `renamed_id`, `group`, `patient_id`, `nCount_RNA`, `nFeature_RNA`, `percent.mito`, `cell_type`, `malignant_state`, `tnk_state`, `myeloid_state_coarse`, `myeloid_state_fine`, `included_in_prop`, `substate_excluded_reason`, `depth_Alt`, `depth_Total`, `Alt_freq`, `is_alt`, `is_malignant` |
| `cell_level/umap_coordinates.tsv.gz` | 132,343 | `cell_id`, `umap_1`, `umap_2` (the 30-component UMAP the major cell types were annotated on, `umap_pc30` in the full object), `umap_corr_1`, `umap_corr_2` (`umap_pc20_correlation`) |
| `cell_level/malignant_state_scores.tsv.gz` | 81,304 | `cell_id`, `G2M`, `S`, `OPC_like`, `AC_like`, `MES_like`, `OC_like`, `state`, `state_score` |
| `sample_level/qc_summary.tsv` | 19 | `sid`, `renamed_id`, `group`, `patient_id`, `n_cells`, `median_nCount`, `median_nFeature`, `median_percent_mito` — the cohort as released, after QC; the pre-QC counts and the per-sample thresholds are not in the released objects |
| `sample_level/cell_type_composition.tsv` | 129 | `sid`, `renamed_id`, `group`, `patient_id`, `category`, `n`, `fraction` — one row per sample and cell type present in it |
| `sample_level/malignant_state_composition.tsv` | 113 | same shape over the six malignant states (G2M and S separate), plus `included_in_prop`; see the note below |
| `sample_level/tnk_state_composition.tsv` | 172 | same shape, the 11 T/NK states |
| `sample_level/myeloid_state_composition.tsv` | 151 | same shape, the 8 coarse myeloid states |
| `sample_level/myeloid_state_fine_composition.tsv` | 225 | same shape, the 12 fine myeloid states |
| `sample_level/k27_alt_freq_by_sample.tsv` | 19 | `sid`, `renamed_id`, `n_cells`, `n_reads_old`, `n_reads_new`, `n_alt_ge_thr`, `median_depth` — cells in the sample, cells with at least one read over the H3 K27 position (first and current read extraction), cells with `Alt_freq >= 0.4`, and the median read depth of the covered cells |
| `panel_data/panel_gene_expression_nonzero.tsv.gz` | 799,006 | `cell_id`, `gene`, `expr` — log-normalised RNA values for the genes in `PANEL_GENES`, non-zero entries only; 4.8 MB |
| `panel_data/cxcl14_dotplot_data.csv` | 36 | `version`, `state`, `Group`, `n_cells`, `avg.exp`, `pct.exp` — CXCL14 per cell state and PL / TL, as `DotPlot(scale = FALSE)` reports it: `avg.exp` is the mean of the de-logged expression, `pct.exp` the percentage of cells expressing it |
| `panel_data/marker_gene_expression_nonzero.tsv.gz` | 386,165 | `cell_key`, `gene`, `expr` |
| `panel_data/tnk_cells.tsv.gz` | 4,201 | `cell_id`, `sid`, `renamed_id`, `group`, `patient_id`, `tnk_state`, `umap_1`, `umap_2`, `in_pooled_bar`, `PDCD1_lognorm` |
| `panel_data/tnk_samples.tsv` | 19 | `sid`, `renamed_id`, `group`, `patient_id`, `n_tnk`, `n_malignant` |
| `panel_data/tnk_marker_dotplot.tsv` | 319 | `tnk_state`, `gene`, `gene_group`, `n_cells`, `avg_exp`, `pct_exp`, `avg_exp_scaled` |
| `panel_data/tnk_signature_scores.tsv` | 440 | `tnk_state`, `signature`, `n_cells`, `mean_score` |
| `panel_data/tnk_signature_genes.tsv` | 1,197 | `signature`, `gene_order`, `gene`, `comb_ES_rank` |
| `panel_data/malignant_state_heatmap_5000cells.tsv.gz` | 5,000 | `cell`, `draw_position`, `sid`, `renamed_id`, `state`, `ATAD2`, `CENPK`, `CENPM`, `CENPU`, `CLSPN`, `DEK`, `DNAJC9`, and 290 more |
| `panel_data/pseudotime_TL6_malignant_cells.tsv` | 8,491 | `cell_id`, `sid`, `renamed_id`, `patient_id`, `tumor_location`, `site`, `cluster_res0.1`, `Lineage1`, `Lineage2`, `Lineage3`, `pt_mean` |
| `panel_data/hox_dotplot_stats_by_sample_celltype.tsv` | 2,166 | `renamed_id`, `celltype`, `gene`, `n_cells`, `mean_log_expr`, `mean_delog_expr`, `pct_exp` |
| `panel_data/pax3_nkx61_dotplot_stats_by_sample_celltype.tsv` | 114 | `renamed_id`, `celltype`, `gene`, `n_cells`, `mean_log_expr`, `mean_delog_expr`, `pct_exp` |
| `panel_data/hox_gene_annotation.tsv` | 38 | `gene`, `location_group`, `order_all_samples`, `order_multisite` |
| `panel_data/sample_resection_order.tsv` | 19 | `renamed_id`, `patient_id`, `group`, `tumor_location`, `resection_order`, `site_number` |
| `metaprograms/program_similarity_matrix.tsv` | 107 | `program`, `PL_1_P6`, `PL_4_P6`, `PL_2_P3`, `PL_11_P8`, `PL_3_P7`, `PL_5_P7`, `PL_9_P5`, `TL_1_P7`, `TL_6_P7`, `TL_8_P3`, `TL_5_P5`, and 96 more |
| `metaprograms/program_metaprogram_assignment.tsv` | 107 | `program`, `sample`, `metaprogram`, `state`, `extraction_order`, `clustered_order`, `block_order`, `cor_S_score`, `cor_G2M_score` |
| `metaprograms/metaprogram_signature_genes.tsv` | 300 | `metaprogram`, `state`, `rank`, `gene` |
| `metaprograms/metaprogram_overlap_liu2022.tsv` | 36 | `metaprogram`, `state`, `liu_state`, `shared_genes`, `row_order`, `col_order` |

The per-sample cluster calls that choose the inferCNV reference are in `../pipeline/06_infercnv/malignant_cluster_calls.csv`, next to the script that reads them.

## Notes on four columns

`included_in_prop` is not redundant with `is_malignant`. 81,304 cells are malignant;
74,090 of them are `included_in_prop`, and the malignant state composition figures use
those. The flag downsamples the one patient with five samples: of the 9,017 malignant
cells of `TL_6` … `TL_10` (`sid` `TL_6_1` … `TL_6_5`), 1,803 are included; every malignant
cell of the other 14 samples is. Quoting 81,304 as the denominator of a composition figure
would be wrong.

In `sample_level/malignant_state_composition.tsv`, `n` and `fraction` are computed over
all malignant cells of a sample, and `included_in_prop` is `TRUE` when any of them is
included, which holds for all 19 samples. The composition the figures use is therefore
not this table: compute it from `cell_level/cell_annotation.tsv.gz`, keeping the cells
with `included_in_prop == TRUE`.

`k27_alt_freq_by_sample.tsv` has both `n_reads_old` and `n_reads_new`. The five `TL_6_*`
samples were missing from the first mutant read extraction and were added later; the two
columns are kept so the difference is visible. Only `n_reads_new` should be used.

`cxcl14_dotplot_data.csv` holds two versions of the cell states: `version = cell.type`
is the state assignment of the released objects, `version = cell.type_0811` an earlier
one that the released objects do not carry, kept for comparison with the earlier version
of the panel. In both, OPC-like and OC-like are merged into `OPC/OC`, AC-like and MES-like
into `AC/MES`, and G2M and S into `Cycling`.
