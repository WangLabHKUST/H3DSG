# Pipeline inputs

Everything the pipeline reads that it does not compute itself. Two directories:

- **`inputs/`** — raw counts, sample tables, manual annotation decisions, reference gene
  sets and frozen inputs. The small files are committed. Not committed (see
  `../.gitignore`): `raw_counts/`, the `.rds` files, `meta/` and `sample_metadata.tsv`
  (sample metadata; the main mapping carries clinical fields), and the read files of the
  samples outside the study in `mut_reads/`.
- **`../objects/`** — local copies of the Seurat objects the analysis steps start from.
  Not committed. `objects/MANIFEST.tsv` is, and lists every file with its size, md5 and
  the path it was copied from.

## inputs/

| Path | Content |
|---|---|
| `raw_counts/<file_id>/` | not committed. 22 10x filtered matrices (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`), one directory per `file_id`: the 19 samples of the study plus `DHY-1`, `DHY-2`, `DHY-4`, which are processed with them but are not part of the 19 |
| `sample_metadata.tsv` | not committed. A copy of `../sample_metadata.tsv`, the 19 samples, kept so that `inputs/` is complete on its own |
| `meta/orig_ident_group_mapping_multi_ID.csv` | not committed. The main mapping, 22 rows, keyed by `orig.ident` (`SP24`, `ZDK-1.filtered_feature_bc_matrix`, …), with the sample locations and the patients' clinical fields |
| `annotations/` | manual decisions as CSV, see below |
| `refsets/` | reference gene sets: `K27scDSG.cellular.states.6sates.202502.xlsx` (used to name the metaprograms), `filbin_gene_set.xlsx` (Filbin cell states), `zhang2021_table_s3.xlsx` (Zhang 2021 T cell metaclusters) |
| `mut_reads/batch_reads/` | per-`file_id` tables of reads covering the H3 K27 position, 24 files. The 18 of the study's samples are identical to `../pipeline/07_k27m_reads/K27M_reads/` and are committed; the six of `DHY-1`, `DHY-2`, `DHY-4`, `MTSC08`, `WX-2`, `WX-7` are not |
| `hox/` | `loc_gene.csv` (gene to anatomical level) and two reference object lists for thalamic and pontine tumours (`location_obj_list2.rds`, 34 genes; `location_objects_full.rds`, all genes, ~700 MB) for the HOX comparison; their spinal entries are dropped |
| `frozen/` | inputs that cannot be recomputed, see below |
| `hg38.cellranger.genes.position.txt` | gene positions for inferCNV |

### annotations/

Cluster-to-label mappings live here, not in the scripts: changing a mapping means changing
the CSV. A cluster number is only valid for the clustering that produced it, so before a
mapping is applied the new clustering is compared with the one the mapping was written
for: identical numbering is applied directly, a renumbering is applied through the old
numbers, an ARI of 0.9 or more inherits the old label cell by cell, and anything lower
stops for manual review.

| File | Content |
|---|---|
| `major_only_PL_6_4_integrated_snn_res.1_anno.csv` | the 33 clusters of the CCA integration at resolution 1 to 8 major cell types; cluster 31 is `myeloid_mixed_T` and is dropped from the final object |
| `tumor_malignant_filter_overrides.csv` | manual overrides of the per-sample malignant cluster calls: `TL_2` cluster 11 to `malig`, `TL_5` cluster 10 to `exclude` |
| `tumor_mp_state_names_plateau_v1_all19.csv` | metaprograms M1 to M6 to state names (G2M, OC-like, S, AC-like, OPC-like, MES-like); the scoring step checks that each metaprogram's best match in the reference agrees with this table and stops otherwise |
| `tcell_nf4000_exclude_res0.3.csv` | T/NK clusters at resolution 0.3 that are removed (4,438 to 4,201 cells) |
| `tcell_nf4000_res0.25_celltype.csv` | the 10 T/NK clusters at resolution 0.25 to 9 states |
| `tcell_nf4000_ifn_subcluster.csv` | the interferon subcluster (`FindSubCluster`, resolution 0.05) to labels |
| `tcell_family_markers.csv` | marker panels per T cell family, for the marker check |
| `zhang2021_metacluster_family.csv` | Zhang 2021 metaclusters to T cell families, for the reference mapping |
| `myeloid_nf2000_exclude_res0.3_orig.csv` | the tumour-mixed myeloid cluster removed in the first resolution 0.3 clustering |
| `myeloid_nf2000_state_fine.csv` | second-pass resolution 0.3 clusters, with cluster 9 split at resolution 0.1, to the 12 fine myeloid states; clusters 8 and 11 are `artifact` |
| `myeloid_nf2000_state_coarse.csv` | the same clusters to the 8 coarse myeloid states |

### frozen/

| File | Rows | Why it is frozen |
|---|---|---|
| `cleaned_0806_cells.txt` | 125,105 cells | the cell list after the 2026-08-06 cleaning, the starting set of the myeloid analysis and the `TL_6` trajectory; the rule that produced it could not be reconstructed. Recorded as `in_cleaned_0806` in the final object |
| `included_in_prop.csv` | 132,610 cells | the flag that downsamples the five `TL_6_*` samples of one patient in the malignant state composition (1,803 of 9,017 kept). The original draw had no seed. 81,304 malignant cells have a value, 74,090 of them `TRUE` |
| `tumor_state_0811.csv` | 132,610 cells | the earlier version of the 10 malignant state columns, stored next to the rescored ones in the final object |
| `tnk_manual_anno.csv` | 4,201 cells | the manual T/NK annotation (`cell`, `sid`, `state`, `state_pretty`) |

## objects/

| Path | Cells | What it is | Read by |
|---|---|---|---|
| `per_sample/<sid>_minf400.rds` | per sample | the 19 per-sample objects after QC, with their own clusters and UMAP | malignant cell filtering, metaprogram QC, inferCNV |
| `integrated/integrated_minf400_only_PL_6_4_cca_anno.rds` | 132,610 | CCA integration of the 19 samples with the major cell types | malignant state scoring, T/NK and myeloid states, `TL_6` trajectory |
| `integrated/integrated_minf400_only_PL_6_4_cca_annotated.rds` | 132,343 | the assembled object: `myeloid_mixed_T` removed, all state columns merged in | HOX comparison, inferCNV |
| `subsets/myeloid_integrated_anno_cleanedpass1_nf2000.rds` | 39,573 | the annotated myeloid object | final myeloid panels |
| `subsets/integrated_minf400_only_PL_6_4_cca_TNK_integrated_nf4000_cleaned_anno.rds` | 4,201 | the annotated T/NK object | manual T/NK annotation |
| `released/scRNA_final.rds`, `scRNA_lymphocyte_final.rds`, `scRNA_myeloid_final.rds` | 132,343 / 4,201 / 39,573 | the three final objects, complete: all reductions (including `umap_pc30` and `umap_uncorrected`), graphs and the integrated `scale.data` | `pipeline/01_preprocessing/02_qc_overview.R`, `pipeline/04` – `08`, `pipeline/99_export_figure_data.R` |
| `released/scRNA_*_final_slim.rds` | same | the slimmed versions described in `../README.md`: `umap_pc30` renamed to `umap`, the default `umap` and `umap_uncorrected` dropped, no `scale.data`, notes in `misc(object)$reduction_notes` (all-cell object only). `scRNA_final_slim_xz.rds` is the all-cell one compressed with xz | not read by the pipeline |
| `released/scRNA_*_final_diet.rds`, `scRNA_*_final_diet_rna.rds` | same | `diet`: RNA counts and data, integrated data, `pca` and `umap`, no graphs; `diet_rna`: the same without the integrated assay | not read by the pipeline |
