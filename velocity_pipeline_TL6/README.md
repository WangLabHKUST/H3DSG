# velocity_pipeline_TL6

RNA velocity of the malignant cells of patient SP73, sampled at five sites (`renamed_id` TL_6 to TL_10). The pipeline
goes from the BAM files through spliced and unspliced counts and scVelo to a UMAP with velocity stream lines, coloured by
`side_name`, `sec_name`, local dorsal fraction (`nn10_prop_Dorsal_umap`) and `root_cells`.

```bash
conda env create -f environment.yml
conda activate velo
python run_pipeline.py                          # steps 1-6; a step whose output exists is skipped
python run_pipeline.py --from 2                 # the looms exist, start at step 2
python run_pipeline.py --from 5 --to 6          # only UMAP and figures
python run_pipeline.py --to 1 --dry-run         # only print the velocyto commands
python run_pipeline.py --force --from 6         # recompute a step and the ones after it
python run_pipeline.py --stems ZDK-inc-edge ZDK-exc-edge   # also the set without the edge sample
```

Each step also runs alone: `python steps/s5_umap.py [--force] [--stems ...]`. All paths are relative to this folder.
velocyto does not build on Windows; steps 2 to 6 do not need it.

## Inputs

None of them is committed, see `.gitignore`. The pipeline does not change them. Paths are set in `config.py` and `samples.tsv`.

| Path | Content | Source |
|---|---|---|
| `data/bam/<sample>/possorted_genome_bam.bam` (and `.bai`) | Cell Ranger BAM | raw sequencing data of the study |
| `../Single_Cell_RNAseq/inputs/raw_counts/<sample>/barcodes.tsv.gz` | filtered barcodes. They only select the cells; spliced and unspliced counts all come from the BAM | Cell Ranger, see `../Single_Cell_RNAseq/inputs/README.md` |
| `data/ref/gencode.v44.annotation.gtf` | gene annotation | GENCODE v44 |
| `data/seurat_tables/meta_subset_excl_Immune_Stromal_OPC-like.csv` | cell metadata (`side_name`, `sec_name`, `State2025`, ...) | exported from Seurat, see below |
| `data/seurat_tables/coords_pca_excl_Immune_Stromal_OPC-like.csv` | `PC_1..PC_30`, the `X_pca` used for velocity | exported from Seurat, see below |
| `data/seurat_tables/coords_umap_subset_excl_Immune_Stromal_OPC-like.csv` | Seurat UMAP. Stored as `X_umap` in the h5ad, not used by the final figures (can be `None`) | exported from Seurat, see below |

`samples.tsv` links `seurat_id` (`orig.ident` in Seurat, TL_6_1 to TL_6_5) to `sample` (the sequencing file ID);
`../Single_Cell_RNAseq/sample_metadata.tsv` links both to `renamed_id`. Samples with `margin=1` are left out of the
`exc-edge` set.

The three Seurat tables hold the 8,491 malignant cells of the five samples, the same cells as
`../Single_Cell_RNAseq/pipeline/08_trajectory_TL6/`, with the corrected PCs of the integrated object. Their `cell_id` columns
must be in the same order (step 2 asserts this).

## Steps and outputs (`results/`, not committed)

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `s1_velocyto.py` | `velocyto run` (default logic, no mask) | `s1_loom/<sample>/<sample>.loom` |
| 2 | `s2_prepare_inputs.py` | adds the `file_id` key (sample name) to the Seurat tables, and alias columns | `s2_inputs/{meta,harmony,umap}_df_combined.csv` |
| 3 | `s3_build_h5ad.py` | loom + tables -> h5ad, merged per analysis set | `s3_h5ad/<sample>.h5ad`, `s3_h5ad/<stem>_combined.h5ad` |
| 4 | `s4_velocity.py` | scVelo: moments, deterministic velocity, velocity graph, confidence, pseudotime (root and end points) | `s4_velocity/<stem>_velocity.h5ad` |
| 5 | `s5_umap.py` | UMAP recomputed on the first 20 PCs of `X_pca` (correlation, nn = 100, min_dist = 0.6, seed 0) | `s5_embedding/<stem>_umap_npc20_correlation_nn100_md0.6.csv` |
| 6 | `s6_figures.py` | stream lines on that UMAP, four colourings; local dorsal fraction | `s6_figures/20260916_flow_on_umap_<tag>_<color>_<stem>.png`, `s6_figures/20260916_nn10_prop_side_<tag>_<stem>.csv` |

`results/logs/s<N>.log` holds the output of each step, and `results/run_info.json` the package versions and parameters.

## Notes

- **Step 1 is heavy.** velocyto needs a BAM sorted by cell barcode. Without `cellsorted_<bam name>` next to the BAM it
  writes one there (about 25 GB per sample here, several hours). A sample that already has a loom is skipped.
- **The UMAP is not the Seurat one.** Step 5 recomputes it from `X_pca`. The velocity graph depends only on the PCA space;
  another UMAP changes only the position of the cells and the projection of the stream lines.
- **The local dorsal fraction depends on the UMAP parameters.** It counts the Dorsal cells among the 10 nearest
  neighbours in the 2D coordinates of step 5, so it changes with `config.UMAP`. It is descriptive only.
- **Labels.** The figures of step 6 carry the original column values (`edge`, `C4-C4`) and titles.
- **Other data.** Change `samples.tsv` and the "Raw inputs" section of `config.py`. `exclude_margin` in `STEM_DEFS` decides
  which samples the `exc-edge` set leaves out.
