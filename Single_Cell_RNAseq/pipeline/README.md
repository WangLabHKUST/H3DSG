# Pipeline: raw counts to the released objects

Run the numbered directories in order. Every script declares its inputs and outputs as
path constants in the first 30 lines and is runnable on its own — nothing is `source()`d
from anywhere else, so a script can be read and checked without following a chain of
includes.

Paths are relative to the script's own directory, so run each script from there. Inputs
come from `../inputs/` (raw counts and small inputs) and `../objects/` (the deposited
Seurat objects; `../inputs/README.md` lists both). Outputs go to `../results/<step>/`,
which is not committed, so a re-run never overwrites the deposited object it is meant to
be compared with. The one exception is step 99, whose output is the committed `../data/`.

In the table, `objects/released/scRNA_final.rds` is written `scRNA_final.rds`.

| Step | Script | Input | Output |
|---|---|---|---|
| 01 | `01_preprocessing/01_preprocess_per_sample.R` | `inputs/raw_counts/<file_id>/` | `results/per_sample/<sid>_minf400.rds`, one object per sample |
| 01 | `01_preprocessing/02_qc_overview.R` | `scRNA_final.rds` | per-sample QC violins, marker panel |
| 02 | `02_integration/01_cca_integration.R` | `objects/per_sample/<sid>_minf400.rds` | `results/integration/scRNA_integrated.rds`; `integrated_snn_res.1_anno`, the major cell types |
| 03 | `03_malignant_states/` | per-sample objects, major cell type annotation | 6 metaprograms, 5 malignant states — **not yet migrated, see below** |
| 04 | `04_lymphocyte_states/01_lymphocyte_states.R` | `scRNA_final.rds` | `results/lymphocyte_states/scRNA_lymphocyte_final.rds`, 11 states in `tnk_state` |
| 05 | `05_myeloid_states/01_myeloid_states.R` | `scRNA_final.rds` | `results/myeloid_states/scRNA_myeloid_final.rds`, `myeloid_state_coarse` (8) and `myeloid_state_fine` (12) |
| 06 | `06_infercnv/01_infercnv_profile.R` | `objects/per_sample/<sid>_minf400.rds`, `scRNA_final.rds`, `malignant_cluster_calls.csv`, `inputs/hg38.cellranger.genes.position.txt` | `results/infercnv/`, one inferCNV run with the cells grouped by sample |
| 07 | `07_k27m_reads/01_k27m_frequency.R` | `K27M_reads/<file_id>.K27M.reads.txt`, `scRNA_final.rds` | `results/k27m_reads/k27_alt_freq_by_sample.tsv`, mutant fraction per sample |
| 08 | `08_trajectory_TL6/01_slingshot_TL6.R` | `scRNA_final.rds` | `results/trajectory_TL6/`, pseudotime per sampling site for patient `TL_6`–`TL_10` |
| 99 | `99_export_figure_data.R` | the three released objects in `objects/released/` | every table in `data/` |
| 99b | `99b_export_major_panels.R` | `objects/released/scRNA_final.rds` | `data/panel_data/marker_gene_expression_nonzero.tsv.gz` |
| 99b | `99b_export_tnk_panels.R` | `objects/released/scRNA_lymphocyte_final.rds`, `data/cell_level/cell_annotation.tsv.gz`, `inputs/refsets/zhang2021_table_s3.xlsx` | `data/panel_data/tnk_*` |
| 99b | `99b_export_states_panels.R` | `results/malignant_states/`, `results/trajectory_TL6/`, `objects/released/scRNA_final.rds` | `data/metaprograms/`, malignant state, pseudotime and HOX tables in `data/panel_data/` |

Step 02 reads the deposited per-sample objects by default; to chain it to your own run
of step 01, set its `PREP_DIR` to `../../results/per_sample`.

Steps 04 to 08 read `scRNA_final.rds`, which is step 02's output with the state columns
of steps 03 to 05 merged in. The cycle is only apparent: step 02 writes
`scRNA_integrated.rds`, the state analyses are run on it, and the merged object is
released as `scRNA_final.rds`. Steps 04 to 08 are checked by running them against the
released object: what they write to `results/` should match the deposited subset objects
in `objects/released/` and the state columns of `scRNA_final.rds`.

## Order of the state analyses

A cell that a state analysis drops still keeps its state in `scRNA_final.rds`: the states
are assigned on the subset object and merged back, so `scRNA_final.rds` is the complete
record and the subset objects are the working objects. `substate_excluded_reason` says
why a cell has no state. The counts differ accordingly — 4,454 T/NK cells in
`scRNA_final.rds`, 4,201 in `scRNA_lymphocyte_final.rds`; 42,232 myeloid, 39,573 kept.

Cell counts to quote in text: 132,610 cells pass the per-sample QC; 132,343 remain after
the 267 cells of cluster 31 (`myeloid_mixed_T`, a mixture of myeloid and T cells) are
removed, and these are the cells of `scRNA_final.rds`; 81,304 of them are malignant.

The malignant composition analyses use the 74,090 malignant cells with
`included_in_prop == TRUE`. That flag downsamples one patient: the five `TL_6_*` samples
are five sites of the same tumour, and only 1,803 of their 9,017 malignant cells are
included; every malignant cell of the other 14 samples is. The subset was drawn without a
seed, so it cannot be redrawn; it is frozen in `inputs/frozen/included_in_prop.csv`.

## H3 K27 mutant reads

`07_k27m_reads/K27M_reads/` has 18 files for 19 samples. `PL_10` (`file_id` `MTSC01`) has
no file because no read covering the H3 K27 position was recovered in that sample. The
per-sample table that step 07 writes counts covered cells only, so `PL_10` is absent from
it; in `data/sample_level/k27_alt_freq_by_sample.tsv` it is listed with `n_reads_new = 0`.
`SP22`, which is `TL_5`, carries H3F3B K27I rather than H3F3A K27M and its file is named
`.K27I.` accordingly.

Step 07 calls a cell mutant when it has at least one alternative read
(`Alt_freq > 0`). The UMAP panel of the mutant reads and the `n_alt_ge_thr` column of
`data/sample_level/k27_alt_freq_by_sample.tsv` use `Alt_freq >= 0.4` instead.

## Run status

`03_malignant_states/` holds the consensus NMF analysis. Its scripts were migrated from the verified analysis and checked for syntax, but not re-run here (about 50 minutes). `results/` holds a local copy of the outputs of the verified run; it is not committed.
