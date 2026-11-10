# Single cell RNA-seq of H3 K27-altered diffuse spinal cord glioma

Code for the single cell RNA-seq analyses of 19 samples from 15 patients with H3 K27M
diffuse spinal cord glioma (H3-DSG): per-sample preprocessing, CCA integration and major
cell type annotation, malignant cell states by consensus NMF, T/NK and myeloid cell
states, copy number profiles by inferCNV, the fraction of cells carrying the H3 K27
mutant transcript, and the pseudotime trajectory across the five sampling sites of one
patient.

## How this directory is organised

Three layers, each with one kind of input:

| Layer | Reads | Writes | Run time |
|---|---|---|---|
| `pipeline/01` – `pipeline/08` | raw count matrices (GEO `GSExxxxxx`) and the deposited objects (Zenodo `10.5281/zenodo.xxxxxxx`) | Seurat objects and tables in `results/` | days, 200+ GB of disk |
| `pipeline/99_export_figure_data.R` | the released Seurat objects | the tables in `data/` | minutes |
| `figures/` | the tables in `data/` | the figure panels | seconds |

One rule keeps the layers apart, and it is the only thing a reader needs to know:

> **No script in `figures/` opens a Seurat object.** Every panel is drawn from a table in
> `data/`, and the tables in `data/` are exported from the released objects, all but two
> of them by `pipeline/99_export_figure_data.R` (`data/README.md` names the two).

So redrawing a panel needs this repository and nothing else — the tables it reads are
committed here, 12 MB in total. Re-deriving the objects the tables came from needs the
counts from GEO and a machine with ~128 GB of RAM. Re-deriving the tables from the
objects, if you want to check that the committed tables are what the objects say, needs
the Zenodo deposit and one script.

| Directory | Content | Committed |
|---|---|---|
| `pipeline/` | the analysis scripts, with the small per-step inputs next to them (`02_integration/major_cluster_annotation.csv`, `06_infercnv/malignant_cluster_calls.csv`, `07_k27m_reads/K27M_reads/`) | yes |
| `figures/` | one script per figure panel | yes |
| `data/` | the figure-ready tables | yes |
| `inputs/` | manual annotation decisions, reference gene sets, frozen inputs, gene positions (`inputs/README.md`) | yes, except `raw_counts/`, `meta/`, `sample_metadata.tsv`, the `.rds` files and the read files of samples outside the study |
| `inputs/raw_counts/` | local copy of the tier 0 count matrices | no |
| `objects/` | local copy of the tier 1 objects | no, only `objects/MANIFEST.tsv` (path, size, md5) |
| `results/` | what a pipeline run writes; holds a local copy of the verified run | no |
| `env/` | R `sessionInfo()` and the Python package list | yes |
| `sequencing_name_mapping.tsv` | single-cell sequencing names to Table S1 patient IDs | yes |
| `other_assay_id_mapping.tsv` | WES and RNA-seq numbers, methylation array IDs and GEO accessions to Table S1 patient IDs | yes |

## Patient identifiers

Patients are identified by the `Sample ID` of Table S1 (`SP01`..`SP86`); `patient_id` columns hold this ID. Samples are labelled by `renamed_id`.

Sequencing names such as `SP24` were given when the samples were sequenced and are **not** Table S1 IDs: sequencing name `SP24` is patient `SP38`, while patient `SP24` is someone else. Convert through `sequencing_name_mapping.tsv`, never by the name.

## Data tiers

| Tier | Content | Size | Where |
|---|---|---|---|
| 0 | raw count matrices, 19 samples (CellRanger / STARsolo filtered output) | 1.5 GB | GEO `GSExxxxxx`; local copy in `inputs/raw_counts/` |
| 1 | the three released Seurat objects; the 19 per-sample objects and the cNMF output they were built from | 7.4 GB full or 4.6 GB slimmed; 8.2 GB per-sample objects; 1.8 GB cNMF output | Zenodo `10.5281/zenodo.xxxxxxx`; local copy in `objects/` |
| 2 | per-cell annotation, per-sample statistics, one table per figure panel | 12 MB | `data/`, in this repository |

Tiers 0 and 1 are too large for version control. The small inputs in `inputs/` (34 MB)
are committed alongside tier 2; the whole repository is about 49 MB.

The three released objects are

| File | Cells | Contents |
|---|---|---|
| `scRNA_final.rds` | 132,343 | all samples; major cell type, malignant state, T/NK state, myeloid state, K27 mutant reads |
| `scRNA_lymphocyte_final.rds` | 4,201 | T/NK cells, 11 states in `tnk_state` |
| `scRNA_myeloid_final.rds` | 39,573 | myeloid cells, 8 coarse and 12 fine states |

Each object exists in a full and a slimmed form (`*_slim.rds`). The pipeline reads the
full form. In `scRNA_final.rds` it carries five reductions:

- `pca`, 30 principal components of the integrated assay;
- `umap_pc30`, on those 30 components with `n.neighbors = 150` and `min.dist = 0.1` —
  the embedding the major cell type annotation was defined on;
- `umap`, the Seurat default layout of the same components (`n.neighbors = 30`,
  `min.dist = 0.3`);
- `umap_uncorrected`, before integration;
- `umap_pc20_correlation`, on 20 components with `metric = "correlation"`, used for the
  K27 mutant read panel.

It also keeps the `scale.data` of the integrated assay. The slimmed form renames
`umap_pc30` to `umap`, which makes it the default of `DimPlot()`, drops the other two
UMAPs and every `scale.data` (call `ScaleData()` if you need it), and records this in
`misc(object)$reduction_notes` (all-cell object only). `inputs/README.md` lists every
form that is kept.

## Quick start

Redraw a panel. Nothing to download:

```r
setwd("figures")
source("umap_major_celltype.R")
```

Scripts and tables are named after what they show, not after a figure number — panel
numbers move between revisions. The panel index in `figures/README.md` is the one place a
figure number appears.

Re-export the tables in `data/` from the released objects. Download the Zenodo deposit
into `objects/` (the layout is in `inputs/README.md`), or leave it elsewhere and point
`SCRNA_RELEASE` at the directory holding the three released objects:

```r
setwd("pipeline")
# Sys.setenv(SCRNA_RELEASE = "/path/to/zenodo/deposit/released")   # only if not in objects/
source("99_export_figure_data.R")
```

Re-run the pipeline from raw counts. Put the GEO matrices in `inputs/raw_counts/`, one
directory per `file_id`, then run the numbered directories of `pipeline/` in order. Every
script uses paths relative to its own directory and is run from there; it reads from
`inputs/` and `objects/` and writes to `results/`, which is not committed, so a re-run
never overwrites a deposited object. Each script has its inputs and outputs in its header
comment, and `pipeline/README.md` has the full table.

## Cohort and sample naming

`sample_metadata.tsv` is the authoritative sample table, 19 rows. Two of its columns are
names for the same sample and the distinction matters when reading the code:

- `sid` — the key used in every file name, directory name and cNMF program name.
- `renamed_id` — the name used in the paper, in every figure and in every exported table.

They differ for six samples: `PL_6_4` is `PL_9`, and `TL_6_1` … `TL_6_5` are
`TL_6` … `TL_10`. The five `TL_6_*` samples are five sites along the cord of one patient,
which is why they share a `sid` stem but not a `renamed_id`. Patients are identified by
the `patient_id` column, never by the prefix of a name.

Gene symbols come from three reference generations: GRCh38-2020-A for 11 samples, a newer
annotation for 7, and hg19 with symbols only for `PL_9`. Comparisons across samples are
made on symbol equality, which costs about 2 genes in 100; symbols are not harmonised.

## Environment

R 4.5.1, Seurat 5.2.1, SeuratObject 5.1.0, inferCNV 1.24.0, slingshot 2.16.0.
Python 3.12.7 with cnmf 1.7.1 for the malignant state analysis. These are the versions of
the environment this release was prepared in. `env/sessionInfo_R.txt` is `sessionInfo()`
with every package the pipeline loads attached, and `env/conda_packages_python.txt` is
`conda list --export` of the Python environment.

The Python environment has pandas 3.0.5 with anndata 0.10.9, which do not work together
as installed: `pd.set_option("future.infer_string", False)` has to be set at module top
level of any script that runs cNMF, so that the worker processes started on Windows pick
it up when they re-import the module.
