# Bulk RNA hypoxia signature scores

Hypoxia scoring of 36 spinal cord glioma cases from bulk RNA-seq, compared between
the two anatomical groups (PL vs TL), against MRI tumour enhancement, and against
mIHC immune-rich cell ratios.

Three published signatures are scored, each with the summary score that the
benchmark of Di Giovannantonio et al. found to work best for it:

| Signature | Genes used | Summary score | Source |
|---|---|---|---|
| Tardon | 19 of 19 | interquartile mean | Calvo Tardon et al., *Biology* 9, 264 (2020), doi:10.3390/biology9090264 |
| Buffa | 48 of 51 | mean | Buffa et al., *Br J Cancer* 102, 428-435 (2010) |
| Ragnum | 30 of 32 | interquartile mean | Ragnum et al., *Br J Cancer* 112, 382-390 (2015) |

Gene lists were taken from the reannotated compendium in Di Giovannantonio et al.,
*Cell Genomics* (2025), doi:10.1016/j.xgen.2025.100764, Table S1a (the Tardon column
holds the 19 hypoxia-upregulated genes of the original 36-gene signature). The
genes that are missing for Buffa and Ragnum are simply absent from the TPM matrix.

## How this directory is organised

Three layers; each layer reads only one kind of input.

| Layer | Reads | Writes | Runtime |
|---|---|---|---|
| `pipeline/` | `raw_data/` (+ the scores table from step 01) | the small tables in `data/`, the same tables with patient ids in `data_full/` | under a minute |
| `figures/` | `data/` only | the PDFs and PNGs in `figures/` | seconds |

**The one rule:**

> **No script under `figures/` ever opens `raw_data/`.** Every figure is drawn from a
> table in `data/`, and every table in `data/` is written by `pipeline/`.

So the GitHub version, which ships no raw data, still redraws every figure.

```bash
Rscript run_all.R            # pipeline then figures (needs raw_data/)
Rscript run_all.R figures    # figures only (what the GitHub version can do)
```

Each script also runs on its own, in any order after step 01, and can be sourced
from RStudio as long as the working directory is somewhere inside the module.

## Steps

| Script | Reads | Writes |
|---|---|---|
| `pipeline/01_hypoxia_scores.R` | TPM matrix, ID map, the three gene lists | `data/hypoxia_scores_per_sample.csv`, `data/score_spearman_correlation.csv`, `data/group_comparison_hypoxia_scores.{csv,txt}` |
| `pipeline/02_mri_enhancement_stats.R` | scores, clinical table | `data/mri/` |
| `pipeline/03_immune_hypoxia_stats.R` | scores, cell-type proportions, ID map | `data/immune_hypoxia/` |
| `pipeline/04_immune_mri_stats.R` | cell-type proportions, ID map, clinical table | `data/immune_mri/` |
| `figures/01_scores_by_group.R` | `data/` | `fig_<signature>_by_group.{pdf,png}` |
| `figures/02_scores_by_mri.R` | `data/mri/` | `fig_<signature>_by_tumor_enhancement.{pdf,png}` |
| `figures/03_immune_hypoxia_scatter.R` | `data/immune_hypoxia/` | `scatter_<score>_vs_<immune>.{pdf,png}`, `scatter_facet_all.{pdf,png}` |
| `figures/04_immune_by_mri.R` | `data/immune_mri/` | `fig_<immune>_by_tumor_enhancement.{pdf,png}` |

`common.R` holds everything shared: paths, cohort constants, the statistics
(interquartile mean, Cliff's delta, Kruskal-Wallis, pairwise Wilcoxon with BH
adjustment, Spearman) and the plot theme. Paths are resolved from its own
location, so the module can be moved or renamed without editing any script.

## Statistics

- **PL vs TL**: two-sided Wilcoxon rank-sum test, with Cliff's delta as the effect
  size. Reported in `data/group_comparison_hypoxia_scores.csv`.
- **MRI tumour enhancement** (Focal / Lept / Without): Kruskal-Wallis across the
  three groups, then pairwise Wilcoxon with Benjamini-Hochberg adjustment; the
  brackets in the figures carry those adjusted p-values.
- **Immune ratios vs scores**: Spearman correlation, BH-adjusted across all
  score x immune-column pairs tested.
- Box plots never draw outliers (`outlier.shape = NA`), because the jitter layer
  already shows every sample; with the default the same sample would appear twice.

## Two versions

| Version | Contents | Size | Purpose |
|---|---|---|---|
| **Archive** | this directory, including `raw_data/` and `data_full/` | ~22 MB (20 MB of it `raw_data/`) | hand the whole folder over; every step from raw input to figure is present |
| **GitHub** | the same code, the tables in `data/`, the figures | ~2 MB | public repository |

There is one copy of the code. The GitHub version is this directory minus what `.gitignore` excludes.

## Patient identifiers

`data/` carries the published sample ids (`SP01`..`SP36`) only. The hospital
`PID` and `Original_ID` columns are written to `data_full/`, which is excluded
from the GitHub version, and the ID map that links them lives in `raw_data/`.
`common.R::write_tiered()` is what enforces this: it writes both copies from the
same data frame, dropping the identifier columns from the `data/` one.

## Requirements

R 4.5.1 with `ggplot2` (4.0.1), `ggsignif` (0.6.4) and `readxl` (1.4.5).

`common.R` carries one workaround at the top: `+<-` entered base R in 4.5.2, and
ggplot2 4.0 is compiled against it, so on R 4.5.1 every `p + layer` would fail
with "could not find function '+<-'". Delete those two lines once R is upgraded.
