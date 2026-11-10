# data -- the small tables the figures are drawn from

The scripts under `figures/` read only this directory, so every table here is
committed and redrawing a figure needs no download. Every table except the three
gene lists is written by `pipeline/`; none was edited by hand, so each can be
regenerated and diffed.

**Sample ids only.** These tables carry `SP01`..`SP36`; the hospital `PID` and
`Original_ID` columns go to `../data_full/`, which the GitHub version does not
ship. The two copies are written from the same data frame by
`common.R::write_tiered()`, so the numbers are identical.

## Signature gene lists (inputs, not outputs)

| File | Genes | Note |
|---|---|---|
| `tardon_genes.csv` | 19 | all 19 found in the TPM matrix |
| `buffa_genes.csv` | 51 | 48 found |
| `ragnum_genes.csv` | 32 | 30 found |

Taken from the reannotated compendium in `../raw_data/`; see that directory's
README for how to regenerate them.

## 01 -- scores and the PL vs TL comparison

| File | Rows | Contents |
|---|---|---|
| `hypoxia_scores_per_sample.csv` | 36 | one row per sample: the three scores, the number of signature genes used for each, the group, and the expression transform |
| `score_spearman_correlation.csv` | 3 | Spearman correlation between the three scores (0.55-0.86) |
| `group_comparison_hypoxia_scores.csv` | 3 | per signature: n, medians, Wilcoxon W and p, 95% CI, Cliff's delta |
| `group_comparison_hypoxia_scores.txt` | - | the same, human-readable, plus the correlation matrix |

## 02 -- MRI tumour enhancement (`mri/`)

| File | Rows | Contents |
|---|---|---|
| `scores_merged_mri_enhancement.csv` | 35 | scores joined to `Tumor Enhancement`, with the short label (Focal / Lept / Without) |
| `enhancement_group_counts.csv` | 3 | n per enhancement group |
| `kruskal_wallis_by_enhancement.csv` | 3 | Kruskal-Wallis per signature |
| `pairwise_wilcox_by_enhancement.csv` | 9 | pairwise Wilcoxon, BH-adjusted; the figure brackets read this |
| `enhancement_group_comparison.txt` | - | human-readable summary |

## 03 -- immune ratios vs scores (`immune_hypoxia/`)

| File | Rows | Contents |
|---|---|---|
| `merged_immune_hypoxia_per_sample.csv` | 29 | the samples with both mIHC and bulk data: immune ratios plus the three scores |
| `immune_hypoxia_spearman.csv` | 9 | rho, p and BH-adjusted p for every score x immune pair, with display labels |

## 04 -- immune ratios vs MRI enhancement (`immune_mri/`)

| File | Rows | Contents |
|---|---|---|
| `merged_immune_mri_per_sample.csv` | 28 | immune ratios joined to the enhancement label |
| `enhancement_group_counts.csv` | 3 | n per enhancement group |
| `kruskal_wallis_immune_by_enhancement.csv` | 3 | Kruskal-Wallis per immune column |
| `pairwise_wilcox_immune_by_enhancement.csv` | 9 | pairwise Wilcoxon, BH-adjusted |
| `enhancement_immune_group_comparison.txt` | - | human-readable summary |

`immune_total` is the sum of the `T_NK`, `Myeloid`, `Neutrophil` and `B_cell`
proportions. The sample counts differ between steps because mIHC covers 31
samples, of which 29 have bulk RNA and 28 have an MRI enhancement label.
