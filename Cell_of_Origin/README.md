# Cell_of_Origin

Positional identity of H3-altered spinal cord gliomas: HOX expression and chromatin state along the rostro-caudal
axis, PAX3 and NKX6-1, mapping onto the dorsoventral progenitor domains of the developing spinal cord, and similarity
to a spatially resolved fetal spinal cord section. The chemokine ranking from the same bulk RNA-seq data is here too.

Run each script from the folder it is in. All paths are relative.

## Scripts

| Script | What it draws or writes | Reads | Runs from what is committed |
|---|---|---|---|
| `hox_bulk_logfc_brain_vs_spine.R` | log2 fold change against brain DMG, per gene, HOX genes in red | `data/` | yes |
| `hox_bulk_heatmap.R` | all 39 HOX genes in brain DMG and spinal tumours, ordered by tumour location | `data/` | yes |
| `hox_atac_k27me3_strips.R` | mean H3K27me3 and ATAC signal at HOX promoters per group | `data/hox_atac_k27me3_promoter_signal.tsv` | yes |
| `hox_pseudobulk_build.R` | writes `data/hox_pseudobulk_by_patient.tsv` (HOX pseudobulk per patient, all cells) | `../Single_Cell_RNAseq/inputs/raw_counts/` | no, needs the 10x matrices; about 15 min |
| `hox_pseudobulk_boxplot.R` | HOX pseudobulk in brain DMG and spinal tumours | `data/` | yes |
| `hox_dotplots/` | position marker gene sets and PAX3 / NKX6-1 by cell type, by group and by sample | see its `00_config.R` | plots yes; statistics need the integrated object |
| `dv_domain_label_transfer.R` | writes `data/dv_domain_cellstate_counts.tsv` (label transfer onto 13 progenitor domains) | reference atlas, integrated object | no, needs both objects; about 30 min |
| `dv_domain_heatmap.R` | mapping probability per progenitor domain, PL and TL | `data/dv_domain_cellstate_counts.tsv` | yes |
| `spatial_tf_correlation.R` | writes `data/spatial_tf_correlation.tsv` (Spearman correlation per Visium spot) | Visium object, integrated object | no, needs both objects; a few minutes |
| `spatial_tf_correlation_plot.R` | the correlation on the section, PL and TL, on one colour scale | `data/spatial_tf_correlation.tsv` | yes |
| `CXCL14_exp_rank.R` | exploratory bulk RNA-seq script; includes the chemokine ranking | `data/` | yes |
| `plot_ATACseq_K27Me3_tracksignals.R` | ATAC and H3K27me3 tracks at the HOX, PAX3 and NKX6-1 loci | bigWig files, not in this repository | no |

`hox_dotplots/` runs in the order 01 to 04, as listed in its `00_config.R`.

## Data

Committed in `data/`: expression tables, gene annotation, the tables written by the scripts above, and sample
annotation. Samples of this study are given as SP IDs (Table S1) in bulk tables and as `renamed_id` (`PL_*`, `TL_*`) in
single-cell tables. Samples of public data sets keep their public IDs.

Not committed (`.gitignore`), to be obtained from the source if the steps above are re-run:

| File | Source |
|---|---|
| `data/HumanSpinalGW4toGW7_GSE171890.integrated.RDS`, `data/HumanSpinalGW4toGW7_GSE171890.MetaData.csv` | developing human spinal cord atlas, Rayon et al., Development 2021, GEO GSE171890 |
| `data/spatial.integrated.RDS` | Visium sections of the fetal human spinal cord, Shi et al., Cell Res 2024, GEO GSE221692; section GSM6893166 is used |
| `../Single_Cell_RNAseq/objects/integrated/` | integrated objects of the 19 tumour samples, see `../Single_Cell_RNAseq/` |

Public brain DMG data: bulk RNA-seq from CGGA and St Jude, single-cell data from GSE210568.

## Notes

- **Groups.** `hox_bulk_logfc_brain_vs_spine.R` and `hox_bulk_heatmap.R` group tumours by the column `TSNE cluster` of
  `data/H3K27altered_project_metadata_updated.tsv` (A = PL, B = TL). This differs from Table S1 for SP02 and SP42. The
  other scripts follow the groups of Table S1.
- **Label transfer.** `dv_domain_label_transfer.R` is the transfer step of the original analysis, reduced to what
  `data/dv_domain_cellstate_counts.tsv` needs. A re-run of this script reproduces every count of that table.
- **TL_4 (C5-T1)** is counted with the cervical samples in the group dot plot; see `hox_dotplots/00_config.R`.
- **Correlation with sequencing depth.** The spot-level correlations in `data/spatial_tf_correlation.tsv` also track the
  number of genes detected per spot, so they show similarity, not an anatomical origin.
