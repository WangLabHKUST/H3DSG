# Volcano plot of the E2F target genes (bulk RNA-seq)

Expression of the HALLMARK_E2F_TARGETS genes in H3-DSG-PL against H3-DSG-TL tumours of the discovery cohort
(36 tumours: 16 H3-DSG-PL, 20 H3-DSG-TL).

```bash
cd Bulk_RNA_E2F_targets_volcano_plot
Rscript E2F_targets_volcano_plot.R
```

## What is plotted

| | |
|---|---|
| x axis | log2(mean TPM in H3-DSG-PL − mean TPM in H3-DSG-TL): the log2 of the **difference** between the two group means. It is not a fold change (ratio). |
| y axis | −log10 P of a two-sided Student t test (equal variances) on TPM |
| red points | x > 1 and P < 0.05 |
| dashed lines | x = −1, x = 1 and P = 0.05 |

TPM is computed from the gene-level read counts (hg19): counts divided by gene length, scaled to 1e6 per tumour
over all genes of the count table.

Of the 200 genes of the gene set, 189 are plotted (162 red, 27 grey):

- 5 genes are not found in the count table under their current symbol (CNOT9, H2AX, H2AZ1, JPT1, MRE11; the table
  lists them under an older symbol) and are not drawn.
- 6 genes have a mean TPM that is not higher in H3-DSG-PL (AK2, ANP32E, BRMS1L, CDKN1A, PRDX4, PRPS1). The log2 of
  their difference is undefined, so they are not drawn.

## Files

| File | Content |
|---|---|
| `E2F_targets_volcano_plot.R` | computes the table and draws the plot |
| `HALLMARK_E2F_TARGETS_genes.tsv` | the 200 genes of the gene set |
| `sample_groups_discovery36.tsv` | subtype of the 36 tumours |
| `E2F_targets_volcano_table.tsv` | output: one row per gene (group means, difference, log2 difference, P, colour) |
| `E2F_targets_volcano.pdf` | output: the plot |

The read counts are read from `../Cell_of_Origin/data/Raw_SP_AB_RNA_count_data_with_gene_length.csv`.

## Relation to the submitted panel

The original plotting script was not kept; this script was rebuilt from the submitted panel. It gives the same
189 points with the same colours. For the labelled genes whose point can be read from the panel without
ambiguity, the y positions agree within 0.03, and the x positions on the panel are a constant 0.03 higher than
computed here (a factor of 1.02 on the difference of the means). The cause of this constant was not identified.

## Gene set

HALLMARK_E2F_TARGETS, MSigDB release 2025.1.Hs (Liberzon et al., Cell Systems 2015), exported from the R package
msigdbr 25.1.1. `gene_symbol` is the current symbol, `db_gene_symbol` the symbol stored in MSigDB.
