# Differential expression between H3-DSG-PL and H3-DSG-TL tumours (bulk RNA-seq, discovery cohort)
# and the two gene lists submitted to Metascape for the functional enrichment networks:
#   genes higher in PL -> GO_network_Pons_like.cys
#   genes higher in TL -> GO_network_Thalamus_like.cys
#
# Input : ../Cell_of_Origin/data/Raw_SP_AB_RNA_count_data_with_gene_length.csv
#             raw read counts, genes x 40 tumours
#         ../Bulk_RNA_volcanoplot_for_enriched_pathways/label_for_twoGroups_SP36_1219.RDS
#             the 36 tumours of the discovery cohort and their subtype (YongG = PL, OldG = TL)
# Output: deseq2_PL_vs_TL_results.tsv   DESeq2 result of every gene tested (log2FoldChange > 0: higher in PL)
#         metascape_input_PL_up.txt     padj < 0.05 and log2FoldChange >  1, one gene symbol per line
#         metascape_input_TL_up.txt     padj < 0.05 and log2FoldChange < -1
#
# Metascape (https://metascape.org): paste one list, species H. sapiens, default settings
# (whole-genome background); the network is exported from there as a Cytoscape session.
#
# Run from this directory: Rscript deseq2_PL_vs_TL_gene_lists.R

suppressPackageStartupMessages(library(DESeq2))

counts <- read.csv("../Cell_of_Origin/data/Raw_SP_AB_RNA_count_data_with_gene_length.csv",
                   strip.white = TRUE, check.names = FALSE)
label  <- readRDS("../Bulk_RNA_volcanoplot_for_enriched_pathways/label_for_twoGroups_SP36_1219.RDS")

mat <- as.matrix(counts[, label$RNA])
rownames(mat) <- counts$Geneid
group <- factor(ifelse(label$group == "YongG", "PL", "TL"), levels = c("TL", "PL"))
stopifnot(sum(group == "PL") == 16, sum(group == "TL") == 20)

mat <- mat[rowSums(mat) >= 10, ]              # genes with fewer than 10 reads in total are not tested
dds <- DESeqDataSetFromMatrix(mat, data.frame(group, row.names = colnames(mat)), ~ group)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds, contrast = c("group", "PL", "TL")))
res <- cbind(gene = rownames(res), res)

sig   <- !is.na(res$padj) & res$padj < 0.05
pl_up <- res[sig & res$log2FoldChange >  1, ]
tl_up <- res[sig & res$log2FoldChange < -1, ]
pl_up <- pl_up[order(-pl_up$log2FoldChange), ]
tl_up <- tl_up[order(tl_up$log2FoldChange), ]

write.table(res, "deseq2_PL_vs_TL_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(pl_up$gene, "metascape_input_PL_up.txt")
writeLines(tl_up$gene, "metascape_input_TL_up.txt")
cat(nrow(res), "genes tested;", nrow(pl_up), "higher in PL;", nrow(tl_up), "higher in TL\n")
