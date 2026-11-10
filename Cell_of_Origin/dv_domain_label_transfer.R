## =============================================================================
## dv_domain_label_transfer.R
##
## Projects the malignant cells of every tumour sample onto the dorsoventral (DV)
## progenitor domains of the developing human spinal cord, by reference-based label
## transfer (Seurat FindTransferAnchors + TransferData), and counts, for each sample,
## how many cells of each malignant state were assigned to each domain.
##
## Run from the Cell_of_Origin/ folder:   Rscript dv_domain_label_transfer.R
## It takes about 30 minutes (SCTransform and anchor finding for 19 samples). Its result
## is committed (data/dv_domain_cellstate_counts.tsv), so the heat map script
## dv_domain_heatmap.R does not need it.
##
## Input
##   data/HumanSpinalGW4toGW7_GSE171890.MetaData.csv
##       cell annotations of the human embryonic spinal cord atlas, public data set
##       GSE171890 (Rayon et al., Development 2021)
##   data/HumanSpinalGW4toGW7_GSE171890.integrated.RDS
##       integrated Seurat object of the same atlas (GSE171890); local copy, not tracked
##   ../Single_Cell_RNAseq/objects/integrated/integrated_minf400_only_PL_6_4_cca_anno_cleaned_mpstate.rds
##       integrated Seurat object of the 19 tumour samples with the malignant state
##       (mp_state) of every malignant cell; local copy, not tracked
##   ../Single_Cell_RNAseq/sample_metadata.tsv
##       sample names (sid used inside the object -> renamed_id used outside) and group
## Output
##   data/dv_domain_cellstate_counts.tsv
##       one row per sample x predicted domain x malignant state, n = number of cells
## =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

## ---- constants --------------------------------------------------------------
# the 13 neural progenitor domains along the dorsoventral axis (dorsal -> ventral)
dv_levels <- c('RP', paste0('dp', 1:6), paste0('p', c(0, 1, 2, 'MN', 3)), 'FP')

f_meta   <- 'data/HumanSpinalGW4toGW7_GSE171890.MetaData.csv'
f_ref    <- 'data/HumanSpinalGW4toGW7_GSE171890.integrated.RDS'
f_query  <- '../Single_Cell_RNAseq/objects/integrated/integrated_minf400_only_PL_6_4_cca_anno_cleaned_mpstate.rds'
f_sample <- '../Single_Cell_RNAseq/sample_metadata.tsv'
f_out    <- 'data/dv_domain_cellstate_counts.tsv'

## ---- 1. reference: annotated cells of the atlas, restricted to the 13 domains ----
mt <- read.delim(f_meta, sep = ',')
mt$cellid <- paste0(mt$orig.ident, '_', substr(mt$X, 1, 16), '-1')
# cells without a fine-grained annotation keep their coarse one
mt$Type_step2[is.na(mt$Type_step2)] <- mt$Type_step1[is.na(mt$Type_step2)]
mt$Type_step2 <- gsub(' ', '_', mt$Type_step2)

ref <- readRDS(f_ref)
ref$Type_step2 <- mt$Type_step2[match(rownames(ref@meta.data), mt$cellid)]
ref <- subset(ref, subset = Type_step2 %in% dv_levels)
DefaultAssay(ref) <- 'integrated'
message('reference cells in the 13 domains: ', ncol(ref))   # 5,460

## ---- 2. query: tumour samples, one label transfer per sample -----------------
query <- readRDS(f_query)
sids  <- unique(query$orig.ident)

df_list <- list()
for (sid in sids) {
  message(sid)
  seurat <- subset(query, subset = orig.ident == sid)
  if ('SCT' %in% names(seurat)) {
    DefaultAssay(seurat) <- 'SCT'
  } else {
    seurat <- SCTransform(seurat, assay = 'RNA', new.assay.name = 'SCT', verbose = TRUE)
  }
  anchors <- FindTransferAnchors(reference = ref, query = seurat, normalization.method = 'SCT',
                                 dims = 1:30, reference.reduction = 'pca')
  predictions <- TransferData(anchorset = anchors, refdata = ref$Type_step2, dims = 1:30)
  seurat <- AddMetaData(seurat, metadata = predictions)
  # predicted domain x malignant state; cells without a state (non-malignant) are dropped by table()
  df_list[[sid]] <- as.data.frame(table(seurat$predicted.id, seurat$mp_state))
}

## ---- 3. long table of counts, samples under their outward names --------------
smp <- read.delim(f_sample, colClasses = 'character')
stopifnot(all(names(df_list) %in% smp$sid))
counts_long <- do.call(rbind, lapply(names(df_list), function(sid) {
  d <- df_list[[sid]]
  data.frame(sample           = smp$renamed_id[match(sid, smp$sid)],
             group            = smp$group[match(sid, smp$sid)],
             predicted_domain = as.character(d$Var1),
             mp_state         = as.character(d$Var2),
             n                = d$Freq,
             stringsAsFactors = FALSE)
}))
write.table(counts_long, f_out, sep = '\t', quote = FALSE, row.names = FALSE)
message('wrote ', f_out, ' (', nrow(counts_long), ' rows)')
