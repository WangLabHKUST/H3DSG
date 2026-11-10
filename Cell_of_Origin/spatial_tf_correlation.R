## =============================================================================
## spatial_tf_correlation.R
##
## Similarity between the transcription-factor (TF) programme of malignant cells and
## every spot of a developing human spinal cord section (10x Visium). For each spot,
## the Spearman correlation between its TF expression and the pseudobulk TF expression
## of the malignant cells of the pons-like (PL) or the thalamus-like (TL) tumours.
##
## Run from the Cell_of_Origin/ folder:   Rscript spatial_tf_correlation.R
## Takes a few minutes (reads a 5 GB Seurat object); its result is committed
## (data/spatial_tf_correlation.tsv), so spatial_tf_correlation_plot.R does not need it.
##
## Input
##   ../Single_Cell_RNAseq/objects/integrated/integrated_minf400_only_PL_6_4_cca_anno_cleaned_mpstate.rds
##       integrated Seurat object of the 19 tumour samples; local copy, not tracked
##   data/TF_names_v_1.01.txt       1,639 human transcription factors (Lambert et al., Cell 2018)
##   data/spatial.integrated.RDS    integrated Seurat object of the Visium sections of the human
##                                  fetal spinal cord (Shi et al., Cell Res 2024; GEO GSE221692);
##                                  local copy, not tracked
##       section used: GSM6893166 (orig.ident GSM6893166_7_GW27, gestational week 27)
##   data/visium-v2_coordinates.txt Visium barcode -> array column / row (Space Ranger reference)
## Output
##   data/spatial_tf_correlation.tsv   one row per spot: barcode, coord_X, coord_Y, rho_PL, rho_TL
## =============================================================================

suppressPackageStartupMessages(library(Seurat))

f_query <- '../Single_Cell_RNAseq/objects/integrated/integrated_minf400_only_PL_6_4_cca_anno_cleaned_mpstate.rds'
f_tfs   <- 'data/TF_names_v_1.01.txt'
f_spat  <- 'data/spatial.integrated.RDS'
f_coord <- 'data/visium-v2_coordinates.txt'
f_out   <- 'data/spatial_tf_correlation.tsv'

slice_id         <- 'GSM6893166_7_GW27'
celltype_col     <- 'integrated_snn_res.1_anno'
malignant_states <- c('S', 'G2M', 'OC-like', 'OPC-like', 'MES-like', 'AC-like')

# mp_state is written with dots (AC.like), the state names above with hyphens (AC-like);
# compare them in one normalised form, otherwise only 'S' and 'G2M' would be kept
canon <- function(x) tolower(gsub('[-_. ]', '', as.character(x)))

## ---- 1. pseudobulk TF profile of the malignant cells of each group -------------
so  <- readRDS(f_query)
tfs <- read.delim(f_tfs, header = FALSE)

so1 <- subset(so, cells = colnames(so)[!is.na(so$mp_state) & canon(so$mp_state) %in% canon(malignant_states)])
n_cells <- table(so1$Group)
message('malignant cells: PL ', n_cells[['PL']], ', TL ', n_cells[['TL']])    # 52,339 / 21,577
message('TFs absent from the tumour data: ', length(setdiff(tfs$V1, rownames(so1[['RNA']]))))

# average of the normalised expression in non-log space (Seurat AverageExpression)
# The cells are chosen by name: subset(so1, subset = Group == group) would read 'group' as the
# metadata column of that name, not as the argument, and silently select every cell.
pseudobulk <- function(grp) {
  cells <- colnames(so1)[so1$Group == grp]
  stopifnot(length(cells) == n_cells[[grp]])
  ae <- AverageExpression(subset(so1, cells = cells), assays = 'RNA', features = tfs$V1, group.by = celltype_col)
  as.data.frame(ae$RNA)       # one column: all malignant cells carry the same annotation
}
tfs_ge_pons     <- pseudobulk('PL')
tfs_ge_thalamus <- pseudobulk('TL')
stopifnot(ncol(tfs_ge_pons) == 1, ncol(tfs_ge_thalamus) == 1, !isTRUE(all.equal(tfs_ge_pons[[1]], tfs_ge_thalamus[[1]])))
rm(so, so1); invisible(gc())

## ---- 2. spots of the Visium section ---------------------------------------------
spe  <- readRDS(f_spat)
spe1 <- subset(spe, subset = orig.ident == slice_id)
rm(spe); invisible(gc())
message('spots in ', slice_id, ': ', ncol(spe1))                               # 2,200

tmp <- as.data.frame(t(FetchData(spe1, vars = tfs$V1)))                        # TF x spot

## ---- 3. Spearman correlation of every spot with the pseudobulk profile ----------
# merge on TF name -> the TFs found in both data sets; no column positions are hard-coded
spot_cor <- function(ref) {
  m <- merge(ref, tmp, by = 0)
  n_ref <- ncol(ref)
  message('TFs in both data sets: ', nrow(m))                                  # 1,562
  as.data.frame(cor(m[, (n_ref + 2):ncol(m)], m[, 2:(n_ref + 1), drop = FALSE], method = 'spearman'))
}
cm_pons     <- spot_cor(tfs_ge_pons)
cm_thalamus <- spot_cor(tfs_ge_thalamus)
stopifnot(nrow(cm_pons) == ncol(spe1), identical(rownames(cm_pons), rownames(cm_thalamus)))

## ---- 4. array coordinates of the spots ------------------------------------------
cod <- read.delim(f_coord, header = FALSE)
cod$V1 <- paste0(cod$V1, '-1')
row.names(cod) <- cod$V1
cod <- cod[, -1]; names(cod) <- c('coord_X', 'coord_Y')

barcode <- substr(rownames(cm_pons), 1, 18)           # 16-nt barcode + '-1'
stopifnot(!anyDuplicated(barcode), all(barcode %in% rownames(cod)))
res <- data.frame(barcode = barcode,
                  coord_X = cod$coord_X[match(barcode, rownames(cod))],
                  coord_Y = cod$coord_Y[match(barcode, rownames(cod))],
                  rho_PL  = cm_pons[[1]],
                  rho_TL  = cm_thalamus[[1]])
write.table(res, f_out, sep = '\t', quote = FALSE, row.names = FALSE)
message('wrote ', f_out, ' (', nrow(res), ' spots)')

## ---- numbers quoted in the legend / methods --------------------------------------
lim <- quantile(c(res$rho_PL, res$rho_TL), c(0.02, 0.98))
message(sprintf('rho PL: %.3f-%.3f, median %.3f', min(res$rho_PL), max(res$rho_PL), median(res$rho_PL)))
message(sprintf('rho TL: %.3f-%.3f, median %.3f', min(res$rho_TL), max(res$rho_TL), median(res$rho_TL)))
message(sprintf('shared colour scale (2nd-98th percentile of the pooled values): %.3f-%.3f', lim[1], lim[2]))
message(sprintf('spots clipped by the scale: PL %.1f%%, TL %.1f%%',
                100 * mean(res$rho_PL < lim[1] | res$rho_PL > lim[2]), 100 * mean(res$rho_TL < lim[1] | res$rho_TL > lim[2])))
message(sprintf('Pearson r(PL, TL) across spots: %.3f', cor(res$rho_PL, res$rho_TL)))
