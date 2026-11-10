# Pseudobulk expression of 39 HOX genes in the single-cell RNA-seq data of 15 patients (19 samples).
# For each patient the raw counts of every HOX gene are summed over all cells of the patient's samples, then divided by
# the total counts of those cells (CPM). Used for Extended Data Fig. 2f (hox_pseudobulk_boxplot.R).
#
# Input:  ../Single_Cell_RNAseq/inputs/raw_counts/<file_id>/   10x raw count matrices of the 19 samples (barcodes, features, matrix)
#         ../Single_Cell_RNAseq/sample_metadata.tsv            sample table: file_id, renamed_id, group (PL / TL), patient_id (SP ID),
#                                                              tumor_location
# Output: data/hox_pseudobulk_by_patient.tsv                   one row per patient and gene
#
# Cells: all cells of the sample with at least 200 detected genes (the filter of Seurat's CreateSeuratObject(min.features = 200));
# no cell type or malignancy filter.
# Merging rule: patients with several samples (the five sampling sites of one patient) are one point. Their counts and their
# total counts are summed first and CPM is calculated afterwards (not the mean of per-sample CPM values).
# Missing genes: the reference of one sample (renamed_id PL_9) has 24 of the 39 HOX genes. The 15 others are NA for that sample,
# not 0 (they were not measured); a patient is NA for a gene only if no sample of the patient has it.
# Group: the group column of sample_metadata.tsv (project metadata); the Source Data of the figure group the patients by Table S1.
#
# Run time: about 15 minutes (reads 1.6 GB of count matrices).
# Run from the Cell_of_Origin/ directory:  Rscript hox_pseudobulk_build.R

suppressPackageStartupMessages({ library(Matrix); library(Seurat) })

raw_dir <- "../Single_Cell_RNAseq/inputs/raw_counts"
MIN_FEATURES <- 200L

# 39 HOX genes (HOXA8 and HOXA12 do not exist in the human genome); order: cluster A-D, then paralog number 1-13
HOX39 <- c(paste0("HOXA", c(1:7, 9:11, 13)), paste0("HOXB", c(1:9, 13)),
           paste0("HOXC", c(4:6, 8:13)),     paste0("HOXD", c(1, 3, 4, 8:13)))
stopifnot(length(HOX39) == 39)
HOX39 <- HOX39[order(substr(HOX39, 1, 4), as.integer(substr(HOX39, 5, 7)))]

meta <- read.delim("../Single_Cell_RNAseq/sample_metadata.tsv", colClasses = "character", check.names = FALSE)
stopifnot(nrow(meta) == 19, all(grepl("^SP[0-9]+$", meta$patient_id)), identical(meta$patient_id, paste0("S", meta$p0x_id)),
          all(meta$group %in% c("PL", "TL")))
have <- file.exists(file.path(raw_dir, meta$file_id, "matrix.mtx.gz"))
stopifnot("raw count matrix missing for some of the 19 samples" = all(have))

# ---- 1. per sample: HOX counts and total counts of the cells with >= 200 genes ----
hox_counts <- function(file_id) {
  mtx <- tryCatch(Read10X(file.path(raw_dir, file_id)),
                  error = function(e) Read10X(file.path(raw_dir, file_id), gene.column = 1))
  if (is.list(mtx)) mtx <- mtx[["Gene Expression"]]
  mtx <- mtx[, Matrix::colSums(mtx > 0) >= MIN_FEATURES, drop = FALSE]
  present <- intersect(HOX39, rownames(mtx))
  cnt <- setNames(rep(NA_real_, length(HOX39)), HOX39)                     # NA: gene is not in the reference
  cnt[present] <- Matrix::rowSums(mtx[present, , drop = FALSE])
  data.frame(gene = HOX39, count = as.numeric(cnt), total_counts = sum(mtx), n_cells = ncol(mtx), stringsAsFactors = FALSE)
}
by_sample <- do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
  message(sprintf("[%2d/19] %s", i, meta$renamed_id[i]))
  cbind(meta[i, c("renamed_id", "patient_id", "group", "tumor_location")], hox_counts(meta$file_id[i]), row.names = NULL)
}))

# ---- 2. per patient: sum counts and total counts of the samples, then CPM ----
merge_samples <- function(d) {
  s <- unique(d[, c("renamed_id", "total_counts", "n_cells")])
  data.frame(patient_id = d$patient_id[1], renamed_ids = paste(unique(d$renamed_id), collapse = ","), group = d$group[1],
             tumor_location = paste(unique(d$tumor_location), collapse = ";"), gene = d$gene[1],
             count = if (all(is.na(d$count))) NA_real_ else sum(d$count, na.rm = TRUE),
             total_counts = sum(s$total_counts), n_cells = sum(s$n_cells), n_samples = nrow(s), stringsAsFactors = FALSE)
}
by_patient <- do.call(rbind, lapply(split(by_sample, paste(by_sample$patient_id, by_sample$gene, sep = "\r")), merge_samples))
by_patient$HOX <- substr(by_patient$gene, 1, 4)
by_patient$int <- as.integer(substr(by_patient$gene, 5, 7))
by_patient$CPM <- 1e6 * by_patient$count / by_patient$total_counts
by_patient <- by_patient[order(by_patient$group, by_patient$patient_id, by_patient$HOX, by_patient$int),
                         c("patient_id", "renamed_ids", "group", "tumor_location", "gene", "HOX", "int", "count", "CPM",
                           "total_counts", "n_cells", "n_samples")]
rownames(by_patient) <- NULL
stopifnot(length(unique(by_patient$patient_id)) == 15, nrow(by_patient) == 15 * 39)
write.table(by_patient, "data/hox_pseudobulk_by_patient.tsv", sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
