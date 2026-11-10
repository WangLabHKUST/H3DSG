# =============================================================================
# Copy number profiles of the malignant cells (inferCNV)
#
# Observations are the tumour cells of each sample, one group per sample. Normal
# references are pooled across samples: oligodendrocytes and up to 100 myeloid cells
# per sample, taken only from clusters that the tumour cell state analysis called
# normal. Oligodendrocyte and T/NK groups that are not used as reference are carried
# along as observation groups, as a sanity check that they come out flat.
#
# Input : ../../objects/per_sample/<sid>_minf400.rds   per-sample objects (counts)
#         ../../objects/released/scRNA_final.rds        cell annotation
#         malignant_cluster_calls.csv                   per-sample cluster calls from the
#                                                       tumour cell state analysis
#                                                       (malig / normal / exclude)
#         ../../inputs/hg38.cellranger.genes.position.txt  gene positions, from the
#                                                       CellRanger hg38 reference
# Output: ../../results/infercnv/, the inferCNV run directory, including infercnv.png
#         and the residual expression matrices
#
# Paths are relative to this directory: run the script from here.
#
# Runtime about 6 h on 10 threads. inferCNV 1.24, R 4.5.1.
# =============================================================================

library(Seurat)
library(Matrix)
library(infercnv)

OBJ        <- '../../objects/released/scRNA_final.rds'
PREP_DIR   <- '../../objects/per_sample'
GENE_ORDER <- '../../inputs/hg38.cellranger.genes.position.txt'
OUT_DIR    <- '../../results/infercnv'
SAMPLES    <- read.delim('../../sample_metadata.tsv', colClasses = 'character')
CALLS      <- read.csv('malignant_cluster_calls.csv', colClasses = c(sid = 'character', cluster = 'character'))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MAX_MYELOID   <- 100   # myeloid cells contributed to the reference pool per sample
MIN_OBS_GROUP <- 20    # smaller oligodendrocyte / T-NK groups are dropped
NUM_THREADS   <- 10
SEED          <- 123
set.seed(SEED)

# 1. roster: which cell plays which role ---------------------------------------
obj  <- readRDS(OBJ)
meta <- obj@meta.data
rm(obj); gc()

meta$call <- CALLS$final_call[match(paste(meta$orig.ident, meta$RNA_snn_res.0.5),
                                    paste(CALLS$sid, CALLS$cluster))]
stopifnot(!anyNA(meta$call))
anno    <- meta$integrated_snn_res.1_anno
cleaned <- as.logical(meta$in_cleaned_0806)

role <- rep(NA_character_, nrow(meta))
role[anno == 'tumor']   <- 'tumor'
role[anno == 'oligo']   <- 'oligo_obs'
role[anno == 'oligo'   & meta$call == 'exclude' & cleaned] <- 'oligo_ref'
role[anno == 'myeloid' & meta$call == 'normal'  & cleaned] <- 'myeloid_pool'
role[anno == 'T/NK'    & meta$call == 'normal'  & cleaned] <- 'tnk_obs'

roster <- data.frame(cell = rownames(meta), sid = as.character(meta$orig.ident),
                     role = role, stringsAsFactors = FALSE)
roster <- roster[!is.na(roster$role), ]

# keep at most MAX_MYELOID myeloid cells per sample, to keep the matrix manageable
myeloid_keep <- unlist(lapply(split(roster$cell[roster$role == 'myeloid_pool'],
                                    roster$sid[roster$role == 'myeloid_pool']),
                              function(x) if (length(x) > MAX_MYELOID) sample(x, MAX_MYELOID) else x),
                       use.names = FALSE)
roster$role[roster$role == 'myeloid_pool'] <- NA
roster$role[roster$cell %in% myeloid_keep] <- 'myeloid_ref'
roster <- roster[!is.na(roster$role), ]

drop_small <- function(df, r) {
  n <- table(df$sid[df$role == r])
  df[!(df$role == r & df$sid %in% names(n)[n < MIN_OBS_GROUP]), ]
}
roster <- drop_small(drop_small(roster, 'oligo_obs'), 'tnk_obs')
print(table(roster$role))

# 2. counts, from the per-sample objects ---------------------------------------
mats <- lapply(SAMPLES$sid, function(sid) {
  cells <- roster$cell[roster$sid == sid]
  if (!length(cells)) return(NULL)
  o <- readRDS(file.path(PREP_DIR, paste0(sid, '_minf400.rds')))
  cnt <- LayerData(o, assay = 'RNA', layer = 'counts')
  # cell names in the integrated object carry the sample prefix
  colnames(cnt) <- paste0(sid, '_', colnames(cnt))
  cnt[, intersect(cells, colnames(cnt)), drop = FALSE]
})
mats  <- Filter(Negate(is.null), mats)
genes <- Reduce(intersect, lapply(mats, rownames))
message('genes shared by all samples: ', length(genes))
counts <- do.call(cbind, lapply(mats, function(m) m[genes, , drop = FALSE]))
rm(mats); gc()

roster <- roster[roster$cell %in% colnames(counts), ]
counts <- counts[, roster$cell, drop = FALSE]
message('counts matrix: ', nrow(counts), ' x ', ncol(counts))

gene_pos <- read.delim(GENE_ORDER, header = FALSE, stringsAsFactors = FALSE)
names(gene_pos)[1] <- 'gene'
gene_pos <- gene_pos[gene_pos$gene %in% genes, ]
genes_txt <- file.path(OUT_DIR, 'gene_order_used.txt')
write.table(gene_pos, genes_txt, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
counts <- counts[rownames(counts) %in% gene_pos$gene, , drop = FALSE]

# 3. groups --------------------------------------------------------------------
group <- rep(NA_character_, nrow(roster))
group[roster$role == 'tumor']       <- roster$sid[roster$role == 'tumor']
group[roster$role == 'oligo_obs']   <- paste0('oligoQ_', roster$sid[roster$role == 'oligo_obs'])
group[roster$role == 'tnk_obs']     <- paste0('ctrl_',   roster$sid[roster$role == 'tnk_obs'])
group[roster$role == 'oligo_ref']   <- 'ref_oligo'
group[roster$role == 'myeloid_ref'] <- 'ref_myeloid'
ref_groups <- c('ref_oligo', 'ref_myeloid')
print(table(group))

anno_file <- file.path(OUT_DIR, 'annotations.txt')
write.table(data.frame(roster$cell, group), anno_file,
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

# 4. run -----------------------------------------------------------------------
options(scipen = 100)
cnv <- CreateInfercnvObject(raw_counts_matrix = as.matrix(counts),
                            annotations_file  = anno_file,
                            delim             = '\t',
                            gene_order_file   = genes_txt,
                            ref_group_names   = ref_groups)
rm(counts); gc()

cnv <- infercnv::run(cnv, cutoff = 0.1, out_dir = OUT_DIR,
                     cluster_by_groups = TRUE, denoise = TRUE,
                     sd_amplifier = 3, noise_logistic = TRUE, HMM = FALSE,
                     num_threads = NUM_THREADS, plot_steps = FALSE,
                     inspect_subclusters = FALSE)
message('done -> ', OUT_DIR)
