# =============================================================================
# Fraction of cells carrying the H3 K27M (or K27I) mutant transcript
#
# Input : K27M_reads/<file_id>.K27M.reads.txt, per-cell read counts over the H3F3A
#         K27 position (SP22 carries H3F3B K27I instead), columns
#         cell_ID, depth_Ref, depth_Alt, depth_Total, MutationType, Variant
#         ../../objects/released/scRNA_final.rds for the cell to sample assignment
#         and the malignant call
# Output: ../../results/k27m_reads/k27_alt_freq_by_sample.tsv  per-sample coverage and
#                                                              mutant cell counts
#         ../../results/k27m_reads/K27_freq_malignant_vs_nonmalignant.pdf
#
# A cell counts as mutant when at least one alternative read is seen
# (Alt_freq = depth_Alt / depth_Total > 0). Cells with no coverage are left out of
# the denominator, so the fractions below are among covered cells only.
#
# Paths are relative to this directory: run the script from here.
# =============================================================================

library(Seurat)
library(ggplot2)

OBJ       <- '../../objects/released/scRNA_final.rds'
READS_DIR <- 'K27M_reads'
OUT_DIR   <- '../../results/k27m_reads'
SAMPLES   <- read.delim('../../sample_metadata.tsv', colClasses = 'character')
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

obj  <- readRDS(OBJ)
meta <- obj@meta.data
rm(obj); gc()
meta$is_malignant <- meta$integrated_snn_res.1_anno == 'tumor'

# 1. attach the K27 read counts, matching on the raw 10x barcode within each sample
meta$depth_Alt <- meta$depth_Total <- NA_real_
bc <- sub('-1$', '', meta$barcode)
for (i in seq_len(nrow(SAMPLES))) {
  sid <- SAMPLES$sid[i]; file_id <- SAMPLES$file_id[i]
  f <- list.files(READS_DIR, pattern = paste0('^', file_id, '\\.K27[MI]\\.reads\\.txt$'), full.names = TRUE)
  if (!length(f)) { message('no K27 reads for ', sid, ' (', file_id, ')'); next }
  reads <- read.csv(f[1], stringsAsFactors = FALSE)
  reads$cell_ID <- sub('-1$', '', reads$cell_ID)
  idx <- which(meta$orig.ident == sid)
  m <- match(bc[idx], reads$cell_ID)
  meta$depth_Alt[idx]   <- reads$depth_Alt[m]
  meta$depth_Total[idx] <- reads$depth_Total[m]
  message(sid, ': ', sum(!is.na(m)), ' of ', length(idx), ' cells covered')
}
meta$Alt_freq <- meta$depth_Alt / meta$depth_Total
meta$is_alt   <- ifelse(meta$Alt_freq > 0, 1L, 0L)

# 2. per-sample summary
covered <- meta[!is.na(meta$is_alt), ]
per_sample <- do.call(rbind, lapply(split(covered, covered$orig.ident), function(d) data.frame(
  sid          = d$orig.ident[1],
  n_covered    = nrow(d),
  n_alt        = sum(d$is_alt),
  median_depth = median(d$depth_Total),
  mal_prop     = mean(d$is_alt[d$is_malignant]),
  nmal_prop    = mean(d$is_alt[!d$is_malignant]))))
per_sample$renamed_id <- SAMPLES$renamed_id[match(per_sample$sid, SAMPLES$sid)]
per_sample$group      <- SAMPLES$group[match(per_sample$sid, SAMPLES$sid)]
print(per_sample, row.names = FALSE)
write.table(per_sample, file.path(OUT_DIR, 'k27_alt_freq_by_sample.tsv'),
            sep = '\t', quote = FALSE, row.names = FALSE)

# 3. malignant versus non-malignant cells, one point per sample
d <- rbind(data.frame(sid = per_sample$sid, variable = 'malignant',     value = per_sample$mal_prop),
           data.frame(sid = per_sample$sid, variable = 'non-malignant', value = per_sample$nmal_prop))
d <- d[!is.na(d$value), ]
p <- ggplot(d, aes(x = variable, y = 100 * value)) +
  geom_boxplot(aes(color = variable), outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  scale_color_manual(values = c(malignant = '#C1272D', `non-malignant` = '#377EB8')) +
  labs(x = NULL, y = '% cells carrying H3 K27M') +
  theme_classic()
ggsave(file.path(OUT_DIR, 'K27_freq_malignant_vs_nonmalignant.pdf'), p, width = 3, height = 3.5)
