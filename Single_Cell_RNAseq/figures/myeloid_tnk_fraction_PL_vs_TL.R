# =============================================================================
# Fraction of myeloid and T/NK cells, PL vs TL
#
# Input : ../data/sample_level/cell_type_composition.tsv
# Output: myeloid_tnk_fraction_PL_vs_TL.pdf / .png in this directory
#
# One point per patient: the mean over the patient's samples of the percentage of the cell
# type among all cells of the sample. TL_6 .. TL_10 are five sites of one patient and
# contribute one point (15 patients: 8 PL, 7 TL). Two-sided Wilcoxon rank sum test between
# the patient values of PL and TL, one test per cell type (wilcox.test, exact: no ties).
# =============================================================================

source("panel_helpers.R")
suppressPackageStartupMessages(library(patchwork))

comp <- read.delim("../data/sample_level/cell_type_composition.tsv", colClasses = c(patient_id = "character"))

samples <- unique(comp[, c("renamed_id", "group", "patient_id")])
sample_pct <- merge(samples, data.frame(category = names(CELLTYPE_LABELS)), by = NULL)
sample_pct$n <- comp$n[match(paste(sample_pct$renamed_id, sample_pct$category), paste(comp$renamed_id, comp$category))]
sample_pct$n[is.na(sample_pct$n)] <- 0L
sample_pct$percent <- 100 * sample_pct$n / ave(sample_pct$n, sample_pct$renamed_id, FUN = sum)

TYPES <- c(myeloid = "Myeloid", "T/NK" = "T/NK")      # raw cell_type -> title

plot_df <- do.call(rbind, lapply(names(TYPES), function(ct) {
  s <- sample_pct[sample_pct$category == ct, ]
  d <- aggregate(percent ~ patient_id + group, data = s, FUN = mean)
  d$celltype <- TYPES[[ct]]
  d
}))
tests <- lapply(TYPES, function(lab) wilcox.test(percent ~ group, data = plot_df[plot_df$celltype == lab, ]))

panels <- lapply(seq_along(TYPES), function(i) {
  lab <- TYPES[[i]]
  group_box(plot_df[plot_df$celltype == lab, ], lab, tests[[i]]$p.value,
            ylab = if (i == 1) "Single cell proportion (%)" else NULL)
})
p <- panels[[1]] + panels[[2]]

save_figure(p, "myeloid_tnk_fraction_PL_vs_TL", width = 4.6, height = 4)
