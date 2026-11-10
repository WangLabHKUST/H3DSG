# =============================================================================
# Fraction of malignant cells, PL vs TL
#
# Input : ../data/sample_level/cell_type_composition.tsv
# Output: malignant_fraction_PL_vs_TL.pdf / .png in this directory
#
# One point per patient: the mean over the patient's samples of the percentage of malignant
# cells among all cells of the sample. TL_6 .. TL_10 are five sites of one patient and
# contribute one point (15 patients: 8 PL, 7 TL). Two-sided Wilcoxon rank sum test between
# the patient values of PL and TL (wilcox.test, exact: no ties).
# =============================================================================

source("panel_helpers.R")

comp <- read.delim("../data/sample_level/cell_type_composition.tsv", colClasses = c(patient_id = "character"))

# percentage of the sample's cells in each cell type, 0 where the type is absent
samples <- unique(comp[, c("renamed_id", "group", "patient_id")])
sample_pct <- merge(samples, data.frame(category = names(CELLTYPE_LABELS)), by = NULL)
sample_pct$n <- comp$n[match(paste(sample_pct$renamed_id, sample_pct$category), paste(comp$renamed_id, comp$category))]
sample_pct$n[is.na(sample_pct$n)] <- 0L
sample_pct$percent <- 100 * sample_pct$n / ave(sample_pct$n, sample_pct$renamed_id, FUN = sum)

# the measurement of this panel, per sample and then per patient
sample_df <- sample_pct[sample_pct$category == "tumor", c("renamed_id", "group", "patient_id", "n", "percent")]
plot_df   <- aggregate(percent ~ patient_id + group, data = sample_df, FUN = mean)
plot_df   <- merge(plot_df, aggregate(renamed_id ~ patient_id, data = sample_df, FUN = length), by = "patient_id")
names(plot_df)[names(plot_df) == "renamed_id"] <- "n_samples"

test <- wilcox.test(percent ~ group, data = plot_df)

p <- group_box(plot_df, "Malignant", test$p.value, ylab = "Single cell proportion (%)",
               breaks = seq(0, 100, 20))
save_figure(p, "malignant_fraction_PL_vs_TL", width = 2.4, height = 4)
