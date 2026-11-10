# TP53 and NF1 mutation status, H3-DSG-PL vs H3-DSG-TL, in the validation cohort.
# Contingency table, Fisher exact test, proportions and in-bar counts are all
# computed from Table S1; nothing is hard coded.
#
# Set cohort_sel to 'Discovery' to run the same comparison on the discovery cohort.
#
# Output: TP53_mut_status_validation_PL_vs_TL.pdf
#         NF1_mut_status_validation_PL_vs_TL.pdf
rm(list = ls())

# run from RStudio, or with Rscript from anywhere
source('required_packages.R')
source('support_function.R')

library(ggplot2)
library(dplyr)

table_s1_path = 'Table S1.xlsx'
cohort_sel = 'Validation'
genes      = c('TP53', 'NF1')
panel_w    = 1.30      # panel size in inches, this is what sets the bar aspect
panel_h    = 1.85
grp_labs   = c(PL = "H3-DSG-PL", TL = "H3-DSG-TL")
fill_cols  = c(mutant = "#fec44f", wildtype = "#2c7fb8")

s1 <- readxl::read_xlsx(table_s1_path, skip = 1)   # row 1 is the table caption
stopifnot(names(s1)[1] == 'Sample ID', nrow(s1) == 86)

## ---- sentinel check, required before any validation PL/TL comparison ----
check_validation_sentinels(s1 %>% filter(`Discovery/Validation` == "Validation"))

dat <- s1 %>%
  filter(`Discovery/Validation` == cohort_sel) %>%
  transmute(
    Sample.ID = trimws(`Sample ID`),
    Group     = factor(ifelse(grepl("^PL", trimws(`Group/Predict Group`)), "PL", "TL"),
                       levels = c("PL", "TL")),
    TP53      = mut_status(TP53),
    NF1       = mut_status(NF1)
  )

cat(sprintf('%s cohort: n = %d\n', cohort_sel, nrow(dat)))
print(table(dat$Group))

## ---- one panel per gene: Fisher test and a stacked proportion plot ----
for (gene in genes) {
  d <- dat %>% filter(!is.na(.data[[gene]]))
  n_drop <- nrow(dat) - nrow(d)

  # rows = Mut/WT, columns = PL/TL
  tab <- table(factor(d[[gene]], levels = c("Mut", "WT")), d$Group)
  ft  <- fisher.test(tab)

  cat(sprintf('\n=== %s (%s cohort) ===\n', gene, cohort_sel))
  cat(sprintf('evaluable n = %d (%d cases with unknown status excluded)\n', nrow(d), n_drop))
  print(tab)
  cat(sprintf("Fisher exact p = %s; OR = %s (95%% CI %s-%s)\n",
              format(ft$p.value, digits = 4),
              format(unname(ft$estimate), digits = 3),
              format(ft$conf.int[1], digits = 3), format(ft$conf.int[2], digits = 3)))

  p <- plot_mut_prop_fisher(tab, gene = gene, group_labels = unname(grp_labs),
                            style = 'compact_inbar', x_text_angle = 45)

  # fix the panel, not the canvas, so the bars match across both genes
  save_fixed_panel(p,
                   paste0(gene, '_mut_status_', tolower(cohort_sel), '_PL_vs_TL.pdf'),
                   panel_width = panel_w, panel_height = panel_h)
}
