# NF1 mutation status, this cohort against published brain DMG cohorts. Two panels:
#   A: SP DMG (all 36 discovery cases) vs brain DMG (published)
#   B: H3-DSG-TL (the TL subgroup of the discovery cohort, n = 20) vs brain
#      thalamus (published)
#
# Counts for this cohort are computed from Table S1 and cross-checked against the
# expected values; the script stops if they disagree. NF1 is complete for all 36
# discovery cases.
#   SP DMG      15 mut / 21 WT
#   H3-DSG-TL   14 mut /  6 WT   (discovery cohort only; the validation TL
#                                 subgroup is 15 mut / 9 WT)
#
# SP36 is NF1 mutant: Table S2 lists a frameshift at chr17:31337558 A>AGAAG from
# the hg38 re-call. The earlier hg19 call had it as wild type, which is why older
# versions of this comparison read 14/22 and 13/7.
#
# The two published groups have no per-sample data in this repository and are
# entered as the published counts:
#   brain DMG        3 mut / 35 WT  (n = 38)
#   brain thalamus   3 mut / 11 WT  (n = 14), the same n as the thalamus total in
#                    the key driver frequency table
#
# Output: NF1_mut_status_SPDMG_vs_BrainDMG.pdf
#         NF1_mut_status_TL_vs_Thalamus.pdf
rm(list = ls())

# run from RStudio, or with Rscript from anywhere
source('required_packages.R')
source('support_function.R')

library(ggplot2)
library(dplyr)

table_s1_path = 'Table S1.xlsx'
gene     = 'NF1'
panel_w  = 0.90        # panel size in inches, this is what sets the bar aspect
panel_h  = 2.05

s1 <- readxl::read_xlsx(table_s1_path, skip = 1)   # row 1 is the table caption
stopifnot(names(s1)[1] == 'Sample ID', nrow(s1) == 86)

dis <- s1 %>%
  filter(`Discovery/Validation` == "Discovery") %>%
  transmute(Group = ifelse(grepl("^PL", trimws(`Group/Predict Group`)), "PL", "TL"),
            NF1   = mut_status(NF1)) %>%
  filter(!is.na(NF1))

## ---- compute from the table, then cross-check ----
sp_dmg <- c(Mut = sum(dis$NF1 == "Mut"), WT = sum(dis$NF1 == "WT"))
tl     <- c(Mut = sum(dis$NF1 == "Mut" & dis$Group == "TL"),
            WT  = sum(dis$NF1 == "WT"  & dis$Group == "TL"))

cat("cross-check (computed from the discovery cohort vs expected):\n")
cat(sprintf("  SP DMG    : got %d mut / %d WT  expected 15 mut / 21 WT  [%s]\n",
            sp_dmg[["Mut"]], sp_dmg[["WT"]],
            if (all(sp_dmg == c(15, 21))) "OK" else "MISMATCH"))
cat(sprintf("  H3-DSG-TL : got %d mut / %d WT  expected 14 mut / 6 WT   [%s]\n",
            tl[["Mut"]], tl[["WT"]],
            if (all(tl == c(14, 6))) "OK" else "MISMATCH"))
stopifnot(all(sp_dmg == c(15, 21)), all(tl == c(14, 6)))

## ---- published comparison groups (no per-sample data available) ----
brain_dmg      <- c(Mut = 3,  WT = 35)   # n = 38
brain_thalamus <- c(Mut = 3,  WT = 11)   # n = 14

comparisons <- list(
  SPDMG_vs_BrainDMG = list(
    tab    = cbind(`SP DMG` = sp_dmg, `Brain DMG` = brain_dmg),
    labels = c("SP DMG", "Brain DMG"),
    tag    = 'SPDMG_vs_BrainDMG'),
  TL_vs_Thalamus = list(
    tab    = cbind(`H3-DSG-TL` = tl, `Thalamus` = brain_thalamus),
    labels = c("H3-DSG-TL", "Thalamus"),
    tag    = 'TL_vs_Thalamus')
)

for (nm in names(comparisons)) {
  tab <- comparisons[[nm]]$tab
  rownames(tab) <- c("Mut", "WT")

  cat(sprintf('\n=== %s ===\n', nm))
  print(tab)
  ft <- fisher.test(tab)
  cat(sprintf("Fisher exact p = %s; OR = %s (95%% CI %s-%s)\n",
              format(ft$p.value, digits = 4), format(unname(ft$estimate), digits = 3),
              format(ft$conf.int[1], digits = 3), format(ft$conf.int[2], digits = 3)))

  p <- plot_mut_prop_fisher(tab, gene = gene, group_labels = comparisons[[nm]]$labels)

  # fix the panel, not the canvas, so the bars match across both panels
  save_fixed_panel(p,
                   paste0(gene, '_mut_status_', comparisons[[nm]]$tag, '.pdf'),
                   panel_width = panel_w, panel_height = panel_h)
}
