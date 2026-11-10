# Overall survival, H3-DSG-PL vs H3-DSG-TL, restricted to similarly treated patients.
# Discovery and validation cohorts pooled, cases with survival data only.
#
# Three subsets, from the Radiotherapy and Chemotherapy columns of Table S1:
#   1) cases with radiotherapy                        RT == Yes
#   2) cases with chemotherapy                        Chemo == Yes
#   3) cases with both radiotherapy and chemotherapy   RT == Yes & Chemo == Yes
#
# Each subset is saved twice: once with the log-rank p drawn on the curves, and
# once without it (the _no_pval files), so the panel can be used either way. The
# group sizes are always in the subtitle, above the legend, which is why the
# ggsurvplot p value annotation is switched off in the second version.
#
# Output: Clinical_survival_comparison_radiotherapy.pdf           (+ _no_pval)
#         Clinical_survival_comparison_chemotherapy.pdf           (+ _no_pval)
#         Clinical_survival_comparison_radio_chemo.pdf            (+ _no_pval)
rm(list = ls())

# run from RStudio, or with Rscript from anywhere
source('required_packages.R')
source('support_function.R')

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

table_s1_path = 'Table S1.xlsx'
legend_labs = c("H3-DSG-PL", "H3-DSG-TL")

s1 <- readxl::read_xlsx(table_s1_path, skip = 1)   # row 1 is the table caption
stopifnot(names(s1)[1] == 'Sample ID', nrow(s1) == 86)

## ---- sentinel check, on the full validation table before any filtering ----
check_validation_sentinels(s1 %>% filter(`Discovery/Validation` == "Validation"))

## ---- pool both cohorts, keep the cases with survival data ----
surv_all <- s1 %>%
  transmute(
    Sample.ID = trimws(`Sample ID`),
    Cohort    = `Discovery/Validation`,
    OS        = suppressWarnings(as.numeric(`OS (Months)`)),
    Sensor    = suppressWarnings(as.numeric(Sensor)),
    Group     = factor(ifelse(grepl("^PL", trimws(`Group/Predict Group`)), "PL", "TL"),
                       levels = c("PL", "TL")),
    RT        = yesno(Radiotherapy),
    Chemo     = yesno(Chemotherapy)
  ) %>%
  filter(!is.na(OS), !is.na(Sensor))

cat(sprintf('pooled cohort with survival data, n = %d (discovery %d + validation %d)\n',
            nrow(surv_all), sum(surv_all$Cohort == "Discovery"),
            sum(surv_all$Cohort == "Validation")))

## ---- the three treatment subsets ----
subsets <- list(
  radiotherapy = list(
    data  = surv_all %>% filter(RT == "Yes"),
    label = 'cases with radiotherapy',
    title = 'Radiotherapy',
    tag   = 'radiotherapy'),
  chemotherapy = list(
    data  = surv_all %>% filter(Chemo == "Yes"),
    label = 'cases with chemotherapy',
    title = 'Chemotherapy',
    tag   = 'chemotherapy'),
  radio_chemo  = list(
    data  = surv_all %>% filter(RT == "Yes", Chemo == "Yes"),
    label = 'cases with both radiotherapy and chemotherapy',
    title = 'Radiotherapy + chemotherapy',
    tag   = 'radio_chemo')
)

for (nm in names(subsets)) {
  surdata <- subsets[[nm]]$data
  tag     <- subsets[[nm]]$tag

  n_pl <- sum(surdata$Group == "PL")
  n_tl <- sum(surdata$Group == "TL")

  cat(sprintf('\n=== %s: n = %d (discovery %d + validation %d) ===\n',
              subsets[[nm]]$label, nrow(surdata),
              sum(surdata$Cohort == "Discovery"), sum(surdata$Cohort == "Validation")))
  print(table(surdata$Group))

  logrank <- survdiff(Surv(OS, Sensor) ~ Group, data = surdata)
  pval_logrank <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
  cat(sprintf('log-rank p = %s\n', format(pval_logrank, digits = 4)))

  km_fit <- survfit(Surv(OS, Sensor) ~ Group, data = surdata)
  print(summary(km_fit)$table[, c("records", "events", "median")])

  plot_title <- subsets[[nm]]$title
  # the subtitle carries only the group sizes; whether the log-rank p appears on
  # the curves is controlled by show_pval, so the _no_pval version has no p at all
  plot_subtitle <- sprintf("%s n = %d vs %s n = %d",
                           legend_labs[1], n_pl, legend_labs[2], n_tl)

  for (show_pval in c(TRUE, FALSE)) {

    ggsurv <- ggsurvplot(
      km_fit,
      data = surdata,
      risk.table = TRUE,
      pval = show_pval,
      pval.size = 4,
      conf.int = FALSE,
      palette = c("#990000", "#2b8cbe"),   # PL dark red, TL blue
      xlim = c(0, 80),
      xlab = "Time in months",
      break.time.by = 12,
      ggtheme = theme_classic(),
      risk.table.y.text.col = TRUE,
      risk.table.height = 0.25,
      risk.table.y.text = FALSE,
      ncensor.plot = FALSE,
      ncensor.plot.height = 0.1,
      conf.int.style = "step",
      surv.median.line = "hv",
      # (n=xx) on its own line under each group name; the order follows the
      # factor levels c("PL","TL")
      legend.labs = legend_labs_with_n(legend_labs, c(n_pl, n_tl)),
      title = plot_title,
      subtitle = plot_subtitle
    )

    ggsurv <- customize_labels(
      ggsurv,
      font.caption   = c(14, "plain", "black"),
      font.x         = c(14, "italic", "black"),
      font.y         = c(14, "italic", "black"),
      font.xtickslab = c(12, "plain", "black"),
      font.ytickslab = c(12, "plain", "black")
    )

    # title and subtitle centred; in the ggplot layout both already sit above
    # the top legend
    ggsurv$plot <- ggsurv$plot +
      theme(
        plot.title    = ggtext::element_markdown(size = 13, face = "bold", colour = "black",
                                                 hjust = 0.5, margin = margin(b = 2)),
        plot.subtitle = ggtext::element_markdown(size = 10.5, colour = "black",
                                                 hjust = 0.5, margin = margin(b = 4))
      )

    fname <- paste0('Clinical_survival_comparison_', tag,
                    if (show_pval) '' else '_no_pval', '.pdf')
    ggsave(fname, ggsurv$plot, width = 5, height = 4.4)
    cat(sprintf('  saved %s\n', fname))
  }
}
