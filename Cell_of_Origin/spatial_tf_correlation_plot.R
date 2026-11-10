## =============================================================================
## spatial_tf_correlation_plot.R
##
## Spatial maps of the similarity between the transcription-factor programme of malignant
## cells (pons-like, PL, and thalamus-like, TL, tumours) and the spots of a gestational
## week 27 human spinal cord section (10x Visium, GEO GSM6893166 of GSE221692).
## Both maps share one colour scale, clipped at the 2nd and 98th percentile of the pooled
## values (the colour only saturates, no value is changed).
##
## Run from the Cell_of_Origin/ folder:   Rscript spatial_tf_correlation_plot.R
##
## Input   data/spatial_tf_correlation.tsv   written by spatial_tf_correlation.R
##         (barcode, coord_X, coord_Y, rho_PL, rho_TL)
## Output  spatial_tf_correlation_PL_TL.pdf
## =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

f_in  <- 'data/spatial_tf_correlation.tsv'
f_pdf <- 'spatial_tf_correlation_PL_TL.pdf'

d <- read.delim(f_in, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(d$barcode), !anyNA(d[, c('coord_X', 'coord_Y', 'rho_PL', 'rho_TL')]))

# one colour scale for both maps: 2nd and 98th percentile of the pooled values
lim <- quantile(c(d$rho_PL, d$rho_TL), c(0.02, 0.98))
message(sprintf('spots: %d; colour scale %.3f-%.3f; clipped: PL %.1f%%, TL %.1f%%', nrow(d), lim[1], lim[2],
                100 * mean(d$rho_PL < lim[1] | d$rho_PL > lim[2]), 100 * mean(d$rho_TL < lim[1] | d$rho_TL > lim[2])))

jet.colors <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")

spot_map <- function(rho, title) {
  ggplot(d, aes(x = coord_X, y = coord_Y)) +
    geom_point(aes(fill = .data[[rho]]), size = 1.4, stroke = 0.1, pch = 21, color = '#999999') +
    scale_fill_gradientn(colours = jet.colors, limits = lim, oob = scales::squish) +
    coord_fixed(ratio = 1.5) +
    theme_void() +
    labs(fill = 'Spearman\ncorrelation\ncoefficient', title = title) +
    theme(plot.title = element_text(hjust = 0.5, size = 11), plot.margin = margin(5, 5, 5, 5))
}

p <- spot_map('rho_PL', 'Pons-like') + spot_map('rho_TL', 'Thalamus-like') +
  plot_layout(nrow = 1, guides = 'collect')

ggsave(f_pdf, p, width = 547 / 72, height = 244 / 72)
message('wrote ', f_pdf)
