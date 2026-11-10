library(tidyverse)

# -----------------------------------------------------------------------
# Entropy contour density maps
#
# Input:  ../data/grid_data_entropy_tumor_vs_non_malig_entropy_{Group1,Group2,diff}.csv
#         written by pipeline/02_entropy/02_kde_grid.py
#         (upstream of which is pipeline/02_entropy/01_cell_level_stats.py)
# Output: density_r.{PL,TL}.{auto,max}.png and density_r_diff.png
#         (all written into this directory)
#
# x axis entropy_tumor     = Shannon entropy of the malignant states within the 27.5um
#                            neighbourhood
# y axis non_malig_entropy = Shannon entropy of the non-malignant cell types in the same
#                            neighbourhood
#
# Group naming: Group1/Group2 in the data layer, PL/TL in the figures (titles and file
# names)
#
# This script reads only the three 400-row grid CSVs in ../data/ and never touches the
# raw data, so the GitHub version (which has no raw_data/) can redraw the figures
# directly.
# -----------------------------------------------------------------------

# The directory of this script (resolved automatically under Rscript; when running
# interactively, setwd to this directory first)
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep('^--file=', args, value = TRUE)
  if (length(m) > 0) return(dirname(normalizePath(sub('^--file=', '', m[1]))))
  normalizePath(getwd())
}

script_dir <- get_script_dir()
data_path <- normalizePath(file.path(script_dir, '..', 'data'), mustWork = TRUE)
fig_path  <- script_dir

var_x <- 'entropy_tumor'
var_y <- 'non_malig_entropy'
var_z <- 'density'
pair_key <- paste0(var_x, '_vs_', var_y)

# The data layer always uses Group1/Group2 and the figures always use PL/TL.
# Kept consistent with group_labels in constants.py (R cannot read the .py, so an
# explicit copy lives here).
group_labels <- c('Group1' = 'PL', 'Group2' = 'TL')

# Canonical display names for the axes
axis_labels <- c(
  'entropy_tumor'     = 'Malignant-state entropy',
  'non_malig_entropy' = 'Non-malignant cell-type entropy'
)

# Upper limit of the y axis: very few grid points have a non-malignant entropy above 1.2,
# so they are clipped to keep the colour scale from being compressed
y_upper <- 1.2
# The density cap of the 'max' version, which lets both groups share one colour scale and
# be compared directly
max_density <- 20

read_grid <- function(tag) {
  df <- read_csv(file.path(data_path, sprintf('grid_data_%s_%s.csv', pair_key, tag)),
                 show_col_types = FALSE)
  df[df[[var_y]] < y_upper, ]
}

base_plot <- function(df, title) {
  ggplot(df, aes(x = .data[[var_x]], y = .data[[var_y]], z = .data[[var_z]])) +
    geom_raster(aes(fill = .data[[var_z]]), interpolate = TRUE) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, y_upper), expand = FALSE) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = 'bold'),
      plot.title = element_text(hjust = 0.5, face = 'bold'),
      legend.position = 'right',
      legend.key.height = unit(1.2, 'cm')
    ) +
    labs(title = title, x = axis_labels[[var_x]], y = axis_labels[[var_y]])
}

# ---------------------- 1. density map for each group ----------------------
# Number of contour levels: the two groups have differently shaped distributions, so each
# gets its own degree of subdivision
contour_num_mapping <- c('Group1' = 5, 'Group2' = 10)

for (cutoff in c('auto', 'max')) {
  for (group in c('Group1', 'Group2')) {
    group_disp <- group_labels[[group]]
    df <- read_grid(group)
    if (cutoff == 'max') {
      df[df[[var_z]] > max_density, var_z] <- max_density
    }

    p <- base_plot(df, group_disp) +
      geom_contour(aes(z = .data[[var_z]]),
                   color = 'black', linewidth = 0.05, bins = 15)

    if (cutoff == 'auto') {
      this_max <- max(df[[var_z]])
      p <- p + scale_fill_distiller(
        palette   = 'YlOrRd',
        direction = 1,
        name      = 'Density',
        breaks    = seq(0, this_max, by = round(this_max / 5))
      )
    } else {
      p <- p + scale_fill_distiller(
        palette   = 'YlOrRd',
        direction = 1,
        limits    = c(0, max_density),
        name      = 'Density'
      )
    }

    # The file names use PL/TL too, so nothing has to be translated on delivery
    ggsave(file.path(fig_path, sprintf('density_r.%s.%s.png', group_disp, cutoff)),
           p, width = 6, height = 5.2, dpi = 300)
  }
}

# ---------------------- 2. difference between the groups ----------------------
df_diff <- read_grid('diff')
min_density <- min(df_diff[[var_z]])
max_density_diff <- max(df_diff[[var_z]])

# Colour scale: 0 is fixed to white, the negative side blue and the positive side red.
# There are 10 colours, so there must be exactly 10 values, or the mapping shifts.
diverging_cols <- c('#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0',  # 5 negative
                    'white',                                                # 0
                    '#fddbc7', '#f4a582', '#d6604d', '#b2182b')             # 4 positive

# The break positions keep the relative structure that was hand-tuned on the original
# figure (as fractions of the extremes) but multiply it by this run's actual extremes, so
# a different sample cannot shift them. The original script hard-coded numbers such as
# -10.05 / 1.47.
neg_frac <- c(1.00, 0.70, 0.45, 0.25, 0.08)   # relative to min_density
pos_frac <- c(0.27, 0.51, 0.75, 1.00)         # relative to max_density_diff
col_values <- scales::rescale(
  c(min_density * neg_frac, 0, max_density_diff * pos_frac),
  from = c(min_density, max_density_diff))
stopifnot(length(col_values) == length(diverging_cols))

# Legend ticks: a few evenly spaced steps on the negative side, plus 0 and the positive
# maximum
breaks <- sort(unique(round(c(seq(min_density, 0, length.out = 5),
                             max_density_diff), 2)))

p_diff <- base_plot(df_diff, sprintf('%s - %s', group_labels[['Group1']],
                                              group_labels[['Group2']])) +
  scale_fill_gradientn(
    colours = diverging_cols,
    values  = col_values,
    limits  = c(min_density, max_density_diff),
    breaks  = breaks,
    name    = 'Density'
  )

ggsave(file.path(fig_path, 'density_r_diff.png'), p_diff,
       width = 6, height = 5.2, dpi = 300)

cat('done. figures ->', fig_path, '\n')
