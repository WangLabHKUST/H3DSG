# =============================================================================
# hox_dotplots / dotplot_df.R -- DotPlotDF(): Seurat::DotPlot that takes a data frame
#
# Seurat::DotPlot does two things: (1) computes avg.exp / pct.exp from the object, and
# (2) draws the dot plot. The statistics of this module are computed beforehand by
# calculate_*_gene_set_stats_advanced() (hox_helpers.R) and look like this:
#
#   features.plot        id      avg.exp     pct.exp   orig.ident
#   PAX3          malignant  0.394172865  40.5503616   PL
#   NKX6-1        malignant  0.025936205   3.0823169   PL
#
# so only step (2) of DotPlot is reproduced here, with a data frame as the input.
# Correspondence with Seurat:
#   * default layout as in Seurat: x = gene / gene set, y = group. The figures of this
#     module are transposed: use flip = TRUE.
#   * scale = TRUE  -> z-score of log1p(avg.exp) across groups for each feature, then
#                      clipped to [col.min, col.max]
#     scale = FALSE -> log1p(avg.exp) as it is (this is what Seurat does)
#     log.transform = FALSE -> not even log1p: the colour is the raw avg.exp
#   * pct.exp can be on a 0-100 scale (the scale of hox_helpers.R) or 0-1 (Seurat's
#     internal scale); it is detected automatically or set with pct.is.percent.
#     dot.min is always compared on the 0-1 scale (as in Seurat).
#   * facet.by is an addition of this module (Seurat has none): one panel per sample / group.
#
# Returns a ggplot object; attr(p, "data.plot") is the table that was drawn.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# Local copy of Seurat::MinMax, so that clipping does not need Seurat
.dp_minmax <- function(x, min, max) {
  x[x > max] <- max
  x[x < min] <- min
  x
}

# Display order of a column: a factor keeps its levels; otherwise the order of appearance
# (sort.levels = TRUE: alphabetical, as ggplot2 does for a character column)
.dp_levels <- function(x, sort.levels = FALSE) {
  if (is.factor(x)) return(levels(droplevels(x)))
  u <- unique(as.character(x))
  if (sort.levels) sort(u) else u
}

.dp_as_factor <- function(x, sort.levels = FALSE) {
  factor(as.character(x), levels = .dp_levels(x, sort.levels))
}

DotPlotDF <- function(
  data,
  features = NULL,             # character vector: subset and order; named list: also gives feature.groups
  group.by = "id",             # column used as the group (Seurat's ident)
  facet.by = NULL,             # column used for panels (not in Seurat)
  split.by = NULL,             # as in Seurat: split each group again, colours per split given by cols
  idents = NULL,               # draw only these values of group.by
  feature.col = "features.plot",
  avg.col = "avg.exp",
  pct.col = "pct.exp",
  cols = c("lightgrey", "blue"),
  col.min = -2.5,             # clipping range of the colour value. scale = TRUE: clips the z-score (as in Seurat);
  col.max = 2.5,              # scale = FALSE: clips avg.exp only if one of them is given explicitly
  dot.min = 0,
  dot.scale = 6,               # scalar = Seurat's range c(0, dot.scale); length 2 = the range itself
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA,
  log.transform = TRUE,        # FALSE = colour is the raw avg.exp
  scale.within.facet = FALSE,  # TRUE = z-score within each panel, FALSE = across all panels
  pct.is.percent = NA,         # NA = detect whether the pct column is 0-100 or 0-1
  flip = FALSE,                # TRUE = groups on x, features on y
  sort.levels = FALSE,         # order of a non-factor column: FALSE = order of appearance, TRUE = alphabetical
  feature.groups = NULL,       # named vector feature -> group name; created when features is a list
  facet.nrow = NULL,
  facet.ncol = NULL,
  facet.drop = TRUE,
  facet.scales = "fixed",      # scales of facet_wrap; "free_x" = each panel has its own x ticks
  color.scale = NULL,          # a ggplot2 colour scale; when given, cols is not used
  color.limits = NULL,
  color.oob = scales::squish,
  color.name = NULL,
  size.name = "Percent Expressed",
  size.breaks = NULL,          # breaks of the dot size legend, e.g. seq(0, 60, 20); NULL = chosen by ggplot
  size.labels = NULL,          # labels of the breaks; NULL = the break values
  size.limits = NULL,          # range of the dot size scale, overrides c(scale.min, scale.max);
                               # c(0, NA) makes the 0 break show up
  x.text.angle = 45,
  base.theme = ggplot2::theme_classic()
) {
  stopifnot(is.data.frame(data))
  # with scale = FALSE the colour is not clipped (as in Seurat) unless the caller gives col.min / col.max
  clamp.unscaled <- !missing(col.min) || !missing(col.max)
  need <- c(feature.col, avg.col, pct.col, group.by)
  miss <- setdiff(need, colnames(data))
  if (length(miss)) {
    stop("data lacks these columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!is.null(facet.by) && !facet.by %in% colnames(data)) {
    stop("facet.by column does not exist: ", facet.by, call. = FALSE)
  }
  if (!is.null(split.by) && !split.by %in% colnames(data)) {
    stop("split.by column does not exist: ", split.by, call. = FALSE)
  }
  scale.func <- switch(
    EXPR = scale.by,
    'size'   = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )

  # ---- features: subset / order / feature.groups ------------------------------
  if (is.list(features) ||
      (!is.null(features) && !is.null(names(features)) && any(!is.na(names(features))))) {
    grp <- unlist(lapply(seq_along(features), function(i) {
      rep(names(features)[i], each = length(features[[i]]))
    }))
    if (any(is.na(grp))) {
      warning("Some feature groups are unnamed.", call. = FALSE, immediate. = TRUE)
    }
    features <- unlist(features, use.names = FALSE)
    names(grp) <- features
    if (is.null(feature.groups)) feature.groups <- grp
  }

  df <- as.data.frame(data, stringsAsFactors = FALSE)
  feature.levels <- if (is.null(features)) .dp_levels(df[[feature.col]], sort.levels) else as.character(features)
  df <- df[as.character(df[[feature.col]]) %in% feature.levels, , drop = FALSE]
  if (!nrow(df)) stop("no row is left after filtering by features.", call. = FALSE)

  # ---- groups ---------------------------------------------------------------
  id.levels <- .dp_levels(df[[group.by]], sort.levels)
  if (!is.null(idents)) {
    df <- df[as.character(df[[group.by]]) %in% as.character(idents), , drop = FALSE]
    id.levels <- id.levels[id.levels %in% as.character(idents)]
    if (!nrow(df)) stop("no row is left after filtering by idents.", call. = FALSE)
  }
  df$id <- as.character(df[[group.by]])

  split.colors <- FALSE
  if (!is.null(split.by)) {
    splits <- as.character(df[[split.by]])
    unique.splits <- unique(splits)
    split.colors <- !any(cols %in% rownames(RColorBrewer::brewer.pal.info))
    if (split.colors) {
      if (length(unique.splits) > length(cols)) {
        stop("Need to specify at least ", length(unique.splits), " colors using the cols parameter")
      }
      cols <- cols[seq_along(unique.splits)]
      names(cols) <- unique.splits
    }
    df$id <- paste(df$id, splits, sep = "_")
    id.levels <- paste0(rep(id.levels, each = length(unique.splits)), "_",
                        rep(unique.splits, times = length(id.levels)))
    id.levels <- id.levels[id.levels %in% unique(df$id)]
  }

  # ---- pct: drawn on 0-100; dot.min is compared on 0-1 (as in Seurat) -------------
  pct <- as.numeric(df[[pct.col]])
  is.percent <- if (is.na(pct.is.percent)) {
    isTRUE(max(pct, na.rm = TRUE) > 1)
  } else {
    isTRUE(pct.is.percent)
  }
  pct.frac <- if (is.percent) pct / 100 else pct
  pct.frac[pct.frac < dot.min] <- NA_real_

  data.plot <- data.frame(
    features.plot = as.character(df[[feature.col]]),
    id            = df$id,
    avg.exp       = as.numeric(df[[avg.col]]),
    pct.exp       = pct.frac * 100,
    stringsAsFactors = FALSE
  )
  if (!is.null(facet.by)) data.plot$facet <- df[[facet.by]]

  # ---- cluster.idents: hierarchical clustering of the groups by (avg.exp, pct.exp) ---
  if (cluster.idents) {
    wide <- reshape(
      data.plot[, c("id", "features.plot", "avg.exp", "pct.exp")],
      idvar = "id", timevar = "features.plot", direction = "wide"
    )
    mat <- as.matrix(wide[, setdiff(colnames(wide), "id"), drop = FALSE])
    rownames(mat) <- wide$id
    mat[is.na(mat)] <- 0
    mat <- scale(mat)
    mat[!is.finite(mat)] <- 0
    id.levels <- rownames(mat)[hclust(d = dist(x = mat))$order]
  }

  # ---- colour: z-score(log1p) / log1p / raw ---------------------------------------
  ngroup <- length(unique(data.plot$id))
  if (ngroup == 1 && scale) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled",
            call. = FALSE, immediate. = TRUE)
  } else if (ngroup < 5 && scale) {
    warning("Scaling data with a low number of groups may produce misleading results",
            call. = FALSE, immediate. = TRUE)
  }
  trans <- if (log.transform) log1p else identity
  key <- if (scale.within.facet && !is.null(facet.by)) {
    paste(data.plot$features.plot, as.character(data.plot$facet), sep = "\r")
  } else {
    data.plot$features.plot
  }
  data.plot$avg.exp.scaled <- ave(trans(data.plot$avg.exp), factor(key), FUN = function(v) {
    if (!scale) {
      if (clamp.unscaled) v <- .dp_minmax(v, min = col.min, max = col.max)
      return(v)
    }
    s <- as.numeric(scale(v))
    s[!is.finite(s)] <- 0   # a feature that is identical in all groups has sd = 0: use 0, not NaN
    .dp_minmax(s, min = col.min, max = col.max)
  })

  if (split.colors) {
    binned <- as.numeric(cut(data.plot$avg.exp.scaled, breaks = 20))
    splits.use <- sub(
      paste0(".*_(", paste(sort(unique.splits, decreasing = TRUE), collapse = "|"), ")$"),
      "\\1", data.plot$id
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) colorRampPalette(colors = c("grey", color))(20)[value],
      color = cols[splits.use], value = binned
    )
  }

  # ---- axes, dot size -------------------------------------------------------------
  data.plot$features.plot <- factor(data.plot$features.plot, levels = feature.levels)
  data.plot$id <- factor(data.plot$id, levels = id.levels)
  if (!is.null(facet.by)) data.plot$facet <- .dp_as_factor(df[[facet.by]], sort.levels)
  if (!is.na(scale.min)) data.plot$pct.exp[data.plot$pct.exp < scale.min] <- scale.min
  if (!is.na(scale.max)) data.plot$pct.exp[data.plot$pct.exp > scale.max] <- scale.max
  if (!is.null(feature.groups)) {
    data.plot$feature.groups <- factor(
      feature.groups[as.character(data.plot$features.plot)],
      levels = unique(feature.groups)
    )
  }
  data.plot <- data.plot[order(data.plot$id, data.plot$features.plot), , drop = FALSE]

  # ---- draw -----------------------------------------------------------------------
  color.by <- if (split.colors) "colors" else "avg.exp.scaled"
  aes.x <- if (flip) "id" else "features.plot"
  aes.y <- if (flip) "features.plot" else "id"
  if (is.null(color.name)) {
    color.name <- if (scale) {
      "Scaled Expression"
    } else if (log.transform) {
      "Expression"
    } else {
      "Average Expression"
    }
  }

  plot <- ggplot(data.plot, aes(x = .data[[aes.x]], y = .data[[aes.y]])) +
    geom_point(aes(size = .data[["pct.exp"]], color = .data[[color.by]])) +
    scale.func(range = if (length(dot.scale) == 2) dot.scale else c(0, dot.scale),
               limits = if (is.null(size.limits)) c(scale.min, scale.max) else size.limits,
               breaks = if (is.null(size.breaks)) waiver() else size.breaks,
               labels = if (is.null(size.labels)) waiver() else size.labels,
               name = size.name) +
    base.theme +
    theme(axis.text.x = element_text(angle = x.text.angle,
                                     hjust = if (x.text.angle == 0) 0.5 else 1)) +
    labs(x = NULL, y = NULL)

  # facets: facet.by and feature.groups may both be given; then they are crossed
  has.fg <- !is.null(feature.groups)
  if (!is.null(facet.by) && has.fg) {
    plot <- plot + if (flip) {
      facet_grid(feature.groups ~ facet, scales = "free_y", space = "free_y", drop = facet.drop)
    } else {
      facet_grid(facet ~ feature.groups, scales = "free_x", space = "free_x", drop = facet.drop)
    }
  } else if (!is.null(facet.by)) {
    plot <- plot + facet_wrap(~facet, nrow = facet.nrow, ncol = facet.ncol,
                              drop = facet.drop, scales = facet.scales)
  } else if (has.fg) {
    plot <- plot + if (flip) {
      facet_grid(feature.groups ~ ., scales = "free_y", space = "free_y", switch = "y")
    } else {
      facet_grid(~feature.groups, scales = "free_x", space = "free_x", switch = "y")
    }
    plot <- plot + theme(panel.spacing = unit(1, "lines"), strip.background = element_blank())
  }

  if (split.colors) {
    plot <- plot + scale_color_identity()
  } else if (!is.null(color.scale)) {
    plot <- plot + color.scale
  } else if (length(cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols, name = color.name,
                                         limits = color.limits, oob = color.oob)
  } else if (length(cols) == 2) {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2], name = color.name,
                                        limits = color.limits, oob = color.oob)
  } else {
    plot <- plot + scale_color_gradientn(colors = cols, name = color.name,
                                         limits = color.limits, oob = color.oob)
  }
  attr(plot, "data.plot") <- data.plot
  plot
}
