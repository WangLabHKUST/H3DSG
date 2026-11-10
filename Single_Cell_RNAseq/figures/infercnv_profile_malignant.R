# =============================================================================
# Copy number profile of the malignant cells (inferCNV)
#
# Input : INFERCNV_RUN/run.final.infercnv_obj   the inferCNV run of pipeline/06_infercnv
#         (environment variable INFERCNV_RUN, default ../results/infercnv)
#         ../sample_metadata.tsv                 sample names and the order of the sampling sites
# Output: infercnv_profile_malignant.png in this directory (600 dpi, or INFERCNV_PNG_RES)
#
# EXCEPTION to the rule of this directory: this script does not read ../data/. The residual
# expression of the 9,252 genes in the 83,965 cells drawn is 780 million values, far too large
# for ../data/. The inferCNV run is not re-run here: the script only reads its final object
# and draws it. Source_Data_utils/code/single_cell_infercnv_profile.R exports the drawn values.
#
# What is drawn: the residual expression relative to the normal reference (1 = no change; red =
# gain, blue = loss) of the 81,304 malignant cells, one block per sample, and of the 2,661
# reference cells (oligodendrocytes and myeloid cells) on top. Genes run along the genome,
# chromosomes 1 to 22 left to right. Cells are not clustered within a sample. Samples are
# stacked from the bottom up in the order of their renamed_id_v2 as a text string
# (PL_1, PL_10, PL_11, PL_2 ... TL_6_5), the order of the submitted figure. The side bar is
# coloured by patient, the shades running along the sampling sites of a patient. The colour
# scale is fixed to 0.9155 - 1.0845, the range inferCNV chose for this run.
#
# inferCNV's own drawing functions are used with four changes: black lines between the samples
# of the malignant block and none in the reference block, no dendrogram, patient colours for
# the malignant side bar and one colour for the reference. The patches are made in memory and
# undone afterwards.
# =============================================================================

suppressPackageStartupMessages(library(infercnv))

RUN_DIR <- Sys.getenv("INFERCNV_RUN", "../results/infercnv")
PNG_RES <- as.integer(Sys.getenv("INFERCNV_PNG_RES", "600"))
X_RANGE <- c(0.915549045197481, 1.08445095480252)
REF_COLOR <- "#d0faea"

meta <- read.delim("../sample_metadata.tsv", colClasses = "character")
obj  <- readRDS(file.path(RUN_DIR, "run.final.infercnv_obj"))
if (!is.matrix(obj@expr.data)) obj@expr.data <- as.matrix(obj@expr.data)
message("loaded ", nrow(obj@expr.data), " genes x ", ncol(obj@expr.data), " cells")

# ---- colours: one hue per patient, shades along the sampling sites ----------------------
stems <- names(sort(tapply(as.integer(meta$location_sort), meta$stem, min)))
base_cols <- setNames(colorRampPalette(RColorBrewer::brewer.pal(min(length(stems), 12L), "Set3"))(length(stems)), stems)
lighten <- function(col, frac = 0.5) colorRampPalette(c("white", col))(100L)[max(1L, round(100L * (1 - frac)))]
darken  <- function(col, frac = 0.4) colorRampPalette(c(col, "black"))(100L)[max(1L, round(100L * frac))]
obs_colors <- character(0)
for (st in stems) {
  rows <- meta[meta$stem == st, ]
  rows <- rows[order(as.integer(rows$location_sort)), ]
  b <- base_cols[[st]]
  cols <- if (nrow(rows) == 1L) b else colorRampPalette(c(lighten(b, 0.55), b, darken(b, 0.35)))(nrow(rows))
  obs_colors <- c(obs_colors, setNames(cols, rows$renamed_id))
}
assign(".infercnv_custom_obs_colors", obs_colors, envir = .GlobalEnv)

# ---- malignant groups only, named by renamed_id, in the order of renamed_id_v2 ---------------
groups <- intersect(meta$sid, names(obj@observation_grouped_cell_indices))
stopifnot(length(groups) == 19)
v2 <- meta$renamed_id_v2[match(groups, meta$sid)]
groups <- groups[order(v2)]                                         # string order, bottom to top
new_names <- meta$renamed_id[match(groups, meta$sid)]
obj@observation_grouped_cell_indices <- setNames(obj@observation_grouped_cell_indices[groups], new_names)

# keep the cells that are drawn and rebuild the integer indices
all_names <- colnames(obj@expr.data)
ref_names <- unique(unlist(lapply(obj@reference_grouped_cell_indices, function(i) all_names[i])))
obs_names <- unique(unlist(lapply(obj@observation_grouped_cell_indices, function(i) all_names[i])))
obj@expr.data <- obj@expr.data[, all_names[all_names %in% union(ref_names, obs_names)], drop = FALSE]
new_cols <- colnames(obj@expr.data)
reidx <- function(l) lapply(l, function(i) as.integer(match(all_names[i], new_cols)))
obj@observation_grouped_cell_indices <- reidx(obj@observation_grouped_cell_indices)
obj@reference_grouped_cell_indices   <- reidx(obj@reference_grouped_cell_indices)
message("drawn: ", ncol(obj@expr.data), " cells (malignant ", length(obs_names), ", reference ", length(ref_names), ")")

# cells are not reordered within a sample: an identity tree per group
make_identity_hc <- function(cell_names) {
  n <- length(cell_names)
  if (n <= 1L) return(NULL)
  merge_rows <- vector("list", n - 1L); heights <- numeric(n - 1L); counter <- 0L
  build <- function(idxs) {
    m <- length(idxs)
    if (m == 1L) return(idxs)
    mid <- m %/% 2L
    left <- build(idxs[seq_len(mid)]); right <- build(idxs[seq.int(mid + 1L, m)])
    counter <<- counter + 1L
    merge_rows[[counter]] <<- c(left, right)
    heights[counter] <<- max(if (left > 0L) heights[left] else 0L, if (right > 0L) heights[right] else 0L) + 1L
    counter
  }
  build(-seq_len(n))
  structure(list(merge = do.call(rbind, lapply(merge_rows, as.integer)), height = heights, order = seq_len(n),
                 labels = cell_names, method = "ward.D", call = NULL, dist.method = "euclidean"), class = "hclust")
}
obj@tumor_subclusters$hc <- list(); obj@tumor_subclusters$subclusters <- list()
for (g in names(obj@observation_grouped_cell_indices)) {
  idx <- obj@observation_grouped_cell_indices[[g]]
  obj@tumor_subclusters$hc[[g]] <- make_identity_hc(new_cols[idx])
  obj@tumor_subclusters$subclusters[[g]] <- setNames(list(setNames(idx, new_cols[idx])), g)
}

# ---- draw ---------------------------------------------------------------------------------
plot_profile <- function() {
  patch_body <- function(fn_name, pat, rep) {
    fn  <- get(fn_name, envir = asNamespace("infercnv"))
    txt <- deparse(body(fn), width.cutoff = 500L)
    if (!any(grepl(pat, txt, perl = TRUE))) stop("patch does not match: ", fn_name, " / ", pat)
    body(fn) <- parse(text = paste(gsub(pat, rep, txt, perl = TRUE), collapse = "\n"))[[1]]
    utils::assignInNamespace(fn_name, fn, ns = "infercnv")
  }
  has <- function(fn_name, pat)
    any(grepl(pat, deparse(body(get(fn_name, envir = asNamespace("infercnv"))), width.cutoff = 500L), perl = TRUE))
  orig_obs <- get(".plot_cnv_observations", envir = asNamespace("infercnv"))
  orig_ref <- get(".plot_cnv_references",   envir = asNamespace("infercnv"))
  on.exit({ utils::assignInNamespace(".plot_cnv_observations", orig_obs, ns = "infercnv")
            utils::assignInNamespace(".plot_cnv_references",   orig_ref, ns = "infercnv") }, add = TRUE)

  if (has(".plot_cnv_observations", "row_seps\\s*=\\s*c\\(\\)"))        # lines between the samples
    patch_body(".plot_cnv_observations", "row_seps\\s*=\\s*c\\(\\)", "row_seps = obs_seps")
  if (has(".plot_cnv_references", "row_seps\\s*=\\s*ref_seps"))         # none in the reference block
    patch_body(".plot_cnv_references", "row_seps\\s*=\\s*ref_seps", "row_seps = c()")
  patch_body(".plot_cnv_observations", 'dendrogram\\s*=\\s*"row"', 'dendrogram = "none"')
  patch_body(".plot_cnv_observations",
    "row_groupings\\s*<-\\s*cbind\\(row_groupings,\\s*get_group_color_palette.*",
    paste0('row_groupings <- cbind(row_groupings, { .tmp_oc <- get(".infercnv_custom_obs_colors", envir = .GlobalEnv); ',
           '.tmp_oc[obs_annotations_names][hcl_obs_annotations_groups] })'))
  patch_body(".plot_cnv_observations",
    "annotations_legend\\s*<-\\s*cbind\\(obs_annotations_names,\\s*get_group_color_palette.*",
    paste0('annotations_legend <- cbind(obs_annotations_names, { .tmp_oc2 <- get(".infercnv_custom_obs_colors", envir = .GlobalEnv); ',
           '.tmp_oc2[obs_annotations_names] })'))
  patch_body(".plot_cnv_references", "row_groupings\\s*<-\\s*as\\.matrix\\(get_group_color_palette.*",
             sprintf('row_groupings <- as.matrix(rep("%s", length(split_groups)))', REF_COLOR))
  patch_body(".plot_cnv_references", "annotations_legend\\s*<-\\s*cbind\\(name_ref_groups.*",
             sprintf('annotations_legend <- matrix(c("Normal", "%s"), nrow = 1L, ncol = 2L)', REF_COLOR))

  out <- tempfile("infercnv_plot_"); dir.create(out)
  infercnv::plot_cnv(obj, out_dir = out, output_filename = "infercnv_profile_malignant",
                     title = "inferCNV: malignant cells by sample", x.center = 1, x.range = X_RANGE,
                     cluster_by_groups = TRUE, cluster_references = FALSE, hclust_method = "ward.D",
                     output_format = "png", png_res = PNG_RES, color_safe_pal = FALSE, useRaster = TRUE)
  file.copy(file.path(out, "infercnv_profile_malignant.png"), "infercnv_profile_malignant.png", overwrite = TRUE)
}
plot_profile()
