# =============================================================================
# Metaprograms: aggregate the cNMF programs of the 19 samples (plateau_v1_all19)
#
# Input : ../../results/malignant_states/03_cnmf/   top-100 genes of every program
#             (<sid>.programs.genes100.k_<K>.csv) and cell usages (<sid>/<sid>.usages.k_<K>.
#             dt_0_02.consensus.txt) of step 04
#         ../../results/malignant_states/selected_K_plateau.txt   K per sample (step 03b)
#         ../../objects/per_sample/<sid>_minf400.rds   S / G2M scores and QC columns, read only
#             for samples missing from the cache cell_qc_cellcycle.csv
#         ../../sample_metadata.tsv
# Output: ../../results/malignant_states/04_metaprograms/plateau_v1_all19/
#             cellular_states_metaprograms.txt   the 50 signature genes of M1..M6 (main result)
#             program_metaprogram_assignment.csv, program_similarity_matrix.csv,
#             program_vs_qcmetric_correlation.csv, program_similarity_heatmap_data.rds
#             program_similarity_heatmap.pdf (all programs, then the metaprograms),
#             metaprogram_similarity_heatmap.png (the second page alone)
#         ../../results/malignant_states/04_metaprograms/cell_qc_cellcycle.csv   cache, written once
#
# Outline
#   (a) top 100 genes of every program; every cell goes to the program with the highest
#       usage; programs with < 10 dominant cells or < 1% of the sample are dropped
#   (b) program similarity = size of the intersection of the two top-100 gene sets (0-100).
#       It is not a correlation and not a Jaccard index
#   (c) hierarchical clustering (pheatmap default: complete linkage, euclidean distance)
#   (d) recursive cut of the tree (define_metaprograms): first cut into 3 subtrees; a subtree
#       with < 6 programs is dropped (outlier); a subtree with a mean pairwise overlap >= 10
#       genes is one metaprogram; any other subtree is cut in two and examined again.
#       Numbered M1..Mn from left to right in the heatmap
#   (e) signature: how often each gene is in the top 100 of the programs of a metaprogram,
#       without ^RP[LS], ^A[CL]<6 digits> and ^MT- genes; the top 50
#   (f) correlation of the program usages with QC metrics and S / G2M scores, to look for
#       technical artefacts
#
# M1..M6 are metaprograms, unrelated to the numbers of the steps. They are named in
# 06c_compare_reference.R and used in 07_score_states_and_plot.R.
# Runtime about 5 min (the first run computes the cell-cycle scores). R 4.5.1.
# Paths are relative to this directory: run the script from here.
# =============================================================================

PER_SAMPLE_DIR   <- "../../objects/per_sample"
SAMPLE_META_FILE <- "../../sample_metadata.tsv"
RES_DIR     <- "../../results/malignant_states"
DIR_03_CNMF <- file.path(RES_DIR, "03_cnmf")
KFILE       <- file.path(RES_DIR, "selected_K_plateau.txt")
DIR_04_MP   <- file.path(RES_DIR, "04_metaprograms")
MP_VERSION  <- "plateau_v1_all19"
DIR_MP      <- file.path(DIR_04_MP, MP_VERSION)
CC_CACHE    <- file.path(DIR_04_MP, "cell_qc_cellcycle.csv")

config <- list(
  seed              = 42,
  program_min_cells = 10,     # programs with fewer dominant cells are dropped
  program_min_prop  = 0.01,   # programs with a smaller share of dominant cells are dropped
  mp_first_cut      = 3,      # number of subtrees of the first cut of the tree
  mp_min_similarity = 10,     # lower bound of the mean top-100 overlap inside a subtree
  mp_min_programs   = 6,      # fewest programs of a metaprogram
  mp_n_genes        = 50      # signature genes per metaprogram
)

suppressPackageStartupMessages({ library(dplyr); library(pheatmap); library(dendextend); library(purrr) })
dir.create(DIR_MP, recursive = TRUE, showWarnings = FALSE)

SAMPLE_META <- read.delim(SAMPLE_META_FILE, stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
SIDS <- SAMPLE_META$sid

#' sid -> renamed_id (the name used outside the pipeline)
to_renamed <- function(sid) {
  out <- SAMPLE_META$renamed_id[match(as.character(sid), SAMPLE_META$sid)]
  ifelse(is.na(out), as.character(sid), out)
}

#' Seurat v5 objects with several layers: join them back into one
join_layers_if_needed <- function(obj, assay = "RNA") {
  joined <- try(SeuratObject::JoinLayers(obj, assay = assay), silent = TRUE)
  if (inherits(joined, "try-error")) obj else joined
}

# =============================================================================
# Similarity heatmap: drawing code
# program_similarity_heatmap_data.rds stores the inputs of the plots (not the plots), so a
# heatmap can be redrawn from the rds without rerunning the clustering.
# =============================================================================
custom_magma <- c(colorRampPalette(c("white", rev(viridis::magma(323, begin = 0.15))[1]))(10),
                  rev(viridis::magma(323, begin = 0.18)))

# Colours of the Sample annotation bar: the Seurat DimPlot default, scales::hue_pal() (evenly
# spaced hues, in the order the samples appear in the matrix = the order of selected_K)
sample_colors <- function(sids) {
  sids <- unique(as.character(sids))
  setNames(scales::hue_pal()(length(sids)), sids)
}

#' Pack what the plots need.
#' @param cnmf_intersect similarity matrix (size of the top-100 gene intersection)
#' @param ann            program annotation table, row names = program names
#' @param ordered_p      programs of the second plot, ordered by metaprogram
#' @param mp_idx         position of every metaprogram in the heatmap order, for the gaps
#' @param tree_row,tree_col clustering trees (hclust) of the first plot; stored so that a
#'        redraw does not depend on the pheatmap version recomputing distances / clustering
#' @param ann_cols       annotation columns drawn beside the heatmap (ann keeps the others)
heatmap_data <- function(cnmf_intersect, ann, ordered_p, mp_idx,
                         tree_row, tree_col, mp_list = NULL, hm_order = NULL,
                         ann_cols = c("Sample", "Metaprogram"),
                         ann_colors = NULL, meta = list()) {
  list(cnmf_intersect = cnmf_intersect,
       ann            = ann,
       ann_cols       = intersect(ann_cols, names(ann)),
       ann_colors     = ann_colors,
       ordered_p      = ordered_p,
       gaps           = cumsum(map_dbl(mp_idx, length)),
       tree_row       = tree_row,
       tree_col       = tree_col,
       mp_idx         = mp_idx,
       mp_list        = mp_list,
       hm_order       = hm_order,
       meta           = meta)
}

#' annotation columns + colours; only the colours of the columns that are drawn are kept,
#' otherwise pheatmap complains that annotation_colors do not match
.ann_args <- function(d) {
  ac   <- d$ann_cols
  acol <- d$ann_colors[intersect(names(d$ann_colors), ac)]
  list(cols = ac, colors = if (length(acol)) acol else NA)  # pheatmap's "unspecified" is NA, not an empty list
}

#' Page 1: all programs, in the order of the clustering tree
draw_program_heatmap <- function(d, silent = FALSE) {
  a <- .ann_args(d)
  pheatmap(d$cnmf_intersect, color = custom_magma, show_colnames = FALSE, border_color = NA,
           cluster_rows = d$tree_row, cluster_cols = d$tree_col,
           # ward.D2 would pull programs of low similarity together and cut metaprograms
           # apart, so the pheatmap default (complete + euclidean) is used; the tree is
           # computed in the main script and passed in
           fontsize_row = 5,
           annotation_col = d$ann[names(d$cnmf_intersect), a$cols, drop = FALSE],
           annotation_colors = a$colors, silent = silent,
           main = "programs clustered by top100-gene overlap")
}

#' Page 2: only the programs that entered a metaprogram, in blocks by metaprogram
draw_metaprogram_heatmap <- function(d, silent = FALSE) {
  a  <- .ann_args(d)
  op <- d$ordered_p
  pheatmap(d$cnmf_intersect[op, op], cluster_rows = FALSE, cluster_cols = FALSE,
           gaps_row = d$gaps, gaps_col = d$gaps,
           annotation_col = d$ann[op, a$cols, drop = FALSE],
           annotation_colors = a$colors, silent = silent,
           color = custom_magma, show_colnames = FALSE, border_color = NA, fontsize_row = 5,
           # the numeric annotation columns are the CORRELATION (-1..1) of the program usage with
           # that metric, not the metric itself
           main = "metaprograms (annotated with cell-cycle score correlations)")
}

#' The two pages of program_similarity_heatmap.pdf
write_similarity_heatmap_pdf <- function(d, file, width = 10, height = 9) {
  pdf(file, width = width, height = height)
  on.exit(dev.off())
  draw_program_heatmap(d)
  draw_metaprogram_heatmap(d)
  invisible(file)
}

#' Page 2 alone as a PNG; sizes in inches as for the PDF; bg = "white", not transparent
write_metaprogram_heatmap_png <- function(d, file, width = 10, height = 9, dpi = 300) {
  args <- list(filename = file, width = width * dpi, height = height * dpi,
               res = dpi, bg = "white")
  if (isTRUE(capabilities("cairo"))) args$type <- "cairo-png"   # anti-aliasing, if cairo is there
  do.call(png, args)
  on.exit(dev.off())
  draw_metaprogram_heatmap(d)
  invisible(file)
}

set.seed(config$seed)

ktab <- read.table(KFILE, header = FALSE, comment.char = "#", col.names = c("sample", "K"), stringsAsFactors = FALSE)
SELECTED_K <- setNames(as.integer(ktab$K), ktab$sample)
samples <- names(SELECTED_K)
stopifnot(setequal(samples, SIDS))
message("aggregating ", length(samples), " samples -> ", DIR_MP)
print(SELECTED_K)

# =============================================================================
# (a) programs of every sample; drop the ones without cells
# =============================================================================
extract_program <- function(sid) {
  k <- SELECTED_K[[sid]]
  gs <- read.delim(file.path(DIR_03_CNMF, sprintf("%s.programs.genes100.k_%d.csv", sid, k)),
                   sep = ",", row.names = 1, check.names = FALSE)
  stopifnot(ncol(gs) == k)
  names(gs) <- paste0(sid, "_P", names(gs))
  usage <- read.delim(file.path(DIR_03_CNMF, sid, sprintf("%s.usages.k_%d.dt_0_02.consensus.txt", sid, k)),
                      row.names = 1, check.names = FALSE)
  grp <- paste0(sid, "_P", apply(usage, 1, which.max))   # every cell goes to the program with the highest usage
  n <- table(grp); p <- prop.table(n)
  drop <- union(names(n)[n < config$program_min_cells], names(p)[p < config$program_min_prop])
  keep <- setdiff(names(gs), drop)
  message(sprintf("  %-10s K=%d  kept %d/%d%s", sid, k, length(keep), ncol(gs),
                  if (length(drop)) paste0("  (dropped: ", paste(drop, collapse = ","), ")") else ""))
  gs[, keep, drop = FALSE]
}
cnmf_top_genes <- do.call(cbind, lapply(samples, extract_program))
message(ncol(cnmf_top_genes), " programs enter the aggregation")

# =============================================================================
# (b)(c) similarity matrix = size of the intersection of the top 100 genes
# =============================================================================
cnmf_intersect <- as.data.frame(
  sapply(cnmf_top_genes, function(x) sapply(cnmf_top_genes, function(y) length(intersect(x, y)))))
write.csv(cnmf_intersect, file.path(DIR_MP, "program_similarity_matrix.csv"))
ph0 <- pheatmap(cnmf_intersect, color = custom_magma, show_colnames = FALSE,
                border_color = NA, fontsize_row = 5, silent = TRUE)

# =============================================================================
# (d) define the metaprograms by recursive cutting of the tree
# =============================================================================
#' corrected version of dendextend::get_subdendrograms (the original takes the labels wrongly)
get_subdendrograms2 <- function(dend, k, ...) {
  clusters <- cutree(dend, k, ...)
  dend_list <- lapply(unique(clusters), function(cluster.id) find_dendrogram(dend, names(clusters)[which(clusters == cluster.id)]))
  class(dend_list) <- "dendlist"
  dend_list
}
#' mean pairwise similarity of the programs of a subtree (diagonal excluded)
avg_similarity <- function(dend) {
  x <- cnmf_intersect[labels(dend), labels(dend)]
  diag(x) <- NA
  mean(as.matrix(x), na.rm = TRUE)
}
define_metaprograms <- function(tree, K, min_similarity, min_programs) {
  recurse <- function(subtree, min_programs, min_similarity) {
    if (attr(subtree, "members") < min_programs) {
      return(NULL)
    } else if (avg_similarity(subtree) >= min_similarity) {
      metaprograms[[i]] <<- labels(subtree)
      i <<- i + 1
      return(labels(subtree))
    } else {
      lapply(get_subdendrograms2(subtree, 2), recurse, min_programs = min_programs, min_similarity = min_similarity)
    }
  }
  i <- 1; metaprograms <- list()
  lapply(get_subdendrograms2(tree, K), recurse, min_programs = min_programs, min_similarity = min_similarity)
  names(metaprograms) <- seq_along(metaprograms)
  metaprograms
}
mp_list <- define_metaprograms(as.dendrogram(ph0$tree_col), K = config$mp_first_cut,
                               min_similarity = config$mp_min_similarity, min_programs = config$mp_min_programs)
message(length(mp_list), " metaprograms, containing ", length(unlist(mp_list)), "/", ncol(cnmf_top_genes), " programs")
if (length(mp_list) == 0) stop("No metaprogram was found")

# number M1..Mn by position from left to right in the heatmap
hm_order <- ph0$tree_col$labels[ph0$tree_col$order]
mp_idx   <- map(mp_list, ~ which(hm_order %in% .x))
mp_list  <- mp_list[names(sort(map_dbl(mp_idx, 1)))]
mp_idx   <- mp_idx[names(sort(map_dbl(mp_idx, 1)))]
names(mp_list) <- names(mp_idx) <- seq_along(mp_idx)
mp_names <- paste0("M", seq_along(mp_list))

ann <- data.frame(Program = names(cnmf_intersect), stringsAsFactors = FALSE)
rownames(ann) <- ann$Program
ann$Sample <- sub("_P[0-9]+$", "", ann$Program)
ann$Sample_renamed  <- to_renamed(ann$Sample)
ann$Program_renamed <- paste0(ann$Sample_renamed, sub("^.*(_P[0-9]+)$", "\\1", ann$Program))
ann$Metaprogram <- "outlier"
for (j in seq_along(mp_list)) ann$Metaprogram[ann$Program %in% mp_list[[j]]] <- mp_names[j]

# =============================================================================
# (f) QC sanity check: program usage vs QC metrics
# =============================================================================
QC_COLS <- c("S.Score", "G2M.Score", "nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo")
QC_META_COLS <- c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo")

# Everything comes from the per-sample minf400 objects and is cached in CC_CACHE:
# S / G2M scores are computed (CellCycleScoring); nCount / nFeature / mito / ribo are read from
# the object metadata.
cellcycle_scores <- function(sids) {
  need <- c("sid", "cell", "S.Score", "G2M.Score", QC_META_COLS)
  cached <- if (file.exists(CC_CACHE)) read.csv(CC_CACHE, stringsAsFactors = FALSE, check.names = FALSE)[, need] else NULL
  todo <- setdiff(sids, unique(cached$sid))
  if (!length(todo)) return(cached)
  message("  computing S/G2M scores and reading QC metrics: ", length(todo), " samples")
  new <- do.call(rbind, lapply(todo, function(sid) {
    obj <- join_layers_if_needed(readRDS(file.path(PER_SAMPLE_DIR, paste0(sid, "_minf400.rds"))))
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    sc <- Seurat::CellCycleScoring(obj, s.features = Seurat::cc.genes.updated.2019$s.genes,
                                   g2m.features = Seurat::cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
    md <- sc@meta.data
    data.frame(sid = sid, cell = rownames(md), S.Score = md$S.Score, G2M.Score = md$G2M.Score,
               md[, QC_META_COLS], stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
  }))
  cached <- rbind(cached, new)
  write.csv(cached, CC_CACHE, row.names = FALSE)
  cached
}
CC_SCORES <- cellcycle_scores(samples)

qc_of <- function(sid) {
  cc <- CC_SCORES[CC_SCORES$sid == sid, , drop = FALSE]
  m <- cc[, c(QC_META_COLS, "S.Score", "G2M.Score"), drop = FALSE]
  rownames(m) <- cc$cell                          # cell names of the per-sample object, the same as in usage (no sid_ prefix)
  m
}
qc_corr <- function(sid) {
  mpw <- read.delim(file.path(DIR_03_CNMF, sid, sprintf("%s.usages.k_%d.dt_0_02.consensus.txt", sid, SELECTED_K[[sid]])),
                    row.names = 1, check.names = FALSE)
  names(mpw) <- paste0(sid, "_P", names(mpw))
  mqc <- qc_of(sid)
  ov <- intersect(rownames(mqc), rownames(mpw))
  if (length(ov) < 10) { message("  ", sid, ": QC metrics do not match the cells of the usage table, skipped"); return(NULL) }
  mqc <- mqc[, intersect(QC_COLS, names(mqc)), drop = FALSE]
  as.data.frame(t(cor(mqc[ov, ], mpw[ov, ], use = "pairwise.complete.obs"))) %>% tibble::rownames_to_column("Program")
}
cor_qcpw <- do.call(rbind, lapply(samples, qc_corr))
write.csv(cor_qcpw, file.path(DIR_MP, "program_vs_qcmetric_correlation.csv"), row.names = FALSE)
ann <- merge(ann, cor_qcpw, by = "Program", all.x = TRUE)
rownames(ann) <- ann$Program
write.csv(ann, file.path(DIR_MP, "program_metaprogram_assignment.csv"), row.names = FALSE)

# ---- similarity heatmap (save the plot inputs first, then draw) ---------------------------
ordered_p <- hm_order[unname(unlist(mp_idx))]
hm_data <- heatmap_data(
  cnmf_intersect = cnmf_intersect, ann = ann, ordered_p = ordered_p, mp_idx = mp_idx,
  tree_row = ph0$tree_row, tree_col = ph0$tree_col, mp_list = mp_list, hm_order = hm_order,
  ann_cols = c("Sample", "Metaprogram", "S.Score", "G2M.Score"),
  ann_colors = list(Sample = sample_colors(ann[names(cnmf_intersect), "Sample"])),
  meta = list(out_dir = DIR_MP, kfile = KFILE, samples = samples, K = SELECTED_K[samples], mp_names = mp_names,
              mp_config = config[grep("^mp_|^program_min_", names(config))], cnmf_dir = DIR_03_CNMF, created = Sys.time()))
saveRDS(hm_data, file.path(DIR_MP, "program_similarity_heatmap_data.rds"))
write_similarity_heatmap_pdf(hm_data, file.path(DIR_MP, "program_similarity_heatmap.pdf"), width = 10, height = 9)
write_metaprogram_heatmap_png(hm_data, file.path(DIR_MP, "metaprogram_similarity_heatmap.png"), width = 10, height = 9, dpi = 300)

# =============================================================================
# (e) signature genes of every metaprogram
# =============================================================================
signature_of <- function(programs, n = config$mp_n_genes) {
  tb <- table(unlist(cnmf_top_genes[, programs, drop = FALSE]))
  tb <- tb[!grepl("^RP[LS]", names(tb))]            # ribosomal protein genes
  tb <- tb[!grepl("^A[CL][0-9]{6}", names(tb))]      # unnamed clones / lncRNA
  tb <- tb[!grepl("^MT-", names(tb))]                # mitochondrial
  g <- names(sort(tb, decreasing = TRUE))
  length(g) <- n
  g
}
mpgs <- sapply(mp_list, signature_of)
colnames(mpgs) <- mp_names
write.table(mpgs, file.path(DIR_MP, "cellular_states_metaprograms.txt"), quote = FALSE, sep = "\t", row.names = FALSE, na = "")
message("signature genes -> ", file.path(DIR_MP, "cellular_states_metaprograms.txt"))
print(table(ann$Metaprogram))
