# =============================================================================
# Malignant cell selection and cNMF input export, one sample at a time (19 samples)
#
# Input : ../../objects/per_sample/<sid>_minf400.rds    per-sample objects; need the
#                                                       seurat_clusters column and a umap
#         ../../inputs/annotations/tumor_malignant_filter_overrides.csv
#                                                       manual cluster calls (sid, cluster, call)
#         ../../sample_metadata.tsv                     the 19 sample ids
# Output: ../../results/malignant_states/02_malignant/<sid>_filter_malignant.png
#                 cluster calls, kept cells and thresholds for each sample
#         ../../results/malignant_states/02_malignant/cluster_calls.csv, cluster_keep_rate.csv,
#                 auto_filter_summary.tsv, logs/<sid>.log
#                 (cluster_calls.csv is the table kept as 06_infercnv/malignant_cluster_calls.csv)
#         ../../results/malignant_states/03_cnmf_input/<sid>_filtered_malignant.tsv
#                 raw counts of the kept cells, cells in rows and genes in columns,
#                 tab separated: the cNMF input of step 03
#
# Outline (per sample, inside its own object, using the clusters stored in the object)
#   1. FindVariableFeatures (vst, 2000 genes) -> ScaleData (no regression) -> centroid of
#      every cluster on scale.data
#   2. classify_clusters() (utilities below) calls each cluster malig / normal / exclude from
#      marker genes; the manual overrides then replace the calls of the listed clusters
#      (TL_2: cluster 11 -> malig; TL_5: cluster 10 -> exclude)
#   3. maximum Pearson correlation of every cell with the normal centroids
#      (maxscore_tonormal) and with the malignant centroids (maxscore_tomalig)
#   4. T_NORMAL = 0.05; T_MALIG = Otsu threshold on maxscore_tomalig of the first malig cluster
#   5. keep the cells with maxscore_tonormal < T_NORMAL and maxscore_tomalig > T_MALIG that
#      sit in a malig cluster, and export their raw counts
#
# The samples run in parallel (PSOCK workers, one sample per worker, NCORES at most).
# Optional command line arguments restrict the run to the listed sample ids.
#
# Runtime about 8 min on 20 cores; each worker peaks at 4-8 GB. Seurat 5.2.1, R 4.5.1.
#
# Paths are relative to this directory: run the script from here.
# =============================================================================

PER_SAMPLE_DIR   <- "../../objects/per_sample"
OVERRIDE_CSV     <- "../../inputs/annotations/tumor_malignant_filter_overrides.csv"   # sid,cluster,call
SAMPLE_META_FILE <- "../../sample_metadata.tsv"
RES_DIR  <- "../../results/malignant_states"
OUT_DIR  <- file.path(RES_DIR, "02_malignant")
LOG_DIR  <- file.path(OUT_DIR, "logs")
SAVE_TO  <- file.path(RES_DIR, "03_cnmf_input")

T_NORMAL         <- 0.05   # the same in all samples, so the normal side is comparable across them
T_MALIG_FALLBACK <- 0.05   # used when Otsu fails (too few cells / degenerate values)
NUM_HVG  <- 2000
MIN_KEEP <- 50             # with fewer kept cells only the plot is made, no matrix is exported
NCORES   <- 20

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
  library(EBImage)       # otsu()
  library(foreach); library(doParallel)
})

dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SAVE_TO, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Utilities: call every cluster malig / normal / exclude from the DotPlot data
#
# The automatic call is only a starting point: classify_clusters() prints the three scores
# and the margin of every cluster, puts each cluster into the group with the highest score,
# and only marks weak calls in the `flag` column (calls$flagged); these must be checked by
# eye. The one rule that overrides "highest score" is malig_over_exclude (oligo-like
# malignant cells, see below), which writes oligo_like_malig into `flag`.
# =============================================================================

# Marker panel. The order is the DotPlot order: malig -> normal -> exclude
MARKER_SETS <- list(
  malig   = c("SOX2", "PTPRZ1"),
  normal  = c("PTPRC", "ITGAM", "AIF1", "P2RY12", "CD79A",
              "FCGR3B", "CD3D", "NKG7", "VWF", "MYL9"),
  # Oligodendrocytes: a group of their own, not counted as normal. A win of this group does
  # not always exclude the cluster: oligo-like tumour cells with equally high SOX2/PTPRZ1
  # are moved back to malig by malig_over_exclude
  exclude = c("PLP1", "MOG", "MAG")
)

# Plotted but not scored. Neither is tumour specific: APOE is as strong in microglia and
# macrophages, and MKI67 only means proliferation. Putting them in the malig group would pull
# myeloid and proliferating normal cells into potential_malig, but they stay in the plot as
# supporting evidence
PLOT_ONLY <- c("MKI67", "APOE")

# Flatten to the features vector of DotPlot, so that the gene sets and the plotted genes are
# defined once. PLOT_ONLY goes after malig to keep the order tumour -> normal -> oligodendrocyte
marker_genes <- function(sets = MARKER_SETS, plot_only = PLOT_ONLY)
  unique(c(sets$malig, plot_only,
           unname(unlist(sets[setdiff(names(sets), "malig")], use.names = FALSE))))

# Fixed colours, so that the 19 samples are coloured alike (with ggplot's default colours the
# colours would shift whenever a sample lacks one of the classes)
ANNO_COLORS <- c(potential_malig = "#B2182B",
                 obvious_normal  = "#2166AC",
                 excluded        = "#F0A202",   # oligodendrocytes: not normal, excluded on their own
                 ambiguous       = "#BDBDBD")   # not decided

# Spread the cluster-level calls over the cells, for colouring a DimPlot.
#   1. It takes the final three vectors, not the return value of classify_clusters, so the plot
#      made after the manual overrides shows what is actually used downstream.
#   2. Clusters in none of the three vectors become "ambiguous" rather than NA, because NA cells
#      turn into a grey NA group, and those are the cells that most need a look.
# Priority: excluded > obvious_normal > potential_malig
annotate_clusters <- function(clusters, obvious_normal, excluded, potential_malig) {
  clusters <- as.character(clusters)
  x <- rep("ambiguous", length(clusters))
  x[clusters %in% potential_malig] <- "potential_malig"
  x[clusters %in% obvious_normal]  <- "obvious_normal"
  x[clusters %in% excluded]        <- "excluded"
  droplevels(factor(x, levels = names(ANNO_COLORS)))
}

# Cluster names are strings such as '0'..'13'; a plain sort would put '10' before '2'
sort_clusters <- function(x) {
  x <- as.character(unique(x))
  n <- suppressWarnings(as.numeric(x))
  if (!anyNA(n)) x[order(n)] else sort(x)
}

# dot_plot   : the ggplot returned by DotPlot() (or its $data)
# use        : "scaled" uses avg.exp.scaled (z-score of a gene across clusters, comparable
#              scales; default); "raw" uses log1p(avg.exp) (highly expressed genes dominate)
# top_n      : per cluster and group, only the top_n most strongly expressed genes are scored
#              (default 3; NULL / Inf = all). This removes the bias from the size of the sets
#              (all groups contribute top_n genes to the sum), and a cluster does not have to
#              light up the whole panel (T cells have CD3D/NKG7 high but P2RY12/VWF low).
#              The first point needs every group to supply top_n genes: if a group has fewer
#              (malig has only SOX2/PTPRZ1) it adds one term less under sum; then either lower
#              top_n to the smallest set size or use agg_fun = "mean".
# min_pct    : genes with pct.exp below this (default 10, i.e. fewer than 10% of cells
#              expressing) are not scored. Without it, a gene that is zero almost everywhere
#              gets a large positive z-score from any small amount of expression (MYL9 with
#              avg.exp 0.048 in 5.6% of cells had z = +1.76, enough on its own to lift the
#              normal group above a tumour cluster with SOX2/PTPRZ1 in 95%+ of cells).
#              A group with all genes filtered out has no expression in that cluster and is
#              set to NA (= -Inf, cannot win). 0 switches the filter off.
# rank_by    : how the top_n genes are picked. "exp_pct" = avg.exp * pct.exp (default: big and
#              dark dots in the dotplot), "value" = the scored quantity itself. exp_pct picks by
#              A and scores by B, which biases the z-scores of the chosen genes upwards, but
#              equally for the three groups.
# agg_fun    : how to combine within a group, "sum" (default) or "mean". With top_n the two
#              differ by the constant top_n and rank the same.
# min_score  : a winning group with a score below this gets the low_score flag (scaled: 0, i.e.
#              the cluster must be above the average of all clusters on that marker set). It
#              targets false positives such as a sample without oligodendrocytes in which PLP1
#              still has a maximum z-score.
# min_margin : first minus second below this gets the low_margin flag. The default depends on
#              agg_fun: under sum it is top_n times the value under mean.
# force_call : TRUE (default) = every cluster is put into the group with the highest score;
#              min_score / min_margin only write the `flag` column as a hint for the manual
#              check. FALSE = old behaviour: a failed test gives "ambiguous", which is in
#              neither normal nor malig.
# malig_over_exclude :
#              a cluster that exclude wins but whose malig score is also high becomes malig
#              (oligo-like malignant). Glioma cells borrow oligodendrocyte-lineage programs
#              (the OPC-like / oligo-like states), so PLP1/MOG are high and SOX2/PTPRZ1 equally
#              high, and by argmax the cluster would be dropped as oligodendrocytes. The
#              z-scores of SOX2/PTPRZ1 in normal oligodendrocytes are negative and do not pass.
#              TRUE (default) = threshold min_score (scaled: 0; raw: min_score is -Inf, so the
#              mean of the malig column is used, same meaning); a number = that threshold;
#              FALSE = off. Only exclude gets this exit, not normal: otherwise APOE-high
#              myeloid clusters would be pulled into potential_malig.
classify_clusters <- function(dot_plot,
                              sets       = MARKER_SETS,
                              use        = c("scaled", "raw"),
                              top_n      = 3,
                              min_pct    = 10,
                              rank_by    = c("exp_pct", "value"),
                              agg_fun    = c("sum", "mean"),
                              min_score  = NULL,
                              min_margin = NULL,
                              force_call = TRUE,
                              malig_over_exclude = TRUE,
                              verbose    = TRUE) {
  use     <- match.arg(use)
  rank_by <- match.arg(rank_by)
  agg_fun <- match.arg(agg_fun)
  if (is.null(top_n)) top_n <- Inf
  if (is.null(min_margin))
    min_margin <- if (agg_fun == "sum") 0.5 * min(top_n, 3) else 0.5
  d <- if (inherits(dot_plot, "ggplot")) dot_plot$data else dot_plot
  need <- c("features.plot", "id", "avg.exp", "avg.exp.scaled", "pct.exp")
  if (!all(need %in% colnames(d)))
    stop("The $data of a DotPlot is needed, missing columns: ", paste(setdiff(need, colnames(d)), collapse = ", "))
  if (is.null(min_score)) min_score <- if (use == "scaled") 0 else -Inf

  d$features.plot <- as.character(d$features.plot)
  d$id            <- as.character(d$id)
  d$value         <- if (use == "scaled") d$avg.exp.scaled else log1p(d$avg.exp)
  d$rank_key      <- if (rank_by == "exp_pct") d$avg.exp * d$pct.exp / 100 else d$value
  if (use == "raw" && agg_fun == "sum" && !is.finite(top_n) && verbose)
    warning("use='raw' + agg_fun='sum' + top_n=Inf: log1p is always positive, so the group with the most genes almost always wins the sum. ",
            "With raw use agg_fun='mean' or keep top_n.", call. = FALSE)

  present <- unique(d$features.plot)
  # only the scored genes are checked; PLOT_ONLY genes may be absent
  miss    <- setdiff(unname(unlist(sets, use.names = FALSE)), present)
  if (length(miss) && verbose)
    message("[classify_clusters] markers not in the object, ignored: ", paste(miss, collapse = ", "))

  sets_ok <- lapply(sets, intersect, present)
  if (any(lengths(sets_ok) == 0))
    stop("No marker matched for these gene sets: ",
         paste(names(sets_ok)[lengths(sets_ok) == 0], collapse = ", "))

  clusters <- sort_clusters(d$id)
  f <- if (agg_fun == "sum") sum else mean

  # Per cluster x gene set: first drop the genes with pct.exp < min_pct (z-scores of almost
  # unexpressed genes are not reliable), then sort by rank_key, take the top_n and combine their
  # values. Nothing left = the marker set is not expressed in the cluster, which stays NA
  # (treated as -Inf below)
  score <- matrix(NA_real_, length(clusters), length(sets_ok),
                  dimnames = list(clusters, names(sets_ok)))
  pct   <- score
  picked <- matrix(NA_character_, nrow(score), ncol(score), dimnames = dimnames(score))
  n_used <- matrix(0L, nrow(score), ncol(score), dimnames = dimnames(score))
  vals   <- vector("list", nrow(score) * ncol(score))   # values of each cell, sorted by rank_key
  dim(vals) <- dim(score)
  for (ci in seq_along(clusters)) for (si in seq_along(sets_ok)) {
    i <- which(d$id == clusters[ci] & d$features.plot %in% sets_ok[[si]] &
                 d$pct.exp >= min_pct)
    i <- head(i[order(d$rank_key[i], decreasing = TRUE)], top_n)
    n_used[ci, si] <- length(i)
    vals[[ci, si]] <- d$value[i]
    if (!length(i)) next
    score[ci, si]  <- f(d$value[i], na.rm = TRUE)
    pct[ci, si]    <- mean(d$pct.exp[i], na.rm = TRUE)   # for reference only, not used in the call
    picked[ci, si] <- paste(d$features.plot[i], collapse = ",")
  }

  s <- score; s[is.na(s)] <- -Inf
  ord    <- t(apply(s, 1, order, decreasing = TRUE))
  best   <- s[cbind(seq_len(nrow(s)), ord[, 1])]
  second <- if (ncol(s) >= 2) s[cbind(seq_len(nrow(s)), ord[, 2])] else rep(-Inf, nrow(s))
  margin <- best - second

  add_flag <- function(flag, hit, tag)
    ifelse(hit, ifelse(nzchar(flag), paste0(flag, "+", tag), tag), flag)

  # `weak` is what marks a call as not reliable (with force_call = FALSE it turns the call into
  # ambiguous); `flag` is the label for humans and also records "a rule changed the call"
  # (oligo_like_malig)
  weak <- (best < min_score) | (margin < min_margin)
  flag <- add_flag(add_flag(rep("", nrow(s)), best < min_score, "low_score"),
                   margin < min_margin, "low_margin")

  # Under sum the groups must add the same number of terms to be comparable. After the min_pct
  # filter and with unequal set sizes (malig has only SOX2/PTPRZ1), a winner with one term more
  # than the runner-up gets a free point. Test: cut the winner to as many genes as the runner-up
  # has and recompute; if it then loses, flag uneven_n. Only clusters that win by having more
  # terms are hit, clear wins are not
  if (agg_fun == "sum" && ncol(s) >= 2) {
    unfair <- vapply(seq_len(nrow(s)), function(k) {
      v <- vals[[k, ord[k, 1]]]
      n2 <- n_used[k, ord[k, 2]]
      is.finite(best[k]) && length(v) > n2 && sum(head(v, n2)) < second[k]
    }, logical(1))
    flag <- add_flag(flag, unfair, "uneven_n")
  }

  call <- colnames(s)[ord[, 1]]      # the group with the highest score, "the most similar one"

  # Oligo-like malignant: see malig_over_exclude above. Only exclude gets this exit.
  if (!identical(malig_over_exclude, FALSE) && all(c("malig", "exclude") %in% colnames(s))) {
    thr <- if (isTRUE(malig_over_exclude)) {
      # scaled: min_score = 0 means "above the average of all clusters"; raw: min_score is -Inf,
      # which as a threshold would flip unconditionally, so fall back to the mean of the malig
      # column, which has the same meaning
      if (is.finite(min_score)) min_score else mean(s[is.finite(s[, "malig"]), "malig"])
    } else malig_over_exclude
    hit <- call == "exclude" & is.finite(s[, "malig"]) & s[, "malig"] >= thr
    call[hit] <- "malig"
    flag <- add_flag(flag, hit, "oligo_like_malig")
  }

  # A cluster without any score (all NA) cannot be assigned to anything
  call[!is.finite(best)] <- "ambiguous"
  if (!force_call) call[weak] <- "ambiguous"

  # top_genes / pct_top follow the final call, not the group with the highest score: a cluster
  # changed by oligo_like_malig must show the genes that carry its malig score, otherwise the
  # table would say call = malig with top_genes = PLP1,MOG,MAG
  win_col <- ord[, 1]
  m       <- match(call, colnames(score))          # "ambiguous" -> NA, keeps ord[, 1]
  win_col[!is.na(m)] <- m[!is.na(m)]
  win <- cbind(seq_len(nrow(score)), win_col)
  res <- data.frame(cluster = clusters,
                    round(score, 2),
                    call      = call,
                    margin    = round(margin, 2),
                    n_used    = apply(n_used, 1, paste, collapse = "/"),
                    flag      = flag,
                    pct_top   = round(pct[win], 1),
                    top_genes = picked[win],   # the genes that carry the score of the winning group
                    row.names = NULL, stringsAsFactors = FALSE)

  if (verbose) {
    cat(sprintf("\n== classify_clusters (use=%s, rank_by=%s, top_n=%s, min_pct=%s, agg_fun=%s, min_score=%s, min_margin=%s, force_call=%s, malig_over_exclude=%s) ==\n",
                use, rank_by, top_n, min_pct, agg_fun, min_score, min_margin,
                force_call, malig_over_exclude))
    cat(sprintf("   most genes per group (min(top_n, genes in the group)): %s\n",
                paste(names(sets_ok), paste0(pmin(lengths(sets_ok), top_n), "/", lengths(sets_ok)),
                      sep = "=", collapse = "  ")))
    cat(sprintf("   n_used column = genes actually scored per group (%s); genes with pct.exp < %s%% are filtered out\n",
                paste(names(sets_ok), collapse = "/"), min_pct))
    print(res, row.names = FALSE)
    ol <- clusters[grepl("oligo_like_malig", flag) & call == "malig"]
    if (length(ol))
      cat("\n>> oligo-like: exclude scored highest but malig is high too, changed to potential_malig: ",
          paste(ol, collapse = ", "), "\n", sep = "")
    fl <- clusters[weak & call != "ambiguous"]
    if (length(fl))
      cat("\n!! called by the most similar group, but score / margin are not reliable, check these first: ",
          paste(fl, collapse = ", "), "\n", sep = "")
    amb <- clusters[call == "ambiguous"]
    if (length(amb))
      cat("\n!! ambiguous (in neither normal nor potential_malig, must be decided by hand): ",
          paste(amb, collapse = ", "), "\n", sep = "")
    cat("\nCan be pasted back into a script:\n")
    fmt <- function(v) if (length(v)) paste0("c(", paste0("'", v, "'", collapse = ", "), ")") else "character(0)"
    cat("obvious_normal  <- ", fmt(clusters[call == "normal"]),  "\n", sep = "")
    cat("excluded        <- ", fmt(clusters[call == "exclude"]), "\n", sep = "")
    cat("potential_malig <- ", fmt(clusters[call == "malig"]),   "\n\n", sep = "")
  }

  list(obvious_normal  = clusters[call == "normal"],
       excluded        = clusters[call == "exclude"],
       potential_malig = clusters[call == "malig"],   # ambiguous clusters are not put in here
       ambiguous       = clusters[call == "ambiguous"],  # usually empty with force_call = TRUE
       flagged         = clusters[nzchar(flag)],      # called, but not reliable: check first
       table           = res,
       score           = score,
       n_used          = n_used,   # cluster x set: genes actually scored in each cell
       picked          = picked)   # cluster x set: the top_n genes chosen for scoring
}

# =============================================================================
# Samples
# =============================================================================
SAMPLE_META <- read.delim(SAMPLE_META_FILE, stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
SIDS <- SAMPLE_META$sid

.args <- commandArgs(trailingOnly = TRUE)
sids <- if (length(.args)) intersect(SIDS, .args) else SIDS
rds_path <- setNames(file.path(PER_SAMPLE_DIR, paste0(sids, "_minf400.rds")), sids)
if (any(!file.exists(rds_path))) stop("No rds for these samples: ", paste(sids[!file.exists(rds_path)], collapse = ", "))
sids <- sids[order(-file.size(rds_path[sids]))]   # large samples first
message(length(sids), " samples: ", paste(sids, collapse = ", "))

OVERRIDES <- read.csv(OVERRIDE_CSV, colClasses = "character", stringsAsFactors = FALSE)
stopifnot(all(c("sid", "cluster", "call") %in% names(OVERRIDES)),
          all(OVERRIDES$call %in% c("malig", "normal", "exclude")))
message(nrow(OVERRIDES), " manual overrides <- ", OVERRIDE_CSV)

# =============================================================================
# One sample
# =============================================================================
run_one <- function(sid) {
  png_out <- file.path(OUT_DIR, paste0(sid, "_filter_malignant.png"))
  tsv_out <- file.path(SAVE_TO,  paste0(sid, "_filtered_malignant.tsv"))

  # ---- 1. load + scale.data -------------------------------------------------
  seurat <- readRDS(rds_path[sid])
  DefaultAssay(seurat) <- "RNA"
  j <- try(SeuratObject::JoinLayers(seurat, assay = "RNA"), silent = TRUE)
  if (!inherits(j, "try-error")) seurat <- j
  stopifnot("seurat_clusters" %in% names(seurat@meta.data), "umap" %in% Reductions(seurat))
  n_cells <- ncol(seurat)
  seurat$seurat_clusters <- as.character(seurat$seurat_clusters)

  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = NUM_HVG, verbose = FALSE)
  seurat <- ScaleData(seurat, verbose = FALSE)
  m <- LayerData(seurat, assay = "RNA", layer = "scale.data")
  colnames(m) <- colnames(seurat)
  rownames(m) <- VariableFeatures(seurat)

  # ---- 2. centroid of every cluster -----------------------------------------
  cls <- as.character(unique(seurat$seurat_clusters))
  celltype.sig.mat <- sapply(cls, function(k) {
    cells <- intersect(colnames(m), rownames(seurat@meta.data)[seurat$seurat_clusters == k])
    rowMeans(m[, cells, drop = FALSE])
  })
  colnames(celltype.sig.mat) <- cls

  # ---- 3. cluster calls + manual overrides ----------------------------------
  gns <- intersect(marker_genes(), rownames(seurat))
  umap_dim <- DimPlot(seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
  dot_plot <- DotPlot(seurat, features = gns, group.by = "seurat_clusters") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  calls <- classify_clusters(dot_plot, min_margin = 0.1)
  obvious_normal  <- calls$obvious_normal
  excluded        <- calls$excluded
  potential_malig <- calls$potential_malig

  ov <- OVERRIDES[OVERRIDES$sid == sid, , drop = FALSE]
  if (nrow(ov)) {
    bad <- setdiff(ov$cluster, cls)
    if (length(bad)) stop(sid, ": these clusters of the override table do not exist: ", paste(bad, collapse = ", "))
    obvious_normal  <- union(setdiff(obvious_normal,  ov$cluster), ov$cluster[ov$call == "normal"])
    excluded        <- union(setdiff(excluded,        ov$cluster), ov$cluster[ov$call == "exclude"])
    potential_malig <- union(setdiff(potential_malig, ov$cluster), ov$cluster[ov$call == "malig"])
    message("!! manual override ", paste(sprintf("%s->%s", ov$cluster, ov$call), collapse = " "))
  }
  if (!length(potential_malig)) stop(sid, ": no potential_malig cluster")

  seurat$annotation <- annotate_clusters(seurat$seurat_clusters, obvious_normal, excluded, potential_malig)
  umap_dim2 <- DimPlot(seurat, reduction = "umap", group.by = "annotation", label = TRUE,
                       cols = ANNO_COLORS[levels(seurat$annotation)])

  # ---- 4. maximum correlation of every cell with the two sets of centroids ---
  if (!length(obvious_normal)) {
    warning(sid, ": no obvious_normal cluster, the normal side condition has no effect")
    seurat$maxscore_tonormal <- rep(-Inf, ncol(seurat))
  } else {
    seurat$maxscore_tonormal <- apply(cor(m, celltype.sig.mat[, obvious_normal, drop = FALSE], method = "pearson"), 1, max)
  }
  seurat$maxscore_tomalig <- apply(cor(m, celltype.sig.mat[, potential_malig, drop = FALSE], method = "pearson"), 1, max)

  # ---- 5. thresholds: Otsu on maxscore_tomalig of the first malig cluster -----
  k0 <- sort_clusters(potential_malig)[1]
  in_k0 <- seurat$seurat_clusters == k0
  v_malig <- seurat$maxscore_tomalig[in_k0]
  t_malig <- tryCatch({
    if (sum(in_k0) < 20 || diff(range(v_malig)) < 1e-8) stop("too few cells or degenerate values")
    otsu(as.matrix(v_malig), range = range(v_malig), levels = 1024)
  }, error = function(e) {
    warning(sid, ": Otsu failed (", conditionMessage(e), "), falling back to ", T_MALIG_FALLBACK)
    T_MALIG_FALLBACK
  })
  t_normal <- T_NORMAL
  message(sprintf("thresholds T_NORMAL=%.4f  T_MALIG=%.4f (otsu, cluster %s)", t_normal, t_malig, k0))

  hist_df <- rbind(data.frame(score = "maxscore_tonormal", value = seurat$maxscore_tonormal[in_k0]),
                   data.frame(score = "maxscore_tomalig",  value = v_malig))
  hist_df <- hist_df[is.finite(hist_df$value), ]
  hist_df$score <- factor(hist_df$score, levels = c("maxscore_tonormal", "maxscore_tomalig"))
  thr_df <- data.frame(score = factor(c("maxscore_tonormal", "maxscore_tomalig"), levels = levels(hist_df$score)),
                       thr = c(t_normal, t_malig))
  p_hist <- ggplot(hist_df, aes(value)) +
    geom_histogram(bins = 50, fill = "grey85", colour = "white", linewidth = 0.2) +
    geom_vline(data = thr_df, aes(xintercept = thr), linetype = "dashed", colour = "red") +
    facet_wrap(~ score, nrow = 1, scales = "free") +
    labs(title = sprintf("%s - cluster %s (otsu)", sid, k0), x = NULL, y = "cells") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())

  # ---- 6. two-level filter ------------------------------------------------------
  keep <- colnames(seurat)[seurat$maxscore_tonormal < t_normal &
                           seurat$maxscore_tomalig  > t_malig &
                           seurat$seurat_clusters %in% potential_malig]
  n_kept <- length(keep)
  message("kept ", n_kept, " / ", n_cells, sprintf("  (%.1f%%)", 100 * n_kept / n_cells))

  # ---- 7. plot ------------------------------------------------------------------
  seurat$filter_result <- factor(ifelse(colnames(seurat) %in% keep, "kept", "dropped"), levels = c("dropped", "kept"))
  dim_filter <- DimPlot(seurat, reduction = "umap", group.by = "filter_result", order = "kept",
                        cols = c(dropped = "#d9d9d9", kept = "#B2182B")) +
    ggtitle(sprintf("%s - filter_malignant: %d / %d kept", sid, n_kept, n_cells))
  ggsave(png_out, (umap_dim + umap_dim2 + dot_plot) / (dim_filter + p_hist), width = 15, height = 10)

  tb <- table(seurat$seurat_clusters, seurat$filter_result)
  keep_rate <- round(100 * prop.table(tb, 1)[, "kept"], 1)

  # ---- 8. export the cNMF input -------------------------------------------------
  exported <- FALSE
  if (n_kept < MIN_KEEP) {
    warning(sid, ": only ", n_kept, " cells left, no matrix exported")
  } else {
    cnt <- LayerData(seurat, assay = "RNA", layer = "counts")[, keep, drop = FALSE]
    write.table(as.matrix(Matrix::t(cnt)), file = tsv_out, quote = FALSE, sep = "\t")   # cells x genes
    rm(cnt); gc(verbose = FALSE)
    exported <- TRUE
  }

  ct <- calls$table
  ct$sid <- sid
  ct$final_call <- ifelse(ct$cluster %in% potential_malig, "malig",
                   ifelse(ct$cluster %in% obvious_normal,  "normal",
                   ifelse(ct$cluster %in% excluded,        "exclude", "ambiguous")))
  ct$overridden <- ct$call != ct$final_call
  ct <- ct[, c("sid", setdiff(names(ct), "sid"))]

  list(
    summary = data.frame(
      sid = sid, n_cells = n_cells, n_kept = n_kept, pct_kept = round(100 * n_kept / n_cells, 1),
      t_normal = t_normal, t_malig = round(t_malig, 4), k0 = k0, n_k0 = sum(in_k0),
      n_malig_cl = length(potential_malig), n_normal_cl = length(obvious_normal), n_exclude_cl = length(excluded),
      malig_clusters   = paste(sort_clusters(potential_malig), collapse = "|"),
      normal_clusters  = paste(sort_clusters(obvious_normal),  collapse = "|"),
      exclude_clusters = paste(sort_clusters(excluded),        collapse = "|"),
      overridden = paste(sprintf("%s->%s", ov$cluster, ov$call), collapse = "|"),
      flagged = paste(calls$flagged, collapse = "|"), exported = exported, stringsAsFactors = FALSE),
    calls = ct,
    keep_rate = data.frame(sid = sid, cluster = names(keep_rate), pct_kept = as.numeric(keep_rate),
                           n_cells = as.integer(table(seurat$seurat_clusters)[names(keep_rate)]),
                           row.names = NULL, stringsAsFactors = FALSE))
}

# =============================================================================
# Parallel run (PSOCK: every global that run_one uses has to be sent to the workers)
# =============================================================================
EXPORT_VARS <- c("rds_path", "OUT_DIR", "SAVE_TO", "T_NORMAL", "T_MALIG_FALLBACK",
                 "NUM_HVG", "MIN_KEEP", "OVERRIDES", "run_one",
                 "MARKER_SETS", "PLOT_ONLY", "ANNO_COLORS",
                 "marker_genes", "annotate_clusters", "sort_clusters", "classify_clusters")
cl <- makeCluster(min(NCORES, length(sids)))
registerDoParallel(cl)
clusterExport(cl, EXPORT_VARS, envir = environment())
t_all <- Sys.time()
res <- foreach(sid = sids, .errorhandling = "pass",
               .packages = c("Seurat", "ggplot2", "patchwork", "EBImage", "Matrix")) %dopar% {
  lg <- file(file.path(OUT_DIR, "logs", paste0(sid, ".log")), open = "wt")
  sink(lg); sink(lg, type = "message")
  on.exit({ sink(type = "message"); sink(); close(lg) })
  tryCatch(run_one(sid), error = function(e) { message("FAILED: ", conditionMessage(e)); e })
}
names(res) <- sids
stopCluster(cl)

bad <- sids[vapply(res, inherits, logical(1), "error")]
for (s in bad) message("!! ", s, " failed (logs/", s, ".log): ", conditionMessage(res[[s]]))
ok <- setdiff(sids, bad)
bind <- function(field) do.call(rbind, lapply(res[ok], `[[`, field))
summ <- bind("summary"); summ <- summ[order(summ$sid), ]
write.table(summ, file.path(OUT_DIR, "auto_filter_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(bind("calls"),     file.path(OUT_DIR, "cluster_calls.csv"), row.names = FALSE)
write.csv(bind("keep_rate"), file.path(OUT_DIR, "cluster_keep_rate.csv"), row.names = FALSE)
print(summ[, c("sid", "n_cells", "n_kept", "pct_kept", "t_malig", "k0", "overridden", "exported")], row.names = FALSE)
message("total time ", round(difftime(Sys.time(), t_all, units = "mins"), 1), " min")
if (length(bad)) stop("These samples failed: ", paste(bad, collapse = ", "))
