# =============================================================================
# Export every table the figure scripts read
#
# This is the only place where a released Seurat object is turned into the tables in
# ../data/. Nothing in ../figures/ opens an object; if a panel needs something that
# ../data/ does not hold, add it here.
#
# Input : scRNA_final.rds, scRNA_lymphocyte_final.rds, scRNA_myeloid_final.rds
#         from ../objects/released/, or from the directory in the environment
#         variable SCRNA_RELEASE if it is set
#         ../sample_metadata.tsv
# Output: ../data/cell_level/   per-cell annotation, UMAP coordinates, state scores
#         ../data/sample_level/ per-sample statistics and compositions
#         ../data/panel_data/   expression values for the panels that need them
#
# Paths are relative to this directory (pipeline/): run the script from here.
#
# Run time about 10 minutes, almost all of it reading scRNA_final.rds.
# =============================================================================

suppressPackageStartupMessages(library(Seurat))

RELEASE <- Sys.getenv("SCRNA_RELEASE", "../objects/released")
HERE    <- ".."
OUT     <- file.path(HERE, "data")
META    <- read.delim(file.path(HERE, "sample_metadata.tsv"), colClasses = "character")

# Genes whose expression a figure panel needs. A panel that needs a new gene is served by
# adding it here, never by loading an object in a figure script.
PANEL_GENES <- c("CXCL14", "SOX2", "PTPRZ1", "OLIG1", "OLIG2", "PDGFRA", "AQP4", "GFAP",
                 "MBP", "PLP1", "MKI67", "TOP2A", "PTPRC", "CD3D", "CD8A", "CD4", "GZMK",
                 "PDCD1", "P2RY12", "C1QC", "CD163", "VEGFA", "HIF1A")

for (d in c("cell_level", "sample_level", "metaprograms", "panel_data"))
  dir.create(file.path(OUT, d), showWarnings = FALSE, recursive = TRUE)

say  <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n"); flush(stdout()) }
wtsv <- function(df, path, gz = FALSE) {
  f <- file.path(OUT, path)
  if (gz) { con <- gzfile(f, "wb"); write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con) }
  else      write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
  say(sprintf("  %-46s %7d rows  %6.2f MB", path, nrow(df), file.size(f) / 1024^2))
}

# Table of the columns wanted in the per-cell annotation: new name <- name in the object.
# A column that is absent is reported and left out rather than silently dropped.
WANT <- c(sid = "orig.ident", nCount_RNA = "nCount_RNA", nFeature_RNA = "nFeature_RNA",
          percent.mito = "percent.mito", cell_type = "integrated_snn_res.1_anno",
          malignant_state = "State2025", tnk_state = "tnk_state",
          myeloid_state_coarse = "myeloid_state_coarse", myeloid_state_fine = "myeloid_state_fine",
          included_in_prop = "included_in_prop",
          substate_excluded_reason = "substate_excluded_reason",
          depth_Alt = "depth_Alt", depth_Total = "depth_Total",
          Alt_freq = "Alt_freq", is_alt = "is_alt")

say("reading scRNA_final.rds")
o <- readRDS(file.path(RELEASE, "scRNA_final.rds"))
m <- o@meta.data
say("  cells", ncol(o), " meta columns", ncol(m), " reductions:", paste(Reductions(o), collapse = ", "))

missing <- setdiff(WANT, colnames(m))
if (length(missing)) say("  !! absent from meta.data:", paste(missing, collapse = ", "))
WANT <- WANT[WANT %in% colnames(m)]

## ---- 1. per-cell annotation ----------------------------------------------------------
ann <- data.frame(cell_id = rownames(m), m[, WANT, drop = FALSE], check.names = FALSE)
colnames(ann)[-1] <- names(WANT)
j <- match(ann$sid, META$sid)
if (anyNA(j)) stop("sid not in sample_metadata.tsv: ", paste(unique(ann$sid[is.na(j)]), collapse = ", "))
ann$renamed_id <- META$renamed_id[j]
ann$group      <- META$group[j]
ann$patient_id <- META$patient_id[j]
for (col in intersect(c("percent.mito", "Alt_freq"), colnames(ann)))
  ann[[col]] <- round(as.numeric(ann[[col]]), 4)

say("  cell_type levels:", paste(sort(unique(as.character(ann$cell_type))), collapse = ", "))
ann$is_malignant <- as.character(ann$cell_type) %in% c("tumor", "malignant")
say(sprintf("  malignant %d   included_in_prop & malignant %d", sum(ann$is_malignant),
            if ("included_in_prop" %in% colnames(ann))
              sum(ann$is_malignant & ann$included_in_prop %in% c(TRUE, "TRUE")) else NA))

ann <- ann[, c("cell_id", "sid", "renamed_id", "group", "patient_id",
               setdiff(colnames(ann), c("cell_id", "sid", "renamed_id", "group", "patient_id")))]
wtsv(ann, "cell_level/cell_annotation.tsv.gz", gz = TRUE)

## ---- 2. UMAP coordinates -------------------------------------------------------------
# The released object carries umap (= 30 PCs, n.neighbors 150, the annotation embedding)
# and umap_pc20_correlation. In the unslimmed object the first is named umap_pc30 and a
# second, differently parameterised umap also exists, so prefer umap_pc30 when present.
main <- if ("umap_pc30" %in% Reductions(o)) "umap_pc30" else "umap"
say("  main UMAP:", main)
e1 <- Embeddings(o, main)[, 1:2]
e2 <- Embeddings(o, "umap_pc20_correlation")[, 1:2]
stopifnot(identical(rownames(e1), rownames(m)), identical(rownames(e2), rownames(m)))
wtsv(data.frame(cell_id = rownames(e1),
                umap_1 = round(e1[, 1], 4), umap_2 = round(e1[, 2], 4),
                umap_corr_1 = round(e2[, 1], 4), umap_corr_2 = round(e2[, 2], 4)),
     "cell_level/umap_coordinates.tsv.gz", gz = TRUE)

## ---- 3. malignant state scores -------------------------------------------------------
SC <- c("G2M", "S", "OPC_like", "AC_like", "MES_like", "OC_like")
if (all(SC %in% colnames(m))) {
  k <- ann$is_malignant
  sc <- data.frame(cell_id = rownames(m)[k], round(m[k, SC], 4),
                   state = m[k, "State2025"], state_score = round(m[k, "state.score"], 4))
  wtsv(sc, "cell_level/malignant_state_scores.tsv.gz", gz = TRUE)
} else say("  !! state score columns absent:", paste(setdiff(SC, colnames(m)), collapse = ", "))

## ---- 4. per-sample QC summary --------------------------------------------------------
# Pre-QC cell counts are not in the released object; they are in the logs of
# 01_preprocessing. What is exported here is the cohort as released.
qc <- do.call(rbind, lapply(split(ann, ann$sid), function(d) data.frame(
  sid = d$sid[1], renamed_id = d$renamed_id[1], group = d$group[1], patient_id = d$patient_id[1],
  n_cells = nrow(d),
  median_nCount = median(d$nCount_RNA), median_nFeature = median(d$nFeature_RNA),
  median_percent_mito = round(median(d$percent.mito), 3))))
qc <- qc[order(qc$group, qc$renamed_id), ]
wtsv(qc, "sample_level/qc_summary.tsv")

## ---- 5. compositions -----------------------------------------------------------------
# One row per sample and category, counts and the within-sample fraction. A comparison
# plot of a fraction between PL and TL reads one of these and has 19 points.
composition <- function(d, col, extra = NULL) {
  d <- d[!is.na(d[[col]]) & d[[col]] != "", ]
  if (!nrow(d)) return(NULL)
  out <- do.call(rbind, lapply(split(d, d$sid), function(s) {
    tb <- table(factor(as.character(s[[col]])))
    data.frame(sid = s$sid[1], renamed_id = s$renamed_id[1], group = s$group[1],
               patient_id = s$patient_id[1], category = names(tb), n = as.integer(tb),
               fraction = round(as.integer(tb) / nrow(s), 6), row.names = NULL)
  }))
  if (!is.null(extra)) out <- cbind(out, extra[match(out$sid, extra$sid), -1, drop = FALSE])
  out[order(out$group, out$renamed_id, out$category), ]
}
wtsv(composition(ann, "cell_type"), "sample_level/cell_type_composition.tsv")

mal <- ann[ann$is_malignant, ]
if ("included_in_prop" %in% colnames(ann)) {
  inprop <- do.call(rbind, lapply(split(mal, mal$sid), function(s)
    data.frame(sid = s$sid[1], included_in_prop = any(s$included_in_prop %in% c(TRUE, "TRUE")))))
  wtsv(composition(mal, "malignant_state", inprop), "sample_level/malignant_state_composition.tsv")
} else wtsv(composition(mal, "malignant_state"), "sample_level/malignant_state_composition.tsv")

wtsv(composition(ann, "tnk_state"),             "sample_level/tnk_state_composition.tsv")
wtsv(composition(ann, "myeloid_state_coarse"),  "sample_level/myeloid_state_composition.tsv")
wtsv(composition(ann, "myeloid_state_fine"),    "sample_level/myeloid_state_fine_composition.tsv")

## ---- 6. expression of the panel genes ------------------------------------------------
# Per cell, long format, for the genes listed in PANEL_GENES only. This is what lets a
# dot plot or a violin be drawn without the expression matrix: 23 genes over 132,343
# cells is a few MB, the matrix is gigabytes.
have <- intersect(PANEL_GENES, rownames(o[["RNA"]]))
if (length(setdiff(PANEL_GENES, have)))
  say("  !! genes absent from RNA:", paste(setdiff(PANEL_GENES, have), collapse = ", "))
ex <- as.matrix(LayerData(o[["RNA"]], layer = "data")[have, , drop = FALSE])
long <- data.frame(cell_id = rep(colnames(ex), each = length(have)),
                   gene    = rep(have, times = ncol(ex)),
                   expr    = round(as.vector(ex), 4))
wtsv(long[long$expr > 0, ], "panel_data/panel_gene_expression_nonzero.tsv.gz", gz = TRUE)

rm(o, ex, long); invisible(gc())

## ---- 7. the two subset objects: nothing new, only a check ----------------------------
# The subset objects hold no annotation that scRNA_final.rds does not; they are read only
# to confirm that, so that the figure layer can rely on the single per-cell table.
for (sp in list(c("scRNA_lymphocyte_final.rds", "tnk_state"),
                c("scRNA_myeloid_final.rds",    "myeloid_state_coarse"))) {
  s <- readRDS(file.path(RELEASE, sp[1]))
  sm <- s@meta.data
  shared <- intersect(rownames(sm), ann$cell_id)
  agree  <- identical(as.character(sm[shared, sp[2]]),
                      as.character(ann[[sp[2]]][match(shared, ann$cell_id)]))
  say(sprintf("  %-30s %6d cells, %6d in cell_annotation, %s agrees: %s",
              sp[1], ncol(s), length(shared), sp[2], agree))
  rm(s, sm); invisible(gc())
}

say("done")
