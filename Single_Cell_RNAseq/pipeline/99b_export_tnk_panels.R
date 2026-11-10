# =============================================================================
# Export the tables behind the T/NK panels
#
# Companion of 99_export_figure_data.R: the T/NK panels need things that script does not
# write (the T/NK embedding, PDCD1 expression in the T/NK cells, the dot plot summary and
# the reference signature scores), all of which come from the T/NK object. Nothing in
# ../figures/ opens an object; the figure scripts read what is written here.
#
# Input : scRNA_lymphocyte_final.rds (4,201 T/NK cells, 18 samples) from ../objects/released/,
#         or from the directory in the environment variable SCRNA_RELEASE if it is set
#         ../sample_metadata.tsv
#         ../data/cell_level/cell_annotation.tsv.gz  (written by 99_export_figure_data.R; malignant
#         cells per sample)
#         ../inputs/refsets/zhang2021_table_s3.xlsx  (Table S3 of Zheng et al., Science 2021,
#         PMID 34914499; the file name carries the name of the senior author)
# Output: ../data/panel_data/
#           tnk_samples.tsv             one row per sample (19): patient_id, T/NK cells in the object,
#                                       malignant cells; patient_id holds SP IDs
#           tnk_cells.tsv.gz            one row per T/NK cell: state, T/NK UMAP, PDCD1, and
#                                       whether the cell is in the pooled composition bars
#           tnk_marker_dotplot.tsv      dot plot summary: state x marker gene
#           tnk_signature_scores.tsv    mean module score of each reference signature per state
#           tnk_signature_genes.tsv     the genes of each reference signature, as scored
#
# Paths are relative to this directory (pipeline/): run the script from here.
# Run time about 1 minute.
# =============================================================================

suppressPackageStartupMessages(library(Seurat))

RELEASE <- Sys.getenv("SCRNA_RELEASE", "../objects/released")
HERE    <- ".."
OUT     <- file.path(HERE, "data", "panel_data")
META    <- read.delim(file.path(HERE, "sample_metadata.tsv"), colClasses = "character")
REF_XLS <- file.path(HERE, "inputs", "refsets", "zhang2021_table_s3.xlsx")

STATE_LEVELS <- c("GZMK+ CD8 T", "CD8 Tex", "CD8 T_m", "CD4 T_m", "CD56-CD16+ NK", "Stress T",
                  "CD56+CD16- NK", "Prolif", "Treg", "Naive T", "IFN")

# marker genes of the dot plot, by the state they mark (the order of the columns)
MARKERS <- list("CD8" = c("CD8A", "CD8B"), "GZMK+ CD8" = "GZMK", "Tex" = c("PDCD1", "HAVCR2", "LAG3"),
                "CD4" = c("CD40LG", "CD4", "IL7R"), "NK16" = c("FGFBP2", "KLRF1", "FCGR3A"),
                "Stress" = c("HSPA1A", "HSPA1B"), "NK56" = c("KRT81", "TMIGD2", "NCAM1"),
                "Prolif" = c("MKI67", "TYMS", "PCLAF"), "Treg" = c("FOXP3", "IL2RA", "RTKN2"),
                "Naive" = c("CCR7", "LEF1", "SELL"), "IFN" = c("ISG15", "MX1", "IFIT1"))

# reference signatures: top SCORE_TOP_N genes of each meta-cluster after removing genes that
# separate patients or clonotypes; signatures with fewer than SCORE_MIN_GENES genes are dropped
SCORE_TOP_N     <- 30L
SCORE_MIN_GENES <- 5L

say  <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n"); flush(stdout()) }
wtsv <- function(df, path, gz = FALSE) {
  f <- file.path(OUT, path)
  if (gz) { con <- gzfile(f, "wb"); write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con) }
  else      write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
  say(sprintf("  %-34s %7d rows  %6.2f MB", path, nrow(df), file.size(f) / 1024^2))
}

is_uninformative <- function(g) {
  grepl("^MT-|^TR[ABDG][VJ][0-9]|^IG[HKL][VDJ][0-9]", g) |
    (grepl("^RP[LS]P?[0-9]+[A-Z]{0,2}[0-9]?$", g) & !grepl("^RPS6K", g)) |
    g %in% c("XIST", "TSIX", "RPS4Y1", "DDX3Y", "UTY", "EIF1AY", "KDM5D", "USP9Y", "ZFY", "TXLNGY", "NLGN4Y", "PRKY")
}

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say("reading scRNA_lymphocyte_final.rds")
o <- readRDS(file.path(RELEASE, "scRNA_lymphocyte_final.rds"))
m <- o@meta.data
say("  cells", ncol(o), " samples", length(unique(m$orig.ident)))
stopifnot(all(as.character(m$tnk_state) %in% STATE_LEVELS), !anyNA(m$tnk_state))
state <- factor(as.character(m$tnk_state), levels = STATE_LEVELS)
sid   <- as.character(m$orig.ident)

# sample -> patient. Patients are the unit of the comparisons: TL_6 to TL_10 are five sites of
# one tumour. Column patient_id holds the SP ID of Table S1, "S" + p0x_id (P38 -> SP38), never a
# hospital number; the mapping was checked against Table S1 (cohort, subtype, tumour location and the
# scRNA-seq flag agree for all 15 patients).
j <- match(sid, META$sid)
if (anyNA(j)) stop("sid not in sample_metadata.tsv: ", paste(unique(sid[is.na(j)]), collapse = ", "))
stopifnot(grepl("^P[0-9]{2}$", META$p0x_id))
patient <- paste0("S", META$p0x_id)[j]

## ---- 0. samples ---------------------------------------------------------------------
# All 19 samples, including PL_6_4 (renamed PL_9): it has 16 T/NK cells, below the 30-cell minimum
# of the T/NK analysis, so it is absent from the T/NK object (n_tnk = 0) but still has malignant cells.
ann <- read.delim(gzfile(file.path(HERE, "data", "cell_level", "cell_annotation.tsv.gz")),
                  colClasses = c(cell_id = "character", sid = "character"))
stopifnot(is.logical(ann$is_malignant))
samples <- data.frame(sid = META$sid, renamed_id = META$renamed_id, group = META$group,
                      patient_id = paste0("S", META$p0x_id),
                      n_tnk = as.vector(table(factor(sid, levels = META$sid))),
                      n_malignant = as.vector(table(factor(ann$sid[ann$is_malignant], levels = META$sid))))
stopifnot(sum(samples$n_tnk) == ncol(o), sum(samples$n_malignant) == sum(ann$is_malignant))
wtsv(samples, "tnk_samples.tsv")

## ---- 1. cells: state, T/NK embedding, PDCD1 ------------------------------------------
emb <- Embeddings(o, "umap")[, 1:2]
stopifnot(identical(rownames(emb), rownames(m)))

# PDCD1 as the analysis used it: log-normalised RNA assay (not the integrated assay)
pdcd1 <- as.vector(LayerData(o, layer = "data", assay = "RNA")["PDCD1", ])
pdcd1_r <- round(pdcd1, 6)
stopifnot(identical(pdcd1 > 0, pdcd1_r > 0), identical(pdcd1 > 1, pdcd1_r > 1))   # rounding changes no call

# Pooled composition bars: the five samples of the multi-sample patient would otherwise make that
# patient 919 of the 4,201 cells, so each of their samples contributes a random 20% of its cells
# (ceiling), drawn with seed 1; every other sample contributes all of its cells. The draw is
# stored per cell so that the bars can be redrawn without repeating it.
set.seed(1)
multi <- names(which(table(patient[!duplicated(sid)]) > 1))     # patients with more than one sample
keep  <- unlist(lapply(split(seq_along(sid), sid), function(i)
  if (patient[i[1]] %in% multi) i[sample.int(length(i), ceiling(0.2 * length(i)))] else i))
in_bar <- seq_along(sid) %in% keep
say("  pooled bars: multi-sample patient(s)", paste(multi, collapse = ","), "| cells", length(sid), "->", length(keep))

cells <- data.frame(cell_id = rownames(m), sid = sid, renamed_id = META$renamed_id[j], group = META$group[j],
                    patient_id = patient, tnk_state = as.character(state),
                    umap_1 = round(emb[, 1], 4), umap_2 = round(emb[, 2], 4),
                    in_pooled_bar = in_bar, PDCD1_lognorm = pdcd1_r)
wtsv(cells, "tnk_cells.tsv.gz", gz = TRUE)

## ---- 2. dot plot summary -------------------------------------------------------------
# DotPlot() on the RNA assay (the default assay of the integrated object is the CCA-corrected
# one): avg.exp is the mean of the de-logged expression, pct.exp the percentage of cells with a
# non-zero value, avg.exp.scaled the z-score of log1p(avg.exp) across the 11 states, clipped to +-2.5.
genes <- unname(unlist(MARKERS))
stopifnot(!anyDuplicated(genes), all(genes %in% rownames(o[["RNA"]])))
dp <- suppressWarnings(DotPlot(o, features = genes, assay = "RNA", group.by = "tnk_state"))$data
dp <- data.frame(tnk_state = as.character(dp$id), gene = as.character(dp$features.plot),
                 gene_group = rep(names(MARKERS), lengths(MARKERS))[match(as.character(dp$features.plot), genes)],
                 n_cells = as.integer(table(state)[as.character(dp$id)]),
                 avg_exp = round(dp$avg.exp, 6), pct_exp = round(dp$pct.exp, 6),
                 avg_exp_scaled = round(dp$avg.exp.scaled, 6))
dp <- dp[order(match(dp$tnk_state, STATE_LEVELS), match(dp$gene, genes)), ]
stopifnot(nrow(dp) == length(STATE_LEVELS) * length(genes))
wtsv(dp, "tnk_marker_dotplot.tsv")

## ---- 3. reference signature scores ---------------------------------------------------
# Meta-cluster signatures of Zheng et al. 2021 (Table S3, sheets CD8 and CD4; genes ranked by
# comb.ES.rank), scored per cell with AddModuleScore on the RNA assay (50 control genes per
# signature gene, seed 1), then averaged within each state.
ref <- do.call(rbind, lapply(c("CD8", "CD4"), function(s) {
  d <- openxlsx::read.xlsx(REF_XLS, sheet = s, startRow = 2)
  d <- data.frame(gene = d[[1]], signature = d[[2]], rank = as.numeric(d[[6]]))
  d[!is.na(d$gene) & !is.na(d$signature) & !is_uninformative(d$gene), ]
}))
ref  <- ref[order(ref$signature, ref$rank), ]
sets <- lapply(split(ref[, c("gene", "rank")], ref$signature), function(d) {
  d <- d[seq_len(min(SCORE_TOP_N, nrow(d))), ]
  d[d$gene %in% rownames(o[["RNA"]]), ]
})
dropped <- names(sets)[vapply(sets, nrow, 1L) < SCORE_MIN_GENES]
say("  signatures", length(sets), "| dropped for < ", SCORE_MIN_GENES, "genes:", paste(dropped, collapse = ", "))
sets <- sets[vapply(sets, nrow, 1L) >= SCORE_MIN_GENES]

o  <- suppressWarnings(AddModuleScore(o, features = unname(lapply(sets, `[[`, "gene")), name = "zh_",
                                      assay = "RNA", ctrl = 50, seed = 1))
sc <- as.matrix(o@meta.data[, paste0("zh_", seq_along(sets))]); colnames(sc) <- names(sets)
n_state <- as.vector(table(state))                       # in STATE_LEVELS order
mean_sc <- rowsum(sc, as.character(state))[STATE_LEVELS, ] / n_state
stopifnot(all(rownames(mean_sc) == STATE_LEVELS))
scores <- data.frame(tnk_state = rep(STATE_LEVELS, times = ncol(mean_sc)),
                     signature = rep(colnames(mean_sc), each = nrow(mean_sc)),
                     n_cells = rep(n_state, times = ncol(mean_sc)),
                     mean_score = round(as.vector(mean_sc), 6))
wtsv(scores, "tnk_signature_scores.tsv")

sig_genes <- do.call(rbind, lapply(names(sets), function(s)
  data.frame(signature = s, gene_order = seq_len(nrow(sets[[s]])), gene = sets[[s]]$gene,
             comb_ES_rank = sets[[s]]$rank)))
wtsv(sig_genes, "tnk_signature_genes.tsv")

say("done")
