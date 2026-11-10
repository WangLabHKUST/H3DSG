# =============================================================================
# hox_dotplots / 02_compute_stats.R -- gene set statistics behind the dot plots
#
# For each group of cells and each cell type group (immune / stromal / malignant):
#   avg_exp = mean over the genes of the set and all cells (zeros included) of the
#             log-normalised expression (RNA "data" layer)
#   pct_exp = percentage of the cells in which at least one gene of the set is detected
#   n_cells = number of cells
#
# Tables written to ../data/hox_dotplots/ (read by the plotting scripts):
#   location_dotplot_by_group.tsv   4 gene sets x 3 cell types in 6 groups: public thalamic and pontine
#                                   DMG cells and this cohort split by location (Cervical / Thoracic)
#                                   and subtype (PL / TL)
#   location_dotplot_by_sample.tsv  4 gene sets x 3 cell types in each of the 19 samples
#   pax3_nkx61_dotplot.tsv          PAX3 and NKX6-1 x 3 cell types in the pontine and thalamic DMG
#                                   cells and in the PL and TL samples
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages(library(Seurat))
source("hox_helpers.R")

sample_tab    <- read_sample_table()
gene_sets_tab <- read.delim(file.path(DATA_DIR, "location_gene_sets.tsv"), stringsAsFactors = FALSE)
gene_sets     <- df_to_genesets(gene_sets_tab, gene_col = "gene", set_loc_col = "gene_set")

analysis_obj <- readRDS(file.path(RESULTS_DIR, "hox_subset.rds"))
stopifnot(all(unlist(gene_sets) %in% rownames(analysis_obj)),
          setequal(levels(analysis_obj$celltype2), CELLTYPES), !anyNA(analysis_obj$celltype2),
          setequal(as.character(analysis_obj$orig.ident), sample_tab$sid))
print(table(analysis_obj$celltype2))

# public brain DMG cells; the spinal cord cells of that data set are not this cohort
brain <- readRDS(BRAIN_LIST_RDS)
brain <- brain[!names(brain) %in% "spinal"]          # keys: thalamic, pontine
stopifnot(all(c("thalamic", "pontine") %in% names(brain)))

stats_by <- function(obj_lst, sets, ...) calculate_multi_samp_gene_set_stats_advanced(
  obj_lst, sets, group.by = "celltype2", pct_method = PCT_METHOD, ...)
cell_counts <- function(obj_lst) {     # one row per object and cell type
  do.call(rbind, lapply(names(obj_lst), function(nm) {
    n <- table(factor(as.character(obj_lst[[nm]]$celltype2), levels = CELLTYPES))
    data.frame(orig.ident = nm, id = names(n), n_cells = as.integer(n), stringsAsFactors = FALSE)
  }))
}
write_tsv <- function(x, file) write.table(x, file.path(DATA_DIR, file), sep = "\t", quote = FALSE,
                                           row.names = FALSE, na = "")

# ---- by group: Fig. location dot plot ------------------------------------------------
# location x subtype groups of this cohort; panel_location carries the rule of 00_config.R
# (TL_4 is counted as Cervical)
group_sids <- list()
for (loc in c("Cervical", "Thoracic")) for (sub in c("PL", "TL")) {
  sids <- sample_tab$sid[sample_tab$panel_location == loc & sample_tab$subtype == sub]
  stopifnot(length(sids) > 0)
  group_sids[[paste(loc, sub, sep = "-")]] <- sids
}
stopifnot(!anyDuplicated(unlist(group_sids)), setequal(unlist(group_sids), sample_tab$sid))   # all 19 samples, once
obj_by_group <- lapply(group_sids, function(s) subset(analysis_obj, cells = colnames(analysis_obj)[analysis_obj$orig.ident %in% s]))
for (nm in names(group_sids))
  message(nm, ": ", length(group_sids[[nm]]), " samples / ", ncol(obj_by_group[[nm]]), " cells (",
          paste(sample_tab$renamed_id[match(group_sids[[nm]], sample_tab$sid)], collapse = ", "), ")")
obj_by_group[["Thalamus"]]  <- brain[["thalamic"]]
obj_by_group[["Brainstem"]] <- brain[["pontine"]]
group_order <- c("Thalamus", "Brainstem", "Cervical-PL", "Cervical-TL", "Thoracic-PL", "Thoracic-TL")

by_group <- merge(stats_by(obj_by_group, gene_sets), cell_counts(obj_by_group), by = c("orig.ident", "id"))
by_group <- data.frame(group = by_group$orig.ident, celltype = by_group$id, gene_set = by_group$features.plot,
                       n_cells = by_group$n_cells, avg_exp = by_group$avg.exp, pct_exp = by_group$pct.exp,
                       stringsAsFactors = FALSE)
by_group$source  <- ifelse(by_group$group %in% c("Thalamus", "Brainstem"), "public brain DMG data", "this study")
by_group$samples <- vapply(by_group$group, function(g) {
  if (g %in% names(group_sids)) paste(sample_tab$renamed_id[match(group_sids[[g]], sample_tab$sid)], collapse = ", ") else ""
}, "")
by_group <- by_group[order(match(by_group$group, group_order), match(by_group$celltype, CELLTYPES),
                           match(by_group$gene_set, GENE_SET_LEVELS)), , drop = FALSE]
stopifnot(nrow(by_group) == 6 * 3 * 4, !anyNA(by_group[, c("n_cells", "avg_exp", "pct_exp")]))
write_tsv(by_group, "location_dotplot_by_group.tsv")

# ---- by sample: each of the 19 samples ----------------------------------------------------
obj_by_sample <- setNames(lapply(sample_tab$sid, function(s) subset(analysis_obj, subset = orig.ident == s)), sample_tab$sid)
by_sample <- merge(stats_by(obj_by_sample, gene_sets), cell_counts(obj_by_sample), by = c("orig.ident", "id"))
by_sample <- data.frame(sid = by_sample$orig.ident, celltype = by_sample$id, gene_set = by_sample$features.plot,
                        n_cells = by_sample$n_cells, avg_exp = by_sample$avg.exp, pct_exp = by_sample$pct.exp,
                        stringsAsFactors = FALSE)
k <- match(by_sample$sid, sample_tab$sid)
by_sample <- cbind(sample_tab[k, c("sid", "renamed_id", "patient_id", "subtype", "tumor_location", "panel_order")],
                   by_sample[, c("celltype", "gene_set", "n_cells", "avg_exp", "pct_exp")])
by_sample <- by_sample[order(by_sample$panel_order, match(by_sample$celltype, CELLTYPES),
                             match(by_sample$gene_set, GENE_SET_LEVELS)), , drop = FALSE]
stopifnot(nrow(by_sample) == 19 * 3 * 4, !anyNA(by_sample),
          sum(by_sample$n_cells[by_sample$gene_set == "Thoracic"]) == ncol(analysis_obj))
write_tsv(by_sample, "location_dotplot_by_sample.tsv")
tl4 <- by_sample[by_sample$sid == "TL_4" & by_sample$celltype == "stromal" & by_sample$gene_set == "Thoracic", ]
message("TL_4, stromal cells, Thoracic gene set: ", round(tl4$pct_exp / 100 * tl4$n_cells), " of ", tl4$n_cells, " cells express it")

# ---- PAX3 and NKX6-1 ---------------------------------------------------------------------------
pax_sets <- list("PAX3" = "PAX3", "NKX6-1" = "NKX6-1")
pax_obj  <- readRDS(file.path(RESULTS_DIR, "pax3_nkx61_objects.rds"))
pax_order <- c("pontine", "thalamic", "PL", "TL")
stopifnot(setequal(names(pax_obj), pax_order))
pax <- merge(stats_by(pax_obj, pax_sets, pos.only = FALSE), cell_counts(pax_obj), by = c("orig.ident", "id"))
pax <- data.frame(group = pax$orig.ident, celltype = pax$id, gene = pax$features.plot,
                  n_cells = pax$n_cells, avg_exp = pax$avg.exp, pct_exp = pax$pct.exp, stringsAsFactors = FALSE)
pax <- pax[order(match(pax$group, pax_order), match(pax$celltype, CELLTYPES), match(pax$gene, names(pax_sets))), , drop = FALSE]
stopifnot(nrow(pax) == 4 * 3 * 2, !anyNA(pax))
write_tsv(pax, "pax3_nkx61_dotplot.tsv")
message("02_compute_stats done")
