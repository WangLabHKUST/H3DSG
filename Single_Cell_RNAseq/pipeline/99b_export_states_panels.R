# =============================================================================
# Export the tables behind the malignant-state, metaprogram, pseudotime and HOX panels
#
# 99_export_figure_data.R writes the tables shared by many panels. This script writes the
# ones that come from the malignant-state analysis (step 03), the TL_6 trajectory (step 08)
# and the HOX / PAX3 / NKX6-1 dot plots. Nothing in ../figures/ opens an object: every panel
# reads one of the tables written here.
#
# Input : ../results/malignant_states/04_metaprograms/plateau_v1_all19/   (step 03, script 06)
#             program_similarity_heatmap_data.rds, cellular_states_metaprograms.txt
#         ../inputs/annotations/tumor_mp_state_names_plateau_v1_all19.csv
#         ../inputs/refsets/filbin_gene_set.xlsx    reference states of Liu et al. 2022
#         ../inputs/hox/loc_gene.csv                gene -> anatomical level
#         ../results/trajectory_TL6/tables/tl6_malignant_cells.csv                  (step 08)
#         ../objects/released/scRNA_final.rds       (sections 5 and 6 only)
#         ../sample_metadata.tsv
# Output: ../data/metaprograms/   program similarity matrix, program -> metaprogram table,
#                                 signature genes, overlap with the Liu et al. states
#         ../data/panel_data/     heat map matrix of 5,000 malignant cells, per-cell
#                                 pseudotime, HOX / PAX3 / NKX6-1 dot plot statistics,
#                                 resection order of the samples
#
# Paths are relative to this directory (pipeline/): run the script from here. Three
# environment variables point at the inputs when they are not where the pipeline puts them:
#   STATES_DIR    directory with the step 03 metaprogram files (default as above)
#   TRAJ_CELLS    per-cell table written by step 08
#   SCRNA_FINAL   the released all-cell object
#   SKIP_OBJECT=1 skip sections 5 and 6, which read scRNA_final.rds
#
# Runtime about 5 minutes, almost all of it reading scRNA_final.rds. Seurat 5.2.1, R 4.5.1.
# =============================================================================

suppressPackageStartupMessages({ library(Matrix); library(readxl) })

STATES_DIR  <- Sys.getenv("STATES_DIR", "../results/malignant_states/04_metaprograms/plateau_v1_all19")
TRAJ_CELLS  <- Sys.getenv("TRAJ_CELLS", "../results/trajectory_TL6/tables/tl6_malignant_cells.csv")
SCRNA_FINAL <- Sys.getenv("SCRNA_FINAL", "../objects/released/scRNA_final.rds")
STATE_NAMES <- "../inputs/annotations/tumor_mp_state_names_plateau_v1_all19.csv"
LIU_REF     <- "../inputs/refsets/filbin_gene_set.xlsx"
LOC_GENE    <- "../inputs/hox/loc_gene.csv"
META_FILE   <- "../sample_metadata.tsv"
OUT_MP      <- "../data/metaprograms"
OUT_PANEL   <- "../data/panel_data"
for (d in c(OUT_MP, OUT_PANEL)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

STATE_ORDER <- c("S", "G2M", "OPC-like", "AC-like", "MES-like", "OC-like")
HEATMAP_N   <- 5000      # malignant cells drawn in the state heat map
HEATMAP_SEED <- 123      # seed of that random draw
# HOX genes, anterior to posterior, in the order the multi-site dot plot lists them
HOX_DISPLAY_ORDER <- c("HOXA1", "HOXB1", "HOXD1", "HOXA2", "HOXB2", "HOXA3", "HOXB3", "HOXD3",
                       "HOXA4", "HOXB4", "HOXC4", "HOXD4", "HOXA5", "HOXB5", "HOXC5", "HOXA6",
                       "HOXB6", "HOXC6", "HOXA7", "HOXB7", "HOXB8", "HOXC8", "HOXD8", "HOXA9",
                       "HOXB9", "HOXC9", "HOXD9", "HOXA10", "HOXC10", "HOXD10", "HOXA11",
                       "HOXC11", "HOXD11", "HOXC12", "HOXD12", "HOXA13", "HOXB13", "HOXC13",
                       "HOXD13")
THALAMUS_GENES <- c("ZIC3", "ZIC1", "ZIC4", "OTX1")   # replace the thalamic set of loc_gene.csv

say  <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n"); flush(stdout()) }
wtsv <- function(df, f, gz = FALSE) {
  if (gz) { con <- gzfile(f, "wb"); write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = ""); close(con) }
  else write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  say(sprintf("  %-62s %7d rows %4d cols %7.2f MB", sub("^\\.\\./data/", "", f), nrow(df), ncol(df), file.size(f) / 1024^2))
}

META <- read.delim(META_FILE, colClasses = "character")
# patient-level grouping uses the SP ID: "S" + p0x_id (P38 -> SP38)
META$sp_id <- paste0("S", META$p0x_id)
state_map <- read.csv(STATE_NAMES, stringsAsFactors = FALSE)           # metaprogram M1..M6 -> state
state_of  <- setNames(state_map$state, state_map$metaprogram)

# =============================================================================
# 1. Metaprograms: similarity matrix, assignment, signature genes, overlap with Liu et al.
# =============================================================================
say("1. metaprogram tables")
hm  <- readRDS(file.path(STATES_DIR, "program_similarity_heatmap_data.rds"))
sim <- as.matrix(hm$cnmf_intersect)
ann <- hm$ann[colnames(sim), ]
stopifnot(identical(rownames(sim), colnames(sim)), nrow(sim) == 107,
          identical(rownames(ann), colnames(sim)), all(sim == t(sim)))
extraction_order <- seq_len(ncol(sim))                    # order of the programs as extracted (sample by sample)
clustered <- hm$hm_order                                  # order on the hierarchical clustering tree
block     <- hm$ordered_p                                 # order in the metaprogram-block heat map (outliers dropped)
stopifnot(setequal(clustered, rownames(sim)), all(block %in% rownames(sim)), length(block) == 84)

pname <- setNames(ann$Program_renamed, ann$Program)       # program name with the renamed sample id
sim_out <- sim[clustered, clustered]
dimnames(sim_out) <- list(pname[clustered], pname[clustered])
wtsv(data.frame(program = rownames(sim_out), sim_out, check.names = FALSE),
     file.path(OUT_MP, "program_similarity_matrix.tsv"))

asg <- data.frame(
  program           = pname[rownames(ann)],
  sample            = ann$Sample_renamed,
  metaprogram       = ann$Metaprogram,
  state             = unname(state_of[ann$Metaprogram]),
  extraction_order  = extraction_order,
  clustered_order   = match(rownames(ann), clustered),
  block_order       = match(rownames(ann), block),
  cor_S_score       = round(ann$S.Score, 6),
  cor_G2M_score     = round(ann$G2M.Score, 6),
  row.names = NULL)
asg <- asg[order(asg$clustered_order), ]
wtsv(asg, file.path(OUT_MP, "program_metaprogram_assignment.tsv"))

mpg <- read.delim(file.path(STATES_DIR, "cellular_states_metaprograms.txt"), check.names = FALSE,
                  colClasses = "character", na.strings = "")
stopifnot(identical(colnames(mpg), paste0("M", 1:6)), nrow(mpg) == 50)
sig <- do.call(rbind, lapply(colnames(mpg), function(m)
  data.frame(metaprogram = m, state = unname(state_of[m]), rank = seq_len(nrow(mpg)), gene = mpg[[m]])))
wtsv(sig, file.path(OUT_MP, "metaprogram_signature_genes.tsv"))

# overlap with the reference states of Liu et al. 2022: size of the intersection of the two
# gene sets (50 signature genes against the 50 genes of the reference state); the three
# OPC-like programs of the reference are merged into one set
genes_of <- function(x) unique(trimws(as.character(x[!is.na(x) & trimws(x) != ""])))
raw <- as.data.frame(read_excel(LIU_REF, sheet = 1))
opc_cols <- c("OPC-like-1", "OPC-like-2", "OPC-like-3")
ref <- lapply(setdiff(names(raw), opc_cols), function(n) genes_of(raw[[n]]))
names(ref) <- setdiff(names(raw), opc_cols)
ref[["OPC-like"]] <- genes_of(unlist(raw[, opc_cols], use.names = FALSE))
ref_cols <- c("S", "G2M", "AC-like", "MES-like", "OPC-like", "OC-like")
ovl <- sapply(ref[ref_cols], function(r) sapply(mpg, function(m) length(intersect(genes_of(m), r))))
rows <- order(apply(ovl, 1, which.max), -apply(ovl, 1, max))     # rows ordered so that the diagonal is high
ovl <- ovl[rows, , drop = FALSE]
ovl_long <- data.frame(
  metaprogram = rep(rownames(ovl), times = ncol(ovl)),
  state       = rep(unname(state_of[rownames(ovl)]), times = ncol(ovl)),
  liu_state   = rep(colnames(ovl), each = nrow(ovl)),
  shared_genes = as.vector(ovl),
  row_order   = rep(seq_len(nrow(ovl)), times = ncol(ovl)),
  col_order   = rep(seq_len(ncol(ovl)), each = nrow(ovl)))
wtsv(ovl_long, file.path(OUT_MP, "metaprogram_overlap_liu2022.tsv"))

# =============================================================================
# 2. Sample table: resection order (cervical to thoracic) and site number
# =============================================================================
say("2. sample table")
# Tumour location -> spinal level interval: M = 0, C1..C8 = 1..8, T1..T12 = 9..20, L1..L5 = 21..25;
# Dorsal / Ventral / edge suffixes are dropped
loc_span <- function(loc) {
  n <- vapply(toupper(strsplit(as.character(loc), "-")[[1]]), function(t) {
    if (identical(t, "M")) return(0)
    m <- regmatches(t, regexec("^([CTLS])([0-9]+)$", t))[[1]]
    if (!length(m)) return(NA_real_)
    unname(c(C = 0, T = 8, L = 20, S = 25)[m[2]] + as.numeric(m[3]))
  }, numeric(1))
  n <- n[!is.na(n)]
  stopifnot(length(n) > 0)
  c(min(n), max(n))
}
sp <- vapply(META$tumor_location, loc_span, numeric(2))
META$resection_order <- match(META$sid, META$sid[order(sp[2, ], sp[1, ], META$sid)])   # caudal end, cranial end, sid
META$site_number <- ave(as.numeric(META$location_sort), META$sp_id, FUN = function(x) rank(x))
samples <- data.frame(renamed_id = META$renamed_id, patient_id = META$sp_id, group = META$group,
                      tumor_location = META$tumor_location, resection_order = META$resection_order,
                      site_number = as.integer(META$site_number))
samples <- samples[order(samples$resection_order), ]
wtsv(samples, file.path(OUT_PANEL, "sample_resection_order.tsv"))

# =============================================================================
# 3. Per-cell pseudotime of the TL_6 malignant cells (step 08)
# =============================================================================
say("3. pseudotime")
pt <- read.csv(TRAJ_CELLS, check.names = FALSE, stringsAsFactors = FALSE)
j <- match(pt$orig.ident, META$sid)
stopifnot(!anyNA(j), nrow(pt) == 8491)
pt_out <- data.frame(cell_id = pt$cell_id, sid = pt$orig.ident, renamed_id = META$renamed_id[j],
                     patient_id = META$sp_id[j], tumor_location = META$tumor_location[j],
                     site = pt$loc_name, cluster_res0.1 = pt[["cluster_res0.1"]],
                     Lineage1 = round(pt$Lineage1, 5), Lineage2 = round(pt$Lineage2, 5),
                     Lineage3 = round(pt$Lineage3, 5), pt_mean = round(pt$pt_mean, 5))
wtsv(pt_out, file.path(OUT_PANEL, "pseudotime_TL6_malignant_cells.tsv"))

if (Sys.getenv("SKIP_OBJECT") == "1") { say("SKIP_OBJECT=1, done"); quit(save = "no") }

# =============================================================================
# 4. Read the released object once
# =============================================================================
suppressPackageStartupMessages(library(Seurat))
say("reading", SCRNA_FINAL)
obj <- readRDS(SCRNA_FINAL)
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = "RNA", layers = c("counts", "data"), dimreducs = NULL, graphs = NULL)
md <- obj@meta.data
stopifnot(ncol(obj) == 132343, !any(md$integrated_snn_res.1_anno == "myeloid_mixed_T"))
sid_of <- as.character(md$orig.ident)
renamed_of <- META$renamed_id[match(sid_of, META$sid)]
stopifnot(!anyNA(renamed_of))
say("  cells", ncol(obj), " genes", nrow(obj))

# =============================================================================
# 5. Heat map of the malignant states: 5,000 random malignant cells x 297 signature genes
#    (the draw made by the plotting code, seed 123); values as drawn by DoHeatmap
# =============================================================================
say("5. state heat map")
mods <- mpg; colnames(mods) <- unname(state_of[colnames(mods)])
features <- lapply(mods, function(col) { v <- col[!is.na(col) & nzchar(trimws(col))]; unique(v)[1:50] })
features_present <- lapply(features, function(g) intersect(g, rownames(obj)))[STATE_ORDER]
gns <- unique(unlist(features_present))
# the scaling is done on all 132,343 cells, then the malignant cells are taken
obj <- ScaleData(obj, features = gns, verbose = FALSE)
tumor_cells <- colnames(obj)[obj@meta.data$integrated_snn_res.1_anno == "tumor"]
tumor_obj <- subset(obj, cells = tumor_cells)
tumor_obj$State2025 <- factor(as.character(tumor_obj$State2025), levels = STATE_ORDER)
stopifnot(ncol(tumor_obj) == 81304, !anyNA(tumor_obj$State2025), length(gns) == 297)

set.seed(HEATMAP_SEED)
heat_cells <- sample(colnames(tumor_obj), size = HEATMAP_N, replace = FALSE)
state_colors <- setNames(c("#d53e4f", "#d53e4f", "#f46d43", "#fdae61", "#fcc5c0", "#fee08b"), STATE_ORDER)
dhp <- DoHeatmap(tumor_obj, features = gns, group.by = "State2025", cells = heat_cells,
                 angle = 30, size = 3, group.colors = state_colors)
d <- dhp$data                                             # Cell, Feature, Expression, Identity as drawn
stopifnot(all(c("Cell", "Feature", "Expression", "Identity") %in% colnames(d)))
cell_levels <- levels(d$Cell); feat_levels <- levels(d$Feature)
mat <- matrix(NA_real_, length(cell_levels), length(feat_levels), dimnames = list(cell_levels, feat_levels))
mat[cbind(as.integer(d$Cell), as.integer(d$Feature))] <- d$Expression
# DoHeatmap inserts empty spacer columns between the state blocks (all values NA); they are not cells
is_cell <- rowSums(!is.na(mat)) > 0
stopifnot(sum(is_cell) == HEATMAP_N, setequal(cell_levels[is_cell], heat_cells),
          length(feat_levels) == length(gns), setequal(feat_levels, gns),
          !anyNA(mat[is_cell, ]))
# state of each cell from the object (the Identity column of DoHeatmap's data is not reliable: Seurat
# indexes it with a factor); the cell levels are already ordered by state
ident <- as.character(tumor_obj$State2025)[match(cell_levels, colnames(tumor_obj))]
# draw_position: column of the cell in the drawn heat map, spacers included.
# Gene columns are in the top-to-bottom order of the heat map.
heat <- data.frame(cell_id = cell_levels[is_cell], draw_position = which(is_cell),
                   sid = as.character(tumor_obj$orig.ident)[match(cell_levels[is_cell], colnames(tumor_obj))],
                   state = ident[is_cell], round(mat[is_cell, rev(feat_levels), drop = FALSE], 3), check.names = FALSE)
heat$renamed_id <- META$renamed_id[match(heat$sid, META$sid)]
stopifnot(!is.unsorted(match(heat$state, STATE_ORDER)), identical(as.vector(table(heat$state)[STATE_ORDER]),
          as.vector(table(factor(tumor_obj$State2025[match(heat$cell_id, colnames(tumor_obj))], STATE_ORDER)))))
# cell key = sample (renamed ID) + 16-nt barcode; the cell names of the object carry sequencing-run names
bc <- sub("^.*_([ACGT]{16})(-1)?$", "\\1", heat$cell_id)
stopifnot(all(nchar(bc) == 16))
heat$cell <- paste(heat$renamed_id, bc, sep = "_")
stopifnot(!anyDuplicated(heat$cell))
heat <- heat[, c("cell", "draw_position", "sid", "renamed_id", "state", rev(feat_levels))]
rownames(heat) <- NULL
say("  drawn columns incl. spacers: ", length(cell_levels), "; seed ", HEATMAP_SEED)
wtsv(heat, file.path(OUT_PANEL, "malignant_state_heatmap_5000cells.tsv.gz"), gz = TRUE)
say("  heat map: ", nrow(heat), " cells x ", length(gns), " genes; states: ",
    paste(sprintf("%s=%d", names(table(heat$state)), table(heat$state)), collapse = " "))
rm(tumor_obj, dhp, d, mat); invisible(gc())

# =============================================================================
# 6. Dot plot statistics: HOX genes and PAX3 / NKX6-1 by sample and cell type
#    cell type: tumor = malignant; B, myeloid, neutrophil, T/NK = immune; oligo, vascular = stromal
# =============================================================================
say("6. dot plot statistics")
anno <- as.character(obj@meta.data$integrated_snn_res.1_anno)
ct2 <- ifelse(anno == "tumor", "malignant",
       ifelse(anno %in% c("B", "myeloid", "neutrophil", "T/NK"), "immune",
       ifelse(anno %in% c("oligo", "vascular"), "stromal", NA)))
stopifnot(!anyNA(ct2), identical(as.vector(table(factor(ct2, c("immune", "stromal", "malignant")))),
                                  c(48591L, 2448L, 81304L)))
X <- LayerData(obj[["RNA"]], layer = "data")

loc <- read.csv(LOC_GENE, stringsAsFactors = FALSE, check.names = FALSE)
loc <- loc[loc$Location != "Thalamus", ]
loc <- rbind(loc, data.frame(gene_name = THALAMUS_GENES, Location = "Thalamus"))
loc <- loc[!duplicated(loc$gene_name), ]
hox_genes <- intersect(unique(c(loc$gene_name, HOX_DISPLAY_ORDER)), rownames(X))
gene_tab <- data.frame(gene = hox_genes,
                       location_group = loc$Location[match(hox_genes, loc$gene_name)],
                       order_all_samples = match(hox_genes, loc$gene_name),
                       order_multisite = match(hox_genes, HOX_DISPLAY_ORDER))
wtsv(gene_tab, file.path(OUT_PANEL, "hox_gene_annotation.tsv"))

# Per sample and cell type, for each gene: mean log-normalised expression (mean_log_expr, the
# value the dot plots colour by), mean of the de-logged expression (mean_delog_expr, which is
# Seurat::DotPlot's avg.exp) and the percentage of cells with a value > 0 (pct_exp)
dot_stats <- function(genes) {
  out <- list()
  for (s in META$sid) for (ct in c("immune", "stromal", "malignant")) {
    idx <- which(sid_of == s & ct2 == ct)
    if (!length(idx)) next
    x <- X[genes, idx, drop = FALSE]
    out[[length(out) + 1]] <- data.frame(
      renamed_id = META$renamed_id[META$sid == s], celltype = ct, gene = genes, n_cells = length(idx),
      mean_log_expr = as.numeric(Matrix::rowMeans(x)),
      mean_delog_expr = as.numeric(Matrix::rowMeans(expm1(x))),
      pct_exp = as.numeric(Matrix::rowMeans(x > 0)) * 100)
  }
  do.call(rbind, out)
}
wtsv(dot_stats(hox_genes), file.path(OUT_PANEL, "hox_dotplot_stats_by_sample_celltype.tsv"))
wtsv(dot_stats(c("PAX3", "NKX6-1")), file.path(OUT_PANEL, "pax3_nkx61_dotplot_stats_by_sample_celltype.tsv"))
say("done")
