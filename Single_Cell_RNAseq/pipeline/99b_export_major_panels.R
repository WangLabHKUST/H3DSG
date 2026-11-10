# =============================================================================
# Export the tables that the major cell type panels need and 99_export_figure_data.R
# does not write
#
# Why a second script: 99_export_figure_data.R writes the expression of the genes in its
# PANEL_GENES, and only five of the 17 marker genes of the cell type marker dot plot are in
# that list. Rather than reopen 99_export_figure_data.R, the missing genes are exported here,
# in the same long format.
#
# Input : scRNA_final.rds from ../objects/released/, or from the directory in the
#         environment variable SCRNA_RELEASE if it is set (SCRNA_OBJECT names another file in
#         it, e.g. scRNA_final_diet_rna.rds, which holds the same RNA data layer)
#         ../data/panel_data/panel_gene_expression_nonzero.tsv.gz  (only to check the overlap)
#         ../sample_metadata.tsv
# Output: ../data/panel_data/marker_gene_expression_nonzero.tsv.gz
#           cell_key, gene, expr   log-normalised RNA values of the 17 marker genes of the
#           cell type marker dot plot, non-zero entries only, same rounding as
#           panel_gene_expression_nonzero.tsv.gz
#
# cell_key is renamed_id + "_" + the 16-bp barcode (PL_1_AAACCCAAGAAGTATC), unique over the
# cohort. It carries no sequencing name, unlike the cell names of the object; to join it to
# cell_level/cell_annotation.tsv.gz build the same key there from renamed_id and the barcode
# in cell_id (the 16 letters A, C, G, T).
#
# The five genes that both files hold (PTPRZ1, SOX2, PTPRC, CD3D, PLP1) must agree exactly;
# the script stops if they do not.
#
# Paths are relative to this directory (pipeline/): run the script from here.
# Run time about 3 minutes (reading the object).
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

RELEASE <- Sys.getenv("SCRNA_RELEASE", "../objects/released")
OBJECT  <- Sys.getenv("SCRNA_OBJECT", "scRNA_final.rds")
OUT     <- "../data/panel_data"
META    <- read.delim("../sample_metadata.tsv", colClasses = "character")

# the dot plot of the cell type markers, top to bottom
MARKER_GENES <- c("PTPRZ1", "SOX2", "PTPRC", "AIF1", "CD14", "FCGR3B", "S100A9", "CSF3R", "CD3D",
                  "IL32", "NKG7", "CD79A", "MS4A1", "COL1A1", "VWF", "PLP1", "MAG")

say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n"); flush(stdout()) }

say("reading", OBJECT)
o <- readRDS(file.path(RELEASE, OBJECT))
stopifnot(all(MARKER_GENES %in% rownames(o[["RNA"]])))
x <- LayerData(o[["RNA"]], layer = "data")[MARKER_GENES, , drop = FALSE]
cells <- colnames(x)
sid <- as.character(o@meta.data[cells, "orig.ident"])
stopifnot(all(sid %in% META$sid))
cell_key <- paste0(META$renamed_id[match(sid, META$sid)], "_", sub("^.*_([ACGT]{16})(-1)?$", "\\1", cells))
stopifnot(!anyDuplicated(cell_key), all(grepl("^[A-Za-z0-9_]+_[ACGT]{16}$", cell_key)))
rm(o); invisible(gc())

# long format, non-zero after rounding to 4 decimals, ordered by cell and then gene
trip <- as.data.frame(Matrix::summary(as(x, "TsparseMatrix")))
trip$expr <- round(trip$x, 4)
trip <- trip[trip$expr > 0, ]
trip <- trip[order(trip$j, trip$i), ]
long <- data.frame(cell_key = cell_key[trip$j], gene = MARKER_GENES[trip$i], expr = trip$expr)
say("non-zero entries:", nrow(long), " cells:", length(cells))

# the genes that panel_gene_expression_nonzero.tsv.gz already holds must be identical
old_file <- file.path(OUT, "panel_gene_expression_nonzero.tsv.gz")
if (file.exists(old_file)) {
  old <- read.delim(gzfile(old_file), colClasses = c("character", "character", "numeric"))
  shared <- intersect(unique(old$gene), MARKER_GENES)
  old$cell_key <- cell_key[match(old$cell_id, cells)]          # the old file is keyed by the object's cell names
  if (anyNA(old$cell_key)) stop("cells of panel_gene_expression_nonzero.tsv.gz not in the object")
  a <- long[long$gene %in% shared, ]; b <- old[old$gene %in% shared, ]
  key <- function(d) d[order(d$gene, d$cell_key), ]
  a <- key(a); b <- key(b)
  ok <- nrow(a) == nrow(b) && identical(a$cell_key, b$cell_key) && identical(a$gene, b$gene) && all(a$expr == b$expr)
  say("genes in both files:", paste(shared, collapse = ", "), "- identical:", ok)
  if (!ok) stop("marker genes differ from panel_gene_expression_nonzero.tsv.gz")
}

stopifnot(!any(grepl("DHY|ZDK", long$cell_key)))
f <- file.path(OUT, "marker_gene_expression_nonzero.tsv.gz")
con <- gzfile(f, "wb"); write.table(long, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)
say(sprintf("%s  %d rows  %.2f MB", basename(f), nrow(long), file.size(f) / 1024^2))
say("done")
