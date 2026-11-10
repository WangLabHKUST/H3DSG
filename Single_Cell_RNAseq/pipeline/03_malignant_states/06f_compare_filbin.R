# =============================================================================
# Metaprograms against the Filbin reference states (second reference)
#
# Input : ../../inputs/refsets/filbin_gene_set.xlsx   8 gene sets
#         ../../results/malignant_states/04_metaprograms/plateau_v1_all19/cellular_states_metaprograms.txt
# Output: ../../results/malignant_states/04_metaprograms/plateau_v1_all19/
#             metaprogram_vs_reference_filbin_overlap.csv
#             metaprogram_vs_reference_filbin.pdf, metaprogram_vs_reference_filbin.png
#
# OPC-like-1/2/3 of the Filbin file are merged (union of the genes) into one OPC-like, only so
# that the reference has the same number of states as there are metaprograms. The three share
# only 1-3 genes pairwise (OPC-like-2 is in fact a ribosomal program and OPC-like-3 an
# immediate-early-gene stress program), so this is a rough comparison only.
# Rows are ordered by the column of their maximum, which puts the best matches on the diagonal.
#
# Runtime a few seconds. R 4.5.1.
# Paths are relative to this directory: run the script from here.
# =============================================================================

REF_FILBIN_XLSX <- "../../inputs/refsets/filbin_gene_set.xlsx"
RES_DIR    <- "../../results/malignant_states"
DIR_04_MP  <- file.path(RES_DIR, "04_metaprograms")
MP_VERSION <- "plateau_v1_all19"
DIR_MP     <- file.path(DIR_04_MP, MP_VERSION)

suppressPackageStartupMessages({ library(pheatmap); library(readxl) })

OPC_COLS  <- c("OPC-like-1", "OPC-like-2", "OPC-like-3")
COL_ORDER <- c("S", "G2M", "AC-like", "MES-like", "OPC-like", "OC-like")
.genes <- function(x) unique(trimws(as.character(x[!is.na(x) & trimws(x) != ""])))

raw <- as.data.frame(read_excel(REF_FILBIN_XLSX, sheet = 1))
ref <- lapply(setdiff(names(raw), OPC_COLS), function(n) .genes(raw[[n]]))
names(ref) <- setdiff(names(raw), OPC_COLS)
ref[["OPC-like"]] <- .genes(unlist(raw[, OPC_COLS], use.names = FALSE))
ref <- ref[c(intersect(COL_ORDER, names(ref)), setdiff(names(ref), COL_ORDER))]

mp <- read.delim(file.path(DIR_MP, "cellular_states_metaprograms.txt"), check.names = FALSE, colClasses = "character")
ovl <- as.matrix(sapply(ref, function(r) sapply(mp, function(m) length(intersect(.genes(m), r)))))
ovl <- ovl[order(apply(ovl, 1, which.max), -apply(ovl, 1, max)), , drop = FALSE]   # rows ordered by the column of their maximum -> diagonal
write.csv(ovl, file.path(DIR_MP, "metaprogram_vs_reference_filbin_overlap.csv"))
print(ovl)

draw <- function() {
  pheatmap(ovl, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = FALSE, border_color = "white",
           cellwidth = 46, cellheight = 34, angle_col = 45,
           main = sprintf("%s vs Filbin reference states\n(shared genes; OPC-like = union of OPC-like-1/2/3)", MP_VERSION))
}
pdf(file.path(DIR_MP, "metaprogram_vs_reference_filbin.pdf"), width = 2 + 0.9 * ncol(ovl), height = 2 + 0.6 * nrow(ovl)); draw(); dev.off()
png(file.path(DIR_MP, "metaprogram_vs_reference_filbin.png"), width = 200 + 90 * ncol(ovl), height = 200 + 60 * nrow(ovl),
    res = 130, bg = "white"); draw(); dev.off()
message("-> ", file.path(DIR_MP, "metaprogram_vs_reference_filbin.{pdf,png}"))
