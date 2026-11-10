# =============================================================================
# Metaprograms against the K27scDSG reference states (the reference that names them)
#
# Input : ../../inputs/refsets/K27scDSG.cellular.states.6sates.202502.xlsx
#             reference gene sets, 6 states x 50 genes (version 2025-02; the earlier 2024-11
#             version still contained ribosomal and unnamed-clone genes, which lowers the
#             overlap systematically)
#         ../../results/malignant_states/04_metaprograms/plateau_v1_all19/cellular_states_metaprograms.txt
# Output: ../../results/malignant_states/04_metaprograms/plateau_v1_all19/
#             metaprogram_vs_reference_overlap.csv
#             metaprogram_vs_reference.pdf, metaprogram_vs_reference.png
#         ../../results/malignant_states/04_metaprograms/reference_states.tsv   (tsv copy of the reference)
#
# The overlap of two 50-gene sets is the number of genes they share (the same measure as the
# program similarity of 06_metaprograms.R). The best match of each metaprogram gives its state
# name, recorded in inputs/annotations/tumor_mp_state_names_plateau_v1_all19.csv, which
# 07_score_states_and_plot.R checks against this table. The plot title keeps the original
# "K27scDMG" label of this reference.
#
# Runtime a few seconds. R 4.5.1.
# Paths are relative to this directory: run the script from here.
# =============================================================================

REF_K27_XLSX <- "../../inputs/refsets/K27scDSG.cellular.states.6sates.202502.xlsx"
RES_DIR    <- "../../results/malignant_states"
DIR_04_MP  <- file.path(RES_DIR, "04_metaprograms")
MP_VERSION <- "plateau_v1_all19"
DIR_MP     <- file.path(DIR_04_MP, MP_VERSION)

suppressPackageStartupMessages({ library(pheatmap); library(readxl) })
if (!file.exists(REF_K27_XLSX)) stop("Reference gene set not found: ", REF_K27_XLSX)

ref <- as.data.frame(read_excel(REF_K27_XLSX, sheet = 1))
ref[] <- lapply(ref, function(x) trimws(as.character(x)))
message("reference: ", ncol(ref), " states x ", nrow(ref), " genes -> ", paste(names(ref), collapse = ", "))
write.table(ref, file.path(DIR_04_MP, "reference_states.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

mp <- read.delim(file.path(DIR_MP, "cellular_states_metaprograms.txt"), check.names = FALSE, colClasses = "character")
mp[] <- lapply(mp, trimws)
ovl <- as.matrix(sapply(ref, function(r) sapply(mp, function(m) length(intersect(na.omit(m[m != ""]), na.omit(r[r != ""]))))))
write.csv(ovl, file.path(DIR_MP, "metaprogram_vs_reference_overlap.csv"))
print(ovl)

draw <- function() {
  pheatmap(ovl, cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.0f", fontsize_number = 11,
           number_color = "black", border_color = "white",
           color = colorRampPalette(c("white", "#FFE9B0", "#F98E52", "#B2182B"))(100),
           breaks = seq(0, max(5, max(ovl)), length.out = 101),
           cellwidth = 46, cellheight = 34, angle_col = 45,
           main = sprintf("%s vs K27scDMG reference states\n(shared genes / 50)", MP_VERSION))
}
pdf(file.path(DIR_MP, "metaprogram_vs_reference.pdf"), width = 2 + 0.9 * ncol(ovl), height = 2 + 0.6 * nrow(ovl)); draw(); dev.off()
png(file.path(DIR_MP, "metaprogram_vs_reference.png"), width = 200 + 90 * ncol(ovl), height = 200 + 60 * nrow(ovl), res = 130); draw(); dev.off()
message("-> ", file.path(DIR_MP, "metaprogram_vs_reference.{pdf,png}"))
