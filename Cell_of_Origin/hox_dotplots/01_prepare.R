# =============================================================================
# hox_dotplots / 01_prepare.R -- read the integrated object once and keep what the statistics need
#
#   1. integrated object (19 samples, 132,343 cells): cell type group celltype2
#      (immune / stromal / malignant) from the annotation integrated_snn_res.1_anno
#   2. gene sets of the four anatomical levels (loc_gene.csv; the Thalamus set is replaced by
#      ZIC3, ZIC1, ZIC4 and OTX1) -> ../data/hox_dotplots/location_gene_sets.tsv
#   3. the object restricted to the genes of the four sets -> results/hox_subset.rds
#   4. PAX3 and NKX6-1 only: the public thalamic and pontine DMG cells (all genes of
#      location_objects_full.rds) and the H3-DSG-PL / H3-DSG-TL cells of this cohort
#      -> results/pax3_nkx61_objects.rds
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages(library(Seurat))
source("hox_helpers.R")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

sample_tab <- read_sample_table()

# ---- 1. integrated object + cell type groups -----------------------------------------
stopifnot(file.exists(INTEGRATED_RDS))
message("loading ", INTEGRATED_RDS)
merged <- readRDS(INTEGRATED_RDS)
DefaultAssay(merged) <- "RNA"
for (a in setdiff(Assays(merged), "RNA")) merged[[a]] <- NULL   # only the RNA data layer is read
stopifnot(ncol(merged) == 132343,
          setequal(unique(as.character(merged$orig.ident)), sample_tab$sid),
          !any(merged$integrated_snn_res.1_anno == "myeloid_mixed_T"))
merged$orig.ident <- factor(as.character(merged$orig.ident), levels = sample_tab$sid)

anno <- as.character(merged$integrated_snn_res.1_anno)
stopifnot(all(anno %in% names(CELLTYPE_MAP)))
merged$celltype2 <- factor(unname(CELLTYPE_MAP[anno]), levels = CELLTYPES)
print(table(merged$celltype2, useNA = "ifany"))
stopifnot(!anyNA(merged$celltype2),
          identical(as.vector(table(merged$celltype2)), c(48591L, 2448L, 81304L)))

# ---- 2. gene sets ----------------------------------------------------------------------
gene_loc <- read.csv(LOC_GENE_CSV, stringsAsFactors = FALSE, check.names = FALSE)
thalamus_markers <- c("ZIC3", "ZIC1", "ZIC4", "OTX1")
stopifnot(all(thalamus_markers %in% rownames(merged)))
gene_sets_tab <- gene_loc[gene_loc$Location != "Thalamus", , drop = FALSE]
gene_sets_tab <- rbind(gene_sets_tab, data.frame(gene_name = thalamus_markers, Location = "Thalamus",
                                                 stringsAsFactors = FALSE))
gene_sets_tab <- gene_sets_tab[!duplicated(gene_sets_tab$gene_name), , drop = FALSE]
stopifnot(setequal(gene_sets_tab$Location, GENE_SET_LEVELS), all(gene_sets_tab$gene_name %in% rownames(merged)))
write.table(data.frame(gene = gene_sets_tab$gene_name, gene_set = gene_sets_tab$Location),
            file.path(DATA_DIR, "location_gene_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("gene sets: ", paste(sprintf("%s = %d", names(table(gene_sets_tab$Location)), table(gene_sets_tab$Location)),
                             collapse = ", "))

# ---- 3. the genes of the four sets ------------------------------------------------------
saveRDS(merged[gene_sets_tab$gene_name, ], file.path(RESULTS_DIR, "hox_subset.rds"))

# ---- 4. PAX3 / NKX6-1 ---------------------------------------------------------------------
pax_genes <- c("PAX3", "NKX6-1")
stopifnot(all(pax_genes %in% rownames(merged)))
merged_pax <- merged[pax_genes, ]
rm(merged); invisible(gc())

stopifnot(file.exists(BRAIN_FULL_RDS))
brain <- readRDS(BRAIN_FULL_RDS)
brain <- brain[!names(brain) %in% "spinal"]          # the spinal cord cells of that data set are not this cohort
## the all-gene version has no celltype2: take it from celltype (Immune / Stromal / Malignant)
obj_list <- lapply(brain, function(o) {
  if (!"celltype2" %in% colnames(o@meta.data)) o$celltype2 <- tolower(as.character(o$celltype))
  stopifnot(all(o$celltype2 %in% CELLTYPES))
  subset(o, features = intersect(pax_genes, rownames(o)))
})
rm(brain); invisible(gc())
for (g in c("PL", "TL")) {
  cells_g <- colnames(merged_pax)[merged_pax$orig.ident %in% sample_tab$sid[sample_tab$subtype == g]]
  obj_list[[g]] <- subset(merged_pax, cells = cells_g)
}
saveRDS(obj_list, file.path(RESULTS_DIR, "pax3_nkx61_objects.rds"))
message("PAX3 / NKX6-1 objects: ", paste(sprintf("%s = %s cells", names(obj_list),
                                                 vapply(obj_list, function(o) format(ncol(o)), "")), collapse = ", "))
message("01_prepare done")
