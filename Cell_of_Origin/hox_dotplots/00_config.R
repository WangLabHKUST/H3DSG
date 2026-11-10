# =============================================================================
# hox_dotplots / 00_config.R -- paths, constants and sample tables shared by the scripts
#
# Dot plots of the position marker gene sets (thalamic, brainstem, cervical and thoracic
# HOX / ZIC / OTX1 genes) and of PAX3 / NKX6-1 in immune, stromal and malignant cells:
#
#   01_prepare.R                    reads the integrated object once, keeps what the statistics need
#   02_compute_stats.R              average expression / percentage of expressing cells -> ../data/hox_dotplots/*.tsv
#   03_plot_location_dotplots.R     hox_location_dotplot_by_group.pdf, hox_location_dotplot_by_sample.pdf
#   04_plot_pax3_nkx61.R            pax3_nkx61_malignant_dotplot.pdf, pax3_nkx61_celltype_dotplot.pdf
#
# Run each script from this folder (it uses relative paths and never changes the working
# directory), in this order:
#
#   Rscript 01_prepare.R           # reads the 5.6 GB integrated object, about 2 min
#   Rscript 02_compute_stats.R
#   Rscript 03_plot_location_dotplots.R
#   Rscript 04_plot_pax3_nkx61.R
#
# The plotting scripts read only the tables in ../data/hox_dotplots/ and open no object, so
# 03 and 04 can be run without 01 and 02.
#
# Inputs (read only):
#   ../../Single_Cell_RNAseq/objects/integrated/integrated_minf400_only_PL_6_4_cca_annotated.rds
#                                   the integrated object of the 19 samples (132,343 cells)
#   ../../Single_Cell_RNAseq/sample_metadata.tsv
#   ../../Single_Cell_RNAseq/inputs/hox/loc_gene.csv                  gene -> anatomical level
#   ../../Single_Cell_RNAseq/inputs/hox/location_obj_list2.rds        public brain DMG data (GSE210568),
#   ../../Single_Cell_RNAseq/inputs/hox/location_objects_full.rds     thalamic and pontine tumours; the
#                                   second file holds all genes (used for PAX3 / NKX6-1), the first 34 genes
#   ../data/hox_dotplots/sample_location.tsv                          location of each sample (explicit table)
# =============================================================================

# ---- paths (relative to this folder) -----------------------------------------
SCRNA_DIR       <- file.path("..", "..", "Single_Cell_RNAseq")
INTEGRATED_RDS  <- file.path(SCRNA_DIR, "objects", "integrated", "integrated_minf400_only_PL_6_4_cca_annotated.rds")
SAMPLE_META_TSV <- file.path(SCRNA_DIR, "sample_metadata.tsv")
LOC_GENE_CSV    <- file.path(SCRNA_DIR, "inputs", "hox", "loc_gene.csv")
BRAIN_LIST_RDS  <- file.path(SCRNA_DIR, "inputs", "hox", "location_obj_list2.rds")
BRAIN_FULL_RDS  <- file.path(SCRNA_DIR, "inputs", "hox", "location_objects_full.rds")

DATA_DIR    <- file.path("..", "data", "hox_dotplots")   # tables read by the plotting scripts
RESULTS_DIR <- "results"                                  # intermediate objects, not committed

# ---- constants ------------------------------------------------------------------
CELLTYPES <- c("immune", "stromal", "malignant")
# integrated_snn_res.1_anno (annotation of the integrated object) -> three groups
CELLTYPE_MAP <- c(myeloid = "immune", neutrophil = "immune", "T/NK" = "immune", B = "immune",
                  vascular = "stromal", oligo = "stromal",
                  tumor = "malignant")
GENE_SET_LEVELS <- c("Thalamus", "Brainstem", "Cervical", "Thoracic")   # anatomical level marked by a gene set
PCT_METHOD  <- "any"   # percentage = cells with at least one gene of the set detected
COLOR_MAX   <- 0.16    # upper end of the colour bar of the by-sample figure (average expression)
SPLIT_CMAX  <- 0.12    # average expression is clipped at this value in the by-group figure

# ---- location of each sample -----------------------------------------------------
# ../data/hox_dotplots/sample_location.tsv gives the anatomical group of each sample:
# Thoracic if the caudal end of the resected tumour is in a thoracic segment (T1 or below),
# Cervical otherwise (the cranial-most tumours start at the medulla, M-C5). It is checked
# against the tumour location of sample_metadata.tsv in read_sample_table().
#
# One sample is counted in a different group in the by-group figure (fixed rule, no switch):
# TL_4. Its tumour spans C5-T1, so by the caudal end it is a thoracic sample, but it only just
# reaches T1 and none of its 34 stromal cells expresses a gene of the Thoracic set (0 of 34).
# Counted with the thoracic samples, those 34 zero cells would dilute the Thoracic-TL stromal
# dot (44 stromal cells: 2.3% with TL_4; 10 cells, 10% with TL_1 alone), so TL_4 is counted
# with the cervical samples here (Cervical-TL, which then has 10 samples; Thoracic-TL is TL_1
# alone). The by-sample figure shows every sample on its own and is not affected.
PANEL_LOCATION_OVERRIDE <- c(TL_4 = "Cervical")

# ---- panel order of the by-sample figure -----------------------------------------
# The tumour location is read as an interval of spinal segments and the samples are
# ordered by caudal end, then cranial end, then sample id.
# M (medulla) = 0, C1..C8 = 1..8, T1..T12 = 9..20, L1..L5 = 21..25; the suffixes
# Dorsal / Ventral / edge are ignored.
loc_span <- function(loc) {
  n <- vapply(toupper(strsplit(as.character(loc), "-")[[1]]), function(t) {
    if (identical(t, "M")) return(0)
    m <- regmatches(t, regexec("^([CTLS])([0-9]+)$", t))[[1]]
    if (!length(m)) return(NA_real_)
    unname(c(C = 0, T = 8, L = 20, S = 25)[m[2]] + as.numeric(m[3]))
  }, numeric(1))
  n <- n[!is.na(n)]
  if (!length(n)) stop("cannot read a spinal segment from: ", loc, call. = FALSE)
  c(min(n), max(n))
}

# The 19 samples with their group, tumour location, anatomical group and panel group.
# Returns one row per sample, in the panel order of the by-sample figure (column panel_order).
read_sample_table <- function() {
  meta <- read.delim(SAMPLE_META_TSV, colClasses = "character")
  loc  <- read.delim(file.path(DATA_DIR, "sample_location.tsv"), colClasses = "character")
  stopifnot(nrow(meta) == 19, !anyDuplicated(meta$sid), setequal(meta$sid, loc$sid),
            all(loc$location %in% c("Cervical", "Thoracic")))
  tab <- data.frame(sid = meta$sid, renamed_id = meta$renamed_id, patient_id = meta$patient_id,
                    subtype = meta$group, tumor_location = meta$tumor_location,
                    location = loc$location[match(meta$sid, loc$sid)], stringsAsFactors = FALSE)
  stopifnot(all(tab$subtype %in% c("PL", "TL")))
  span <- vapply(tab$tumor_location, loc_span, numeric(2))
  stopifnot(identical(tab$location == "Thoracic", unname(span[2, ] > 8)))   # caudal end in T1 or below
  stopifnot(all(names(PANEL_LOCATION_OVERRIDE) %in% tab$sid))
  tab$panel_location <- tab$location
  tab$panel_location[match(names(PANEL_LOCATION_OVERRIDE), tab$sid)] <- unname(PANEL_LOCATION_OVERRIDE)
  ord <- order(span[2, ], span[1, ], tab$sid)                              # caudal end, cranial end, sample id
  tab <- tab[ord, , drop = FALSE]
  tab$panel_order <- seq_len(nrow(tab))
  rownames(tab) <- NULL
  tab
}

# ---- colour bar of the dot plots -----------------------------------------------------
# cmax = NA: the colour bar follows the data; otherwise it is fixed to c(0, cmax) and
# larger values are squished to the top colour
avg_color_scale <- function(name = "Average Expression", cmax = COLOR_MAX) {
  if (is.na(cmax)) return(ggplot2::scale_color_gradient2(low = "lightgray", high = "#1a08fe", name = name))
  ggplot2::scale_color_gradient2(low = "lightgray", high = "#1a08fe", name = name,
                                 limits = c(0, cmax), oob = scales::squish,
                                 breaks = ggplot2::waiver(), labels = ggplot2::waiver())
}

# Theme without any text, tick marks kept: the by-sample figure is drawn without labels
# (the sample labels are added when the figure is assembled).
# The sub-elements axis.text.x / axis.text.y are blanked explicitly: DotPlotDF sets
# axis.text.x itself, and blanking only the parent element does not override it.
nude_theme <- function() ggplot2::theme(
  axis.title       = ggplot2::element_blank(),
  axis.title.x     = ggplot2::element_blank(),
  axis.title.y     = ggplot2::element_blank(),
  axis.text        = ggplot2::element_blank(),
  axis.text.x      = ggplot2::element_blank(),
  axis.text.y      = ggplot2::element_blank(),
  strip.text       = ggplot2::element_blank(),
  strip.text.x     = ggplot2::element_blank(),
  strip.text.y     = ggplot2::element_blank(),
  strip.background = ggplot2::element_blank(),
  legend.position  = "none",
  plot.title       = ggplot2::element_blank(),
  plot.subtitle    = ggplot2::element_blank(),
  plot.caption     = ggplot2::element_blank())

# Read a table written by 02_compute_stats.R
read_table <- function(file) {
  path <- file.path(DATA_DIR, file)
  if (!file.exists(path)) stop("missing ", path, ": run 02_compute_stats.R first", call. = FALSE)
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}
