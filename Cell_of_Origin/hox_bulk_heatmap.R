# Extended Data Fig. 2a: expression of 39 HOX genes in bulk RNA-seq, log2(RPKM + 1).
# Columns: 3 thalamus DMG, 2 pons DMG (brain) and 40 spinal cord tumours (H3-DSG-PL, H3-DSG-TL), ordered by tumour location.
#
# Input (data/):
#   Spinalcord.H3K27M.expression.RPKM.txt     RPKM of 40 spinal cord tumours (columns are SP IDs), hg19
#   H3K27altered_project_metadata_updated.tsv project metadata; TSNE.cluster A = pons-like (PL), B = thalamus-like (TL)
#   SP45_ID_location_detail.txt               tumour location of each spinal cord tumour (Location, Location_2)
#   CGGA.H3K27M.hg38.expression.RPKM.txt      RPKM of H3 K27M gliomas of the CGGA cohort (public, hg38); three thalamic ones are used
#   StJude.DIPG.RPKM.txt                      RPKM of pons DMGs of the St Jude PeCan cohort (public); two are used
# Output: hox_bulk_heatmap.pdf
#
# Grouping: the row "type" follows the project metadata (TSNE.cluster A / B), as in the submitted figure.
# The Source Data of the figure follow Table S1, which moves SP02 and SP42 from TL to PL.
# Run from the Cell_of_Origin/ directory:  Rscript hox_bulk_heatmap.R

library(pheatmap)

# ---- data ----
ge <- read.delim("data/CGGA.H3K27M.hg38.expression.RPKM.txt", row.names = 1)          # thalamus DMG
df <- read.delim("data/Spinalcord.H3K27M.expression.RPKM.txt")                        # spinal cord tumours
df.cl <- read.delim("data/H3K27altered_project_metadata_updated.tsv")

df.pons     <- df[, names(df) %in% df.cl$ID2[df.cl$TSNE.cluster == "A"]]
df.thalamus <- df[, names(df) %in% df.cl$ID2[df.cl$TSNE.cluster == "B"]]
stopifnot(ncol(df.pons) == 16, ncol(df.thalamus) == 24)                                # all 40 tumours matched

ge.thalamus <- ge[, names(ge) %in% c("CGGA_791", "CGGA_P266", "CGGA_1185", "CGGA_1320", "CGGA_1735")]

ge.pons <- read.delim("data/StJude.DIPG.RPKM.txt")                                    # pons DMG
ge.pons <- ge.pons[!is.na(ge.pons$gene), ]

# ---- expression matrix m: 39 HOX genes x 45 samples ----
gns <- ge.pons$gene[grepl("HOX[ABCD]", ge.pons$gene)]

hox.thalamus_like  <- df.thalamus[gns, ]
hox.pons_like      <- df.pons[gns, ]
hox.thalamus_brain <- ge.thalamus[gns, c("CGGA_1320", "CGGA_1735", "CGGA_791")]
hox.pons_brain     <- ge.pons[ge.pons$gene %in% gns, 1:6]
rownames(hox.pons_brain) <- hox.pons_brain$gene

hoxgns <- c(paste0("HOXA", c(1:7, 9:11, 13)), paste0("HOXB", c(1:9, 13)),
            paste0("HOXC", c(4:6, 8:13)),     paste0("HOXD", c(1, 3, 4, 8:13)))

hox.thalamus_like  <- hox.thalamus_like[hoxgns, ]
hox.pons_like      <- hox.pons_like[hoxgns, ]
hox.thalamus_brain <- hox.thalamus_brain[hoxgns, ]
hox.pons_brain     <- hox.pons_brain[hoxgns, c(4, 5)]                                 # SJHGG102 and SJHGG145

m <- cbind(hox.thalamus_brain, hox.thalamus_like, hox.pons_brain, hox.pons_like)
stopifnot(nrow(m) == 39, ncol(m) == 45, !anyNA(m))

# ---- column annotation and column order ----
m_ann <- data.frame(row.names = names(m),
                    type = c(rep("Brain thalamus", 3),
                             rep("Thalamus-like", 24),
                             rep("Brain pons",     2),
                             rep("Pons-like",     16)))

tumorlocation <- read.delim("data/SP45_ID_location_detail.txt")
m_ann$location  <- tumorlocation$Location[match(rownames(m_ann), tumorlocation$ID2)]
m_ann$location2 <- tumorlocation$Location_2[match(rownames(m_ann), tumorlocation$ID2)]
m_ann$location[m_ann$type == "Brain pons"]      <- "Brain_pons"
m_ann$location[m_ann$type == "Brain thalamus"]  <- "Brain_thalamus"
m_ann$location2[m_ann$type == "Brain pons"]     <- "Brain_pons"
m_ann$location2[m_ann$type == "Brain thalamus"] <- "Brain_thalamus"
m_ann$location2 <- factor(m_ann$location2,
                          levels = c("Brain_thalamus", "Brain_pons", "M-C", "C", "C-T", "M-T1", "T", "T-L"))
stopifnot(!anyNA(m_ann$location), !anyNA(m_ann$location2))

# Columns are sorted by Location_2, then by the exact location. The original script broke ties between tumours with the same
# location by their hospital sample number (not published); that order is fixed here, so the columns are as in the figure.
tie_order <- c("CGGA_1320", "CGGA_1735", "CGGA_791", "SJHGG102", "SJHGG145", "SP51", "SP13", "SP32", "SP07", "SP09", "SP33",
               "SP17", "SP25", "SP30", "SP14", "SP04", "SP10", "SP16", "SP03", "SP22", "SP12", "SP18", "SP43", "SP35", "SP20",
               "SP02", "SP05", "SP41", "SP36", "SP29", "SP15", "SP31", "SP21", "SP06", "SP26", "SP01", "SP11", "SP34", "SP28",
               "SP24", "SP19", "SP23", "SP42", "SP27", "SP08")
stopifnot(setequal(tie_order, rownames(m_ann)))
m_ann <- m_ann[order(m_ann$location2, m_ann$location, match(row.names(m_ann), tie_order)), ]
m_ann$idx <- 1:nrow(m_ann)
# Manual reorder of the original figure: columns 17 and 18 are moved in front of 15 and 16 (within the C tumours).
m_ann <- m_ann[c(1:14, 17, 18, 15:16, 19:45), ]
names(m_ann)[2:3] <- c("location2", "location")   # as in the original script: the exact location and the Location_2 class swap names here

# ---- plot ----
ann_cols <- list(type = c("Brain pons"     = "#f1a23a",
                          "Brain thalamus" = "#95add2",
                          "Pons-like"      = "#990000",
                          "Thalamus-like"  = "#2b8cbe"))

pdf("hox_bulk_heatmap.pdf", width = 481 / 72, height = 283 / 72)                       # page size of the original figure
pheatmap(log2(m[, rownames(m_ann)] + 1),
         color = colorRampPalette(colors = c("darkblue", "#ffff90", "#ffff50", "yellow", "orange", "red"))(50),
         cluster_rows = F, cluster_cols = F, show_colnames = F,
         gaps_row = c(11, 21, 30), border_color = "#dddddd", legend = T,
         annotation_col = subset(m_ann, select = c("type", "location")),
         fontsize_row = 6, treeheight_col = 10, gaps_col = c(3, 5),
         annotation_colors = ann_cols)
dev.off()
