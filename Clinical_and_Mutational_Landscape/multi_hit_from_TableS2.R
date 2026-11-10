# Cells of the mutational landscape where one sample carries two or more variants
# in the same gene, taken from Table S2 (hg38, SAVI).
#
# The landscape draws one mutation type per cell (from Data_clinical_info.csv), so
# a second type in the same cell is invisible. This script finds those cells and
# writes the hidden type out, to be drawn as a right triangle on the cell: the
# lower left half gets the colour of the second type and the upper right half
# keeps the colour the landscape already draws.
#
# It also checks every gene x sample cell both ways and stops if the two sources
# disagree: the type the landscape draws has to be one of the types Table S2 lists
# for that cell, and a cell is either called in both sources or in neither. That
# check is what catches a clinical table that has fallen behind the variant table.
#
# Output: landscape_multi_hit.csv
rm(list = ls())
if (!require('readxl'))
  install.packages('readxl')

# run from RStudio, or with Rscript from anywhere

S2_PATH <- 'Table S2.xlsx'
# these have to stay in step with alter_fun in Clinical_and_Mutational_landscape.R
MUT_COL <- c(Mut_Missense = 'darkturquoise', Mut_Frameshift = 'sienna4',
             Mut_Splice = 'hotpink1', Mut_Stop_gain = 'black',
             Mut_Inframe = 'darkturquoise')
hex <- function(cols) toupper(rgb(t(col2rgb(cols)), maxColorValue = 255))

## ---- Table S2 -> one row per variant ----
# Gene_Name can be a compound name such as "NF1, AC..." and is split on commas
s2 <- as.data.frame(readxl::read_excel(S2_PATH, skip = 1))
stopifnot(c('Patient ID', 'Gene_Name', 'Effect', 'Amino_Acid_Change',
            'Position (hg38)', 'Tumor_freq') %in% names(s2))
genes_split <- strsplit(as.character(s2$Gene_Name), '\\s*,\\s*')
v <- s2[rep(seq_len(nrow(s2)), lengths(genes_split)), ]
v$gene <- unlist(genes_split)
v$P    <- sub('^SP', 'P', v[['Patient ID']])

# SnpEff effect -> landscape type, using the first transcript's effect
# (truncating > splice site > inframe > missense > splice region)
effect_type <- function(e) {
  e <- sub(',.*', '', e)
  ifelse(grepl('frameshift', e), 'Mut_Frameshift',
  ifelse(grepl('stop_gained', e), 'Mut_Stop_gain',
  ifelse(grepl('splice_(donor|acceptor)', e), 'Mut_Splice',
  ifelse(grepl('inframe', e), 'Mut_Inframe',
  ifelse(grepl('missense', e), 'Mut_Missense',
  ifelse(grepl('splice_region', e), 'Mut_Splice', NA))))))
}
v$type <- effect_type(v$Effect)

## ---- the landscape data ----
clin  <- read.csv('Data_clinical_info.csv', check.names = FALSE, row.names = 1,
                  colClasses = 'character')
genes <- rownames(clin)[which(rownames(clin) == 'H3-3A'):which(rownames(clin) == 'FGFR1')]
v <- v[v$gene %in% genes & !is.na(v$type), ]

## ---- check every cell, then collect the multi-hit ones ----
bad <- character()
rows <- list()
for (g in genes) for (p in colnames(clin)) {
  drawn <- clin[g, p]
  drawn <- if (is.na(drawn) || drawn %in% c('', 'NA')) NA else drawn
  vv <- v[v$gene == g & v$P == p, ]
  types <- unique(vv$type)
  if (is.na(drawn) != (nrow(vv) == 0) || (!is.na(drawn) && !drawn %in% types)) {
    bad <- c(bad, sprintf('%s %s: landscape=%s, Table S2=%s',
                          g, p, drawn, paste(types, collapse = '/')))
  }
  if (nrow(vv) >= 2) {
    other <- setdiff(types, drawn)
    if (length(other) > 1)
      stop(sprintf('%s %s has more than two types, a triangle cannot show them: %s',
                   g, p, paste(types, collapse = '/')))
    rows[[length(rows) + 1]] <- data.frame(
      sample = p, SP = sub('^P', 'SP', p), gene = g, n_variants = nrow(vv),
      drawn_type = drawn, drawn_color = hex(MUT_COL[[drawn]]),
      triangle_type = if (length(other)) other else NA,
      triangle_color = if (length(other)) hex(MUT_COL[[other]]) else NA,
      variants = paste(sprintf('%s %s (%s, VAF %s%%)',
                               sub(',.*', '', vv$Amino_Acid_Change),
                               sub('Mut_', '', vv$type),
                               vv[['Position (hg38)']], vv$Tumor_freq), collapse = '; '),
      stringsAsFactors = FALSE)
  }
}
if (length(bad))
  stop('the landscape data and Table S2 disagree:\n  ', paste(bad, collapse = '\n  '))
cat(sprintf('check passed: %d genes x %d samples, every type the landscape draws is in Table S2\n',
            length(genes), ncol(clin)))

res <- do.call(rbind, rows)
res <- res[order(match(res$gene, genes), res$sample), ]
write.csv(res, 'landscape_multi_hit.csv', row.names = FALSE)
print(res[, c('sample', 'gene', 'n_variants', 'drawn_type', 'drawn_color',
              'triangle_type', 'triangle_color')], row.names = FALSE)
cat('\nwrote landscape_multi_hit.csv\n')
