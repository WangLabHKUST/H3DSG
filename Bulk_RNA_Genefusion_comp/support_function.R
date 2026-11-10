# convenience functions to add/remove "chr" prefix
addChr <- function(contig) {
  ifelse(contig == "MT", "chrM", paste0("chr", contig))
}
removeChr <- function(contig) {
  sub("^chr", "", sub("^chrM", "MT", contig, perl=T), perl=T)
}

# convenience function to check if a value is between two others
between <- function(value, start, end) {
  value >= start & value <= end
}
fusion_rbind <- function(fusions){
  if (colnames(fusions)[1] == "X.gene1") { # Arriba output
    colnames(fusions)[colnames(fusions) %in% c("X.gene1", "strand1.gene.fusion.", "strand2.gene.fusion.")] <- c("gene1", "strand1", "strand2")
    fusions$display_contig1 <- sub(":[^:]*$", "", fusions$breakpoint1, perl=T)
    fusions$display_contig2 <- sub(":[^:]*$", "", fusions$breakpoint2, perl=T)
    fusions$contig1 <- removeChr(fusions$display_contig1)
    fusions$contig2 <- removeChr(fusions$display_contig2)
    fusions$breakpoint1 <- as.numeric(sub(".*:", "", fusions$breakpoint1, perl=T))
    fusions$breakpoint2 <- as.numeric(sub(".*:", "", fusions$breakpoint2, perl=T))
    fusions$split_reads1 <- fusions$split_reads1
    fusions$split_reads2 <- fusions$split_reads2
    fusions$type <- sub(".*(translocation|duplication|deletion|inversion).*", "\\1", fusions$type, perl=T)
    fusions$fusion_transcript <- gsub("[()^$]", "", fusions$fusion_transcript)
    return(fusions)
  }
}