# =============================================================================
# hox_dotplots / hox_helpers.R -- gene set statistics for the dot plots
#   df_to_genesets                                      gene table -> list of gene sets
#   calculate_gene_set_stats_advanced                   avg.exp / pct.exp of each gene set
#                                                       in each group of one Seurat object
#   calculate_multi_samp_gene_set_stats_advanced        the same for a list of objects
#
# avg.exp = mean of the (log-normalised, "data" layer) expression over all genes of the set
#           and all cells of the group, zeros included.
# pct.exp = percentage of the cells of the group in which at least one gene of the set is
#           detected (pct_method = "any"); for a gene set of one gene this is the usual
#           percentage of cells expressing the gene.
# =============================================================================
library(Seurat)
library(dplyr)
library(rlang)

df_to_genesets <- function(df, gene_col = "gene_name", set_loc_col = "Location") {
  df <- df %>%
    dplyr::group_by(!!rlang::sym(set_loc_col)) %>%
    dplyr::summarise(genes = list(!!rlang::sym(gene_col)), .groups = "drop")
  out <- list()
  for (i in seq_len(nrow(df))) {
    out[[as.character(df[[set_loc_col]][i])]] <- df$genes[[i]]
  }
  out
}

calculate_gene_set_stats_advanced <- function(
    seurat_obj,
    gene_sets,
    group.by = "celltype2",
    layer = "data",
    pct_method = "any",
    avg_method = "mean",
    pos.only = FALSE) {
  results <- data.frame(stringsAsFactors = FALSE)
  expr_layer <- GetAssayData(seurat_obj, layer = layer)
  cell_types <- seurat_obj@meta.data[[group.by]]

  for (set_name in names(gene_sets)) {
    genes_in_set <- intersect(gene_sets[[set_name]], rownames(seurat_obj))
    ## genes absent from this object -> keep the set but report zero expression
    ## (instead of dropping it) so every requested gene still appears downstream
    has_genes <- length(genes_in_set) > 0
    if (has_genes) {
      expr_matrix <- expr_layer[genes_in_set, , drop = FALSE]
    }
    for (cell_type in unique(cell_types)) {
      cells_of_type <- colnames(seurat_obj)[cell_types == cell_type]
      if (!length(cells_of_type)) {
        next
      }
      if (!has_genes) {
        results <- rbind(
          results,
          data.frame(
            features.plot = set_name,
            id = cell_type,
            avg.exp = 0,
            pct.exp = 0,
            stringsAsFactors = FALSE
          )
        )
        next
      }
      expr_subset <- expr_matrix[, cells_of_type, drop = FALSE]

      if (pct_method == "any") {
        pct_exp <- mean(colSums(expr_subset > 0) > 0) * 100
      } else if (pct_method == "all") {
        pct_exp <- mean(colSums(expr_subset > 0) == nrow(expr_subset)) * 100
      } else {
        pct_exp <- mean(colMeans(expr_subset > 0)) * 100
      }

      if (avg_method == "mean") {
        if (pos.only) {
          vals <- as.vector(expr_subset)
          vals <- vals[vals > 0]
          avg_exp <- if (length(vals)) mean(vals) else 0
        } else {
          avg_exp <- mean(as.vector(expr_subset))
        }
      } else if (avg_method == "max") {
        if (pos.only) {
          max_per_cell <- apply(expr_subset, 2, function(x) if (any(x > 0)) max(x[x > 0]) else NA_real_)
          avg_exp <- if (any(!is.na(max_per_cell))) mean(max_per_cell, na.rm = TRUE) else 0
        } else {
          avg_exp <- mean(apply(expr_subset, 2, max))
        }
      } else {
        if (pos.only) {
          sum_per_cell <- apply(expr_subset, 2, function(x) sum(x[x > 0]))
          avg_exp <- mean(sum_per_cell)
        } else {
          avg_exp <- mean(colSums(expr_subset))
        }
      }

      results <- rbind(
        results,
        data.frame(
          features.plot = set_name,
          id = cell_type,
          avg.exp = avg_exp,
          pct.exp = pct_exp,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  results
}

calculate_multi_samp_gene_set_stats_advanced <- function(
    seurat_obj_lst,
    gene_sets,
    group.by = "celltype2",
    layer = "data",
    sample_group.by = "orig.ident",
    pct_method = "any",
    avg_method = "mean",
    pos.only = FALSE) {
  out <- data.frame(check.names = FALSE)
  for (sample_id in names(seurat_obj_lst)) {
    obj_this <- seurat_obj_lst[[sample_id]]
    dot_data <- calculate_gene_set_stats_advanced(
      obj_this,
      gene_sets,
      group.by = group.by,
      layer = layer,
      pct_method = pct_method,
      avg_method = avg_method,
      pos.only = pos.only
    )
    dot_data[[sample_group.by]] <- sample_id
    out <- rbind(out, dot_data)
  }
  out
}
