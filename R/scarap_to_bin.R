#' Retrieves data from SCARAP
#'
#' @noRd
scarap_to_bin <- function(bin_mat, pheno_mat) {
  if (is.character(bin_mat)) {
    bin_mat <- as.data.frame(readr::read_delim(bin_mat, show_col_types = FALSE,
                                               delim = "\t", col_names = c("gene_ID", "sample_ID", "gene_family")))
  } else {
    colnames(bin_mat) <- c("gene_ID", "sample_ID", "gene_family")
  }
  gfs <- unique(bin_mat$gene_family)
  strains <- unique(bin_mat$sample_ID)

  bin_mat_mat <- matrix(0, nrow = length(strains), ncol = length(gfs))
  rownames(bin_mat_mat) <- strains
  colnames(bin_mat_mat) <- gfs

  for (i in 1:nrow(bin_mat_mat)) {
    genes <- bin_mat$gene_family[bin_mat$sample_ID == rownames(bin_mat_mat)[i]]
    bin_mat_mat[i,which(is.element(colnames(bin_mat_mat), genes))] <- 1
  }

  bin_mat <- as.data.frame(bin_mat_mat)

  bin_mat <- bin_mat[match(pheno_mat$ids, rownames(bin_mat)), ]
  return(bin_mat)
}
