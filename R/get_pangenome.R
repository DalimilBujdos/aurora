#' Retrieves data from Roary or Panaroo
#'
#' @noRd
get_pangenome <- function(bin_mat, pheno_mat) {
  if (is.character(bin_mat)) {
    bin_mat <- as.data.frame(readr::read_csv(bin_mat, show_col_types = FALSE))
  }

  annotations <- bin_mat$Annotation
  names(annotations) <- bin_mat$Gene

  fut_row_names <- bin_mat$Gene
  bin_mat <- bin_mat[, -1:-14]
  bin_mat <- bin_mat %>%
    dplyr::mutate_all(~replace(., is.na(.), 0)) %>%
    dplyr::mutate_all(~replace(., .!=0, 1))
  bin_mat <- sapply(bin_mat, as.numeric)
  bin_mat <- t(bin_mat)
  colnames(bin_mat) <- fut_row_names

  bin_mat <- bin_mat[match(pheno_mat$ids, rownames(bin_mat)), ]
  return(list(bin_mat = bin_mat, annotations = annotations))
}
