#' Retrieves data from unitig_counter
#'
#' @noRd
#'
unitigs_to_bin <- function(bin_mat_unitigs, pheno_mat) {
  gz_file <- gzfile(bin_mat_unitigs, "rt")
  raw_text <- readLines(gz_file)
  close(gz_file)

  make_bin <- function(x, pheno_mat) {
    x <-  strsplit(x, " \\| ")[[1]]
    ids <- regmatches(x[2], gregexpr("strain_\\d+", x[2]))
    ids <- unlist(ids)
    final_vec <- rep(0, nrow(pheno_mat))
    final_vec[is.element(pheno_mat$ids, ids)] <- 1
    return(list(name = x[1], vec = final_vec))
  }

  list_new_rows <- lapply(raw_text, make_bin, pheno_mat = pheno_mat)
  names_list <- lapply(list_new_rows, function(x) x$name)
  vecs_list <- lapply(list_new_rows, function(x) x$vec)

  # Combine the 'vec' elements into a matrix and transpose it to have rows as columns
  bin_mat <- do.call(cbind, vecs_list)

  # Create a dataframe using the combined matrix and set column names
  bin_mat <- as.data.frame(bin_mat)
  bin_mat <- apply(bin_mat, 1, as.numeric)
  bin_mat <- t(bin_mat)
  colnames(bin_mat) <- unlist(names_list)
  rownames(bin_mat) <- pheno_mat$ids
  return(as.data.frame(bin_mat))
}
