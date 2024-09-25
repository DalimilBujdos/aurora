#' Retrieves data from tool DRAM
#'
#' @noRd
get_DRAM <- function(filename, pheno_mat) {
  sheets <- readxl::excel_sheets(filename)
  sheets <- sheets[!is.element(sheets, c("rRNA", "tRNA"))] # we dont want to include rRNA and rRNA cols
  df_raw <- data.frame(strain_id = pheno_mat$ids)
  for (i in 1:length(sheets)) {
    df_tmp <- readxl::read_excel(filename, sheet = sheets[i])
    df_tmp <- df_tmp %>% distinct(gene_id, .keep_all = TRUE)
    new_cols <- paste(as.character(df_tmp$gene_id), sheets[i], sep = "_")
    df_tmp <- df_tmp %>%
      dplyr::mutate_all(~replace(., .>0, 1))
    df_tmp <- t(as.data.frame(df_tmp[,-1:-5]))
    colnames(df_tmp) <- new_cols
    # check that the strains that are in pheno_mat are also in DRAM output
    if (i == 1) {
      diff <- setdiff(rownames(df_tmp), pheno_mat$ids)
      if (length(diff) > 0) {
        warning("DRAM output contains strains that are not in pheno_mat.",
                " Strains that are not in pheno_mat were removed.\n")
      }

      diff <- setdiff(pheno_mat$ids, rownames(df_tmp))
      if (length(diff) > 0) {
        messa <- paste("pheno_mat contains strains:", paste(diff, collapse = " "), "that are not in DRAM output.", sep = " ")
        stop(messa, "\n")
      }
    }
    df_tmp <- df_tmp[match(pheno_mat$ids, rownames(df_tmp)), ]
    df_raw <- cbind(df_raw, as.data.frame(df_tmp))
  }

  df_raw <- as.data.frame(sapply(df_raw[, -1], as.numeric))
  rownames(df_raw) <- pheno_mat$ids
  return(df_raw)
}
