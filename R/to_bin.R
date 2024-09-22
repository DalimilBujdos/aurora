#' Initial filtering of the binary matrix
#'
#' @noRd
to_bin <- function(bin_mat,
                   pheno_mat,
                   low_perc_cutoff = 3,
                   upp_perc_cutoff = 99,
                   run_chisq = FALSE,
                   cutoff_chisq = 0.1) {

  # this function takes path to presence/absence file produced by panaroo and
  # returns the same filtered datagrame but the data are encoded as 1/0
  no_col_initially <- ncol(bin_mat)

  no_gf <- nrow(bin_mat)
  cutoff <- no_gf * (upp_perc_cutoff/100)
  select_ID_up <- which(colSums(bin_mat) < cutoff)

  # here do the same but exclude gene families with freq lower than cutoff
  cutoff <- no_gf * (low_perc_cutoff/100)
  non_select_ID_low <- which(colSums(bin_mat) > cutoff)

  bin_mat <- bin_mat[, intersect(select_ID_up, non_select_ID_low)]

  ncol_after_filter <- ncol(bin_mat)

  cat("Initial filtering removed", no_col_initially - ncol_after_filter, "variants. There are",
      ncol(bin_mat), "variants left \n")

  # Choose only the gene families that have chi-sq that can differentiate
  # between the two groups

  if (run_chisq == TRUE) {

    calculate_chisq <- function(x, bin_mat, pheno_mat) {
      # calculates chisq test for presence/absence pattern of a gene family
      conting_mat <- matrix(data = NA, nrow = 2, ncol = nlevels(as.factor(pheno_mat$pheno)))
      rownames(conting_mat) <- c("present", "absent")
      lvls <- levels(as.factor(pheno_mat$pheno))
      colnames(conting_mat) <- lvls

      datatest <- bin_mat[, colnames(bin_mat) == x]

      for (i in 1:length(lvls)) {
        # get ids
        ids <- pheno_mat$ids[pheno_mat$pheno == lvls[i]]
        x_lvl <- datatest[is.element(rownames(bin_mat), ids)]
        # write the contingency table
        conting_mat[1,i] <- sum(x_lvl) # gfs present
        conting_mat[2,i] <- length(x_lvl) - sum(x_lvl) # gfs absent
      }
      result_chisq <- stats::chisq.test(conting_mat, simulate.p.value = TRUE)
      return(result_chisq$p.value)
    }

    p_vals <- sapply(colnames(bin_mat), calculate_chisq,
                     bin_mat = bin_mat,
                     pheno_mat = pheno_mat)
    # if the chisq p_value is above the treshold remove it
    bin_mat <- bin_mat[, p_vals < cutoff_chisq]
    no_removed <- length(which((p_vals < cutoff_chisq) == FALSE))
    cat("Chi-square filter with the treshold p-value =", cutoff_chisq,
        "removed", no_removed, "variants \n")
  }

  bin_mat <- bin_mat[match(pheno_mat$ids, rownames(bin_mat)), ]

  return(list(bin_mat = bin_mat,
              no_removed = no_col_initially - ncol_after_filter,
              chisq_removed = if (exists("no_removed")) no_removed else FALSE))
}
