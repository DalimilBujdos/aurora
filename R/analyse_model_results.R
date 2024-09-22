#' Creates the final data frame
#'
#' @noRd
analyse_model_results <- function(df_final, phenotypes, lookup_tbl) {
  index <- c()
  col_names <- c()
  med_cols <- 1:length(phenotypes)
  iqr_cols <- (length(phenotypes)+1):(length(phenotypes)*2)
  probs_cols <- (length(phenotypes)*2 + 1):(length(phenotypes)*3)
  for (i in 1:length(phenotypes)) {
    med <- paste0(lookup_tbl[i], "_classification_median")
    iqr1 <- paste0(lookup_tbl[i], "_classification_interquartile_range")
    probs1 <- paste0(lookup_tbl[i], "_classification_label")
    col_names <- c(col_names, med, iqr1, probs1)
    index <- c(index, med_cols[i], iqr_cols[i], probs_cols[i])
  }
  df_final <- as.data.frame(df_final[,index])
  colnames(df_final) <- col_names
  df_final <- cbind(as.data.frame(rownames(df_final)), df_final)
  names(df_final)[1] <- "strain_ID"
  df_final$observed_phenotype <- lookup_tbl[match(pheno_mat$pheno, names(lookup_tbl))]

  # which phenotype was predicted
  df_final$predicted_phenotype <- c(rep(NA, nrow(df_final)))
  for (i in 1:nrow(df_final)) {
    prob_labels <- df_final[i, seq(from = 1, to = ncol(df_final), by = 3)[-1]]
    true_phenotypes <- which(prob_labels == TRUE)
    if (length(true_phenotypes) == 0) {
      df_final$predicted_phenotype[i] <- "INCONCLUSIVE"
    } else if (length(true_phenotypes) == 1) {
      df_final$predicted_phenotype[i] <- lookup_tbl[true_phenotypes]
    } else {
      df_final$predicted_phenotype[i] <- paste(lookup_tbl[true_phenotypes],
                                               collapse = " and ")
    }
  }

  return(df_final)
}

