#' Postprocess data from ML models in Threshold calculation phase
#'
#' @noRd
postprocess_fit_miss <- function(rf_result,
                                 misslabel_mat,
                                 phenotypes) {
  # initialise the final df
  df_miss <- matrix(data = NA, nrow = nrow(rf_result), ncol = length(phenotypes)+1)
  df_miss <- as.data.frame(df_miss)
  # puts the result from a classificator into one dataframe and with a correct class label
  for (i in 1:length(phenotypes)) {
    miss_strains <- misslabel_mat[rownames(misslabel_mat) == phenotypes[i], ]
    names(miss_strains) <- rownames(misslabel_mat)
    for (j in 1:length(phenotypes)) {
      if(i == j) {next}
      label <- paste0(phenotypes[i], "to", names(miss_strains)[j])
      index <- which(rownames(rf_result) == miss_strains[j])
      values <- rf_result[index, 1:length(phenotypes)]
      # append results to dataframe
      df_miss[which(is.na(df_miss[,1]))[1], ] <- c(as.numeric(values), label)
    }
  }
  # now put the non-misslabelled strains into the same df
  remove_strains <- as.character(misslabel_mat)
  remove_strains <- remove_strains[!is.na(remove_strains)]
  data_nonmiss <- rf_result[!is.element(rownames(rf_result), remove_strains), ]
  df_miss[which(is.na(df_miss[,1]))[1]:nrow(df_miss), ] <- data_nonmiss

  return(df_miss)
}
