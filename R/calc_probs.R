#' Process the probabilites calculated in Outlier calculation phase
#'
#' @noRd
calc_probs <- function(strain,
                       pheno_mat,
                       probs,
                       fin_proxy_mat,
                       fin_proxy_mat_orig) {
  # function analyses probs from outlier calculation


  # first check that the strains main prob is higher than probs of strains misslabelled into
  # this pheno
  own_pheno <- pheno_mat$pheno[pheno_mat$ids == strain]
  # retrieve own probs of the strain
  own_probs <- probs[[which(names(probs) == strain)]]
  own_probs <- own_probs[, as.numeric(own_pheno)]
  p_vals_other <- rep(NA, length(phenotypes))

  # calculate bonfferoni treshold value
  bonffer_p_val <- 0.05/(nrow(pheno_mat)*2) # for every strain we have are testing tow hypotheses

  for (i in 1:length(phenotypes)) {
    if (own_pheno == phenotypes[i]) {next}
    # retrieve all probs of misslabelled strains
    miss_prob <- fin_proxy_mat[phenotypes[i],own_pheno, ]
    # it matters which label is alphabetically higher
    df_MW <- data.frame(probs = c(own_probs, miss_prob),
                        label = as.factor(c(rep("a_label", length(own_probs)),
                                            rep("b_label", length(miss_prob)))) )
    # statistical test
    MW_result <- wilcox.test(df_MW$probs~df_MW$label, alternative = "greater")
    p_vals_other[i] <- MW_result$p.value
  }
  # in the next step we will check that that probs in other phenos then own_pheno
  # are smaller then probs in fin_proxy_mat_orig
  other_probs <- probs[[which(names(probs) == strain)]]
  p_vals_other2 <- rep(NA, length(phenotypes))
  for (i in 1:length(phenotypes)) {
    if (own_pheno == phenotypes[i]) {next}
    other_probs_pheno <- other_probs[, as.numeric(phenotypes[i])]
    miss_prob <- fin_proxy_mat_orig[as.numeric(phenotypes[i]),as.numeric(own_pheno),]

    df_MW <- data.frame(probs = c(other_probs_pheno, miss_prob),
                        label = as.factor(c(rep("a_label", length(own_probs)),
                                            rep("b_label", length(miss_prob)))) )
    # statistical test
    MW_result <- wilcox.test(df_MW$probs~df_MW$label, alternative = "less")
    p_vals_other2[i] <- MW_result$p.value
  }

  # if both p_vals are below bonffer_p_val then strains should be autochthonous in its own_pheno
  # if one of the p_vals is above bonffer_p_val then the result is INCONCLUSIVE
  # if both p_vals are above bonffer_p_val then the result is that the strain should belong
  # to this pheno
  return_label <- rep(NA, length(phenotypes))
  if (p_vals_other[!is.na(p_vals_other)] < bonffer_p_val &&
      p_vals_other2[!is.na(p_vals_other2)] < bonffer_p_val) {
    return_label[as.numeric(own_pheno)] <- "TRUE"
    return_label[-as.numeric(own_pheno)] <- "FALSE"
    return(return_label)
  }

  return_label[as.numeric(own_pheno)] <- "INCONCLUSIVE"
  for (i in 1:length(phenotypes)) {
    if (own_pheno == phenotypes[i]) {next}
    if (p_vals_other[i] > bonffer_p_val & p_vals_other2[i] > bonffer_p_val) {
      return_label[i] <- "TRUE"
    } else if (p_vals_other[i] < bonffer_p_val & p_vals_other2[i] < bonffer_p_val) {
      return_label[i] <- "FALSE"
    } else {
      return_label[i] <- "INCONCLUSIVE"
    }
  }
  return(return_label)
}


