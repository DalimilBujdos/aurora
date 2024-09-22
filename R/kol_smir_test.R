#' Uses Kolmogorov-Smirnov test to calculate if the strain is autochthonous and allochthonous
#'
#' @noRd
kol_smir_test <- function(strain_probs,
                          df_threshold,
                          pheno_mat,
                          phenotypes,
                          lookup_tbl) {
  # function given strains probabiliteis will run kolmogorov-smirnov test to
  # fund strains that were misslabelled

  # null hyphothesis of kolmogorov smirnov test say that strains
  # x and y were drawn from the same continuous distribution
  str_id <- rownames(strain_probs)[1]
  observed_pheno <- pheno_mat$pheno[pheno_mat$ids == str_id]
  test_stat1 <- c()
  test_stat2 <- c()
  origin_probs <- strain_probs[, as.character(observed_pheno)]
  for (i in 1:length(phenotypes)) {
    if (phenotypes[i] == observed_pheno) {
      # tests if probs of strain in pheno 1 is less then probs of other strains
      # in pheno 1. Probs of pheno 1 are being used for the test. If the p-val is low
      # it means that the strain probs are lower then probs of other strains and thus the
      # strain might be allochthonous in pheno 1
      index <- df_threshold$pheno == as.character(observed_pheno)
      test_probs <- as.numeric(df_threshold[index, as.character(observed_pheno)])
      ks_result <- suppressWarnings(stats::ks.test(test_probs, origin_probs, alternative = "less"))
      test_stat1 <- c(test_stat1, ks_result$p.value)
      test_stat2 <- c(test_stat2, NA)
    } else {
      # test that strain with orig pheno 1 is not strain of phenotype X misslabelled
      # to pheno 1. Use probs in pheno 1 to test this. If the p val is low it means
      # that the probs of strain are deffinetly higher then probs of strains misslabelled
      # from pheno X to pheno 1. This means that the strain is autochthonous is pheno 1
      test_label <- paste0(phenotypes[i], "to",observed_pheno)
      index <- df_threshold$pheno == test_label
      test_probs <- as.numeric(df_threshold[index, as.character(observed_pheno)])
      ks_result <- suppressWarnings(stats::ks.test(origin_probs, test_probs, alternative = "less"))
      test_stat1 <- c(test_stat1, ks_result$p.value)
      # do the same here but test with probs from pheno X. If the p val is low it means
      # that the probs of strain in pheno X are much lower then probs of strains X
      # misslabelled into pheno 1. This means that the strain does not belong to pheno X
      test_probs <- as.numeric(df_threshold[index, as.character(phenotypes[i])])
      nonorigin_probs <- strain_probs[, as.character(phenotypes[i])]
      ks_result <- suppressWarnings(stats::ks.test(test_probs, nonorigin_probs, alternative = "less"))
      test_stat2 <- c(test_stat2, ks_result$p.value)
    }
  }

  # the return object is a vector that state if the strain does belong to each phenotype
  retrun_obj <- c()
  for (i in 1:length(phenotypes)) {
    if (phenotypes[i] == observed_pheno) {
      # if the p val in the original pheno is big it means that strain is autochthonous
      # in the original pheno
      if (test_stat1[i] > 0.05) {
        retrun_obj <- c(retrun_obj, "TRUE")
      } else {
        retrun_obj <- c(retrun_obj, "FALSE")
      }
    } else {
      # if both p vals are above 0.05 then the strain belongs to the phenotype
      # if neither is bigger then it does not belog to the phenotype and
      # if only one is bigger the result is INCONCLUSIVE
      if (test_stat1[i] > 0.05 & test_stat2[i] > 0.05) {
        retrun_obj <- c(retrun_obj, "TRUE")
      } else if (test_stat1[i] < 0.05 & test_stat2[i] < 0.05) {
        retrun_obj <- c(retrun_obj, "FALSE")
      } else {
        retrun_obj <- c(retrun_obj, "INCONCLUSIVE")
      }
    }
  }

  # if test_stat1 and test_stat2 all are below 0.05 it just means that strain
  # belogs to its pheno its just not the most typical strain and thus it has its
  # prob a bit lower. Fix it here
  if (!any(!is.element(retrun_obj, "FALSE"))) {
    retrun_obj[as.character(phenotypes) == as.character(observed_pheno)] <- "TRUE not typical"
  }
  # now make the final verdict
  final <- paste(lookup_tbl[retrun_obj == "TRUE" | retrun_obj == "TRUE not typical"],
                 collapse = " or ")
  if (any(retrun_obj == "INCONCLUSIVE") & final == "") {
    final <- "INCONCLUSIVE"
  } else if (any(retrun_obj == "INCONCLUSIVE")) {
    final <- paste(final, "INCONCLUSIVE", sep = " or ")
  } else if (final == "") {
    final <- "INCONCLUSIVE"
  }
  # order the labels alphabetically
  splitted <- str_split(final, pattern = " or ")[[1]]
  final <- paste(sort(splitted), collapse = " or ")
  # return also observed phenotype
  retrun_obj <- c(retrun_obj, final, lookup_tbl[names(lookup_tbl) == as.character(observed_pheno)])
  return(retrun_obj)
}
