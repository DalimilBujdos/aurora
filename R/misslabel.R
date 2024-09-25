#' Generates new pheno_mat with randomly mislabelled strains and a matrix that shows which strains were mislabelled
#'
#' @noRd
misslabel <- function(pheno_mat,
                      phenotypes,
                      misslabel_no = 1) {
  pheno_mat$pheno <- as.character(pheno_mat$pheno)
  phenotypes <- as.character(phenotypes)
  misslabel_mat <- matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
  rownames(misslabel_mat) <- colnames(misslabel_mat) <- phenotypes
  # misslabe order is row -> column. This means that strain 1,3 was misslabeled
  # from phenotype 1 to phenotype 3.

  pheno_mat_original <- pheno_mat

  for (i in 1:length(phenotypes)) {
    original_strains <- unique(pheno_mat_original$ids[pheno_mat_original$pheno
                                                           == phenotypes[i]])
    available_strains <- pheno_mat$ids[is.element(pheno_mat$ids, original_strains)]
    for (j in 1:length(phenotypes)) {
      if(i == j) {next}
      misslabeled <- base::sample(available_strains, misslabel_no)
      available_strains <- available_strains[!is.element(available_strains, misslabeled)]
      misslabel_mat[i, j] <- misslabeled
      index <- which(is.element(pheno_mat$ids, misslabeled))
      pheno_mat$pheno[index] <- phenotypes[j]
    }
  }


  return(list(pheno_mat_miss = pheno_mat, misslabel_mat = misslabel_mat))
}
