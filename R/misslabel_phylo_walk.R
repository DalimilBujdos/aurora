#' Generates new pheno_mat with mislabelled strains using phylogenetic walk and a matrix that shows which strains were mislabelled
#'
#' @noRd
misslabel_phylo_walk <- function(pheno_mat,
                                 phenotypes,
                                 misslabel_no,
                                 tree) {

  pheno_mat$pheno <- as.character(pheno_mat$pheno)
  phenotypes <- as.character(phenotypes)
  misslabel_mat <- matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
  rownames(misslabel_mat) <- colnames(misslabel_mat) <- phenotypes
  # misslabe order is row -> column. This means that strain 1,3 was misslabeled
  # from phenotype 1 to phenotype 3.

  pheno_mat_original <- pheno_mat

  for (i in 1:length(phenotypes)) {
    # get all strains from the phenotype that were not mislabelled
    missl_strains <- as.character(misslabel_mat)
    missl_strains <- missl_strains[!is.na(missl_strains)]
    # remove strains that have already been mislabelled
    available_strains <- pheno_mat_original$ids[pheno_mat_original$pheno == phenotypes[i]]
    available_strains <- available_strains[!is.element(available_strains, missl_strains)]
    tree_class <- ape::drop.tip(tree,
                                tip = pheno_mat_original$ids[!is.element(pheno_mat_original$ids, available_strains)])
    for (j in 1:length(phenotypes)) {
      if(i == j) {next}
      misslabeled <- phylogenetic_walk(bags = vector("list", 1),
                                       tree_class,
                                       misslabeled_str = FALSE,
                                       bag_size = misslabel_no,
                                       max_per_bag = 1000,
                                       max_misslablel = FALSE,
                                       no_rounds = 1)

      tree_class <- ape::drop.tip(tree_class, tip = misslabeled$bags[[1]])
      available_strains <- available_strains[!is.element(available_strains, misslabeled$bags[[1]])]
      misslabel_mat[i, j] <- misslabeled$bags[[1]]
      index <- which(is.element(pheno_mat$ids, misslabeled$bags[[1]]))
      pheno_mat$pheno[index] <- phenotypes[j]
    }
  }
  return(list(pheno_mat_miss = pheno_mat, misslabel_mat = misslabel_mat))
}
