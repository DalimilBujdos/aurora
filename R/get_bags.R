#' Functions that uses random walk or phylogenetic walk to obtain training datasets
#'
#' @param pheno_mat Data frame that contains unique indexes in the first column
#'                  and the phenotype classes in the second column. The unique indexes
#'                  should contain only letters, numbers and special signs "_", ".".
#' @param dist_mat_or_tree  Provide either distance matrix that will contain phylogenetic
#'                          distances or a phylogenetic tree as a \code{phylo} object.
#' @param bagging The algorithm used to bootstrap the original dataset. Choose either
#'                "phylogenetic_walk" or "random_walk". Default: "phylogenetic_walk"
#' @param misslabel_mat Matrix that indicates which strains were mislabeled. This is argument
#'                      is only used internally in \code{aurora_pheno}. Default: FALSE
#' @param bag_size  The size of the bag for each class. Default: NA. If NA than the
#'                  bag_size is calculated as 5* the number of strains in the class
#'                  with the fewest strains. Provide the size as a number for each class
#'                  i.e., c(50, 50, 50) for a phenotype with three classes.
#' @param no_rf Number of bags that the function returns. Default: 100
#' @param max_per_bag Maximum number of times a strain can be repeated in the bag. Default: NA
#'                    If NA than the max_per_bag is calculated so that none of the strains
#'                    exceeds 20% of the bag of each class.
#'
#' @return The output of this function is a list with strain indexes selected by phylogenetic_walk" or "random_walk".
#'         The list is as long as the \code{no_rf} argument specifies. Additionally the output contains count that shows
#'         how many times was each strain selected in the bag.
#' @export
#'
#' @examples
#' \dontrun{
#'   data(tree_reuteri)
#'   data(pheno_mat_reuteri)
#'
#'   get_bags(pheno_mat = pheno_mat,
#'            tree = tree,
#'            no_rf = 100) # generates 100 bags using phylogenetic_walk
#' }
get_bags <- function(pheno_mat,
                     dist_mat_or_tree,
                     bagging = "phylogenetic_walk",
                     misslabel_mat = FALSE,
                     bag_size = NA,
                     no_rf = 1,
                     max_per_bag = NA) {
  # calculate how big the bags should be
  # it is the size of the smallest class*10
  if (any(is.na(bag_size))) {
    table_pheno_mat <- as.data.frame(table(pheno_mat$pheno))
    bag_size <- rep(min(table_pheno_mat$Freq)*5, nrow(table_pheno_mat))
  } else {
    table_pheno_mat <- as.data.frame(table(pheno_mat$pheno))
  }

  bags <- vector(mode = "list", length = no_rf)
  count_strains_class <- vector(mode = "list", length = nrow(table_pheno_mat))
  # get the maximum number of genome copies in the bag
  if (is.na(max_per_bag)) {
    max_per_bag <- round(bag_size[which.min(bag_size)]*0.2) # should not exceed 20% of the bag
  }

  for (i in 1:length(table(pheno_mat$pheno))) {
    curr_phenotype <- table_pheno_mat$Var1[i]
    if (bagging == "random_walk") {
      # get dist_mat for starins of the class
      index <- is.element(rownames(dist_mat_or_tree),
                          pheno_mat$ids[pheno_mat$pheno == table_pheno_mat$Var1[i]])
      dist_mat_class <- dist_mat_or_tree[index, index]
    }

    if (bagging == "phylogenetic_walk") {
      # get the subtree for starins of the class
      tree_class <- ape::drop.tip(dist_mat_or_tree, tip = pheno_mat$ids[pheno_mat$pheno != curr_phenotype])
      if (!ape::is.rooted(tree_class)) {tree_class <- phangorn::midpoint(tree_class)}

    }

    no_strains_class <- ifelse(bagging == "phylogenetic_walk", length(tree_class$tip.label), nrow(dist_mat_class))

    if (is.matrix(misslabel_mat)) {
      misslabeled_str <- as.character(misslabel_mat[,colnames(misslabel_mat) == curr_phenotype])
      misslabeled_str <- misslabeled_str[!is.na(misslabeled_str)]
      # calculate the maximum number of repetitions a misslabelled strain can have
      # this prevents the misslabelled strains to overtake the whole bag
      max_misslablel <- ceiling(bag_size[i]/no_strains_class) + (round(bag_size[i]*0.1))
      if (max_misslablel == max_per_bag) {max_misslablel <- FALSE} # this is here in case max_misslablel and max_per_bag are the same values
    } else if (is.character(misslabel_mat)) {
      # this is here for CART
      misslabeled_str <- misslabel_mat[is.element(names(misslabel_mat), table_pheno_mat$Var1[i])]
      max_misslablel <- ceiling(bag_size[i]/no_strains_class) + (round(bag_size[i]*0.1))
      if (max_misslablel == max_per_bag) {max_misslablel <- FALSE} # this is here in case max_misslablel and max_per_bag are the same values
    } else {
      misslabeled_str <- FALSE
      max_misslablel <- FALSE
    }
    if (bagging == "random_walk") {
      result_bag <- random_walk(bags,
                                dist_mat_class,
                                misslabeled_str,
                                bag_size[i],
                                max_per_bag,
                                max_misslablel,
                                no_rf)
    } else if (bagging == "phylogenetic_walk") {
      result_bag <- phylogenetic_walk(bags,
                                      tree_class,
                                      misslabeled_str,
                                      bag_size[i],
                                      max_per_bag,
                                      max_misslablel,
                                      no_rf)
    }
    bags <- result_bag$bags
    count_strains_class[[i]] <- result_bag$count
  }
  return(list(bags = bags, count = count_strains_class))
}
