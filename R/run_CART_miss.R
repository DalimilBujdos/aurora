#' Run CART models in Threshold calculation phase
#'
#' @noRd
run_CART_miss <- function(phenotype,
                          bin_mat,
                          pheno_mat,
                          dist_mat_or_tree,
                          bagging,
                          misslabel_no,
                          bag_size,
                          max_per_bag) {
  # runs CART classificator with misslabelled strains
  phenotype1 <- substr(phenotype, 1,1)
  phenotype2 <- substr(phenotype, 2,2)
  # get the bin_mat and pheno_mat for this run
  pheno_mat <- pheno_mat[pheno_mat$pheno == phenotype1 | pheno_mat$pheno == phenotype2, ]

  bin_mat <- bin_mat[is.element(rownames(bin_mat), pheno_mat$ids),
                     colnames(bin_mat) != "ids"]

  if (bagging == "random_walk") {
    # misslabel strains and build a bag
    # if bagging is a random walk then the sampling of the mislabelled strain is random
    oneTtwo <- base::sample(pheno_mat$ids[pheno_mat$pheno == phenotype1], misslabel_no)
    twoTone <- base::sample(pheno_mat$ids[pheno_mat$pheno == phenotype2], misslabel_no)
    pheno_mat$pheno[is.element(pheno_mat$ids, oneTtwo)] <- phenotype2
    pheno_mat$pheno[is.element(pheno_mat$ids, twoTone)] <- phenotype1
    pheno_mat$pheno <- as.character(pheno_mat$pheno)

    dist_mat_or_tree <- dist_mat_or_tree[is.element(rownames(dist_mat_or_tree), pheno_mat$ids),
                         is.element(colnames(dist_mat_or_tree), pheno_mat$ids)]
  }

  if (bagging == "phylogenetic_walk") {
    # misslabel strains and build a bag
    # if bagging is a phylo walk then the sampling of the mislabelled strain is based on phylo walk
    pheno_strains <- pheno_mat$ids[pheno_mat$pheno == phenotype1]
    tree_class <- ape::drop.tip(dist_mat_or_tree,
                                tip = dist_mat_or_tree$tip.label[!is.element(dist_mat_or_tree$tip.label, pheno_strains)])
    oneTtwo <- phylogenetic_walk(bags = vector("list", 1),
                                 tree_class,
                                 misslabeled_str = FALSE,
                                 bag_size = misslabel_no,
                                 max_per_bag = 1000,
                                 max_misslablel = FALSE,
                                 no_rounds = 1)
    oneTtwo <- unique(oneTtwo$bags[[1]])

    pheno_strains <- pheno_mat$ids[pheno_mat$pheno == phenotype2]
    tree_class <- ape::drop.tip(dist_mat_or_tree,
                                tip = dist_mat_or_tree$tip.label[!is.element(dist_mat_or_tree$tip.label, pheno_strains)])
    twoTone <- phylogenetic_walk(bags = vector("list", 1),
                                 tree_class,
                                 misslabeled_str = FALSE,
                                 bag_size = misslabel_no,
                                 max_per_bag = 1000,
                                 max_misslablel = FALSE,
                                 no_rounds = 1)
    twoTone <- unique(twoTone$bags[[1]])
    pheno_mat$pheno[is.element(pheno_mat$ids, oneTtwo)] <- phenotype2
    pheno_mat$pheno[is.element(pheno_mat$ids, twoTone)] <- phenotype1
    pheno_mat$pheno <- as.character(pheno_mat$pheno)

    # get tree for phylo walk
    other_strains <- dist_mat_or_tree$tip.label[!is.element(dist_mat_or_tree$tip.label, pheno_mat$ids)]
    dist_mat_or_tree <- ape::drop.tip(dist_mat_or_tree, tip = other_strains)

  }

  misslabel_mat <- setNames(c(oneTtwo,twoTone),
                            as.character(c(rep(phenotype2, length(oneTtwo)),
                                           rep(phenotype1, length(twoTone)))))
  result_get_bags <- get_bags(pheno_mat,
                              dist_mat_or_tree,
                              bagging,
                              misslabel_mat,
                              bag_size,
                              no_rf = 1,
                              max_per_bag)
  # construct training bin_mat and pheno_mat

  curr_bag_order <- result_get_bags$bags[[1]]
  # construct bin_mat and pheno_mat for this current run
  index_select <- c()
  for (i in 1:length(curr_bag_order)) {
    the_row <- which(rownames(bin_mat) == curr_bag_order[i])
    index_select <- c(index_select, the_row)
  }
  bin_mat_train <- bin_mat[index_select, ]
  pheno_mat_train <- pheno_mat[index_select, ]
  # fit and prune the model then run predict on the test data
  bin_mat_train$pheno <- pheno_mat_train$pheno
  # the problem is that the decision tree overfits due to the fact that phylo_walk
  # results in many repetitions of the mislabelled genome. A possible solution is to
  # set the minbucket argument to max_misslablel+1. this will force the algorithm to
  # not include the mislabelled strains into one leaf node
  no_strains_class <- as.data.frame(table(pheno_mat$pheno))
  no_strains_class <- no_strains_class$Freq[which.min(no_strains_class$Freq)]
  bagsize_local <- length(curr_bag_order)/2
  minbucket_val <- ceiling(bagsize_local/no_strains_class) + (round(bagsize_local *0.1))

  model <- rpart::rpart(pheno ~ .,
                        data = bin_mat_train,
                        method = "class",
                        control = rpart::rpart.control(cp = 0.01, minsplit = 2, minbucket = minbucket_val+1, xval = 4))
  # get complexity parameter (cp) with the one standard deviation rule and prune the tree
  index <- which.min(model$cptable[,4])
  threshold <- model$cptable[index,4] + model$cptable[index,5]
  if (threshold != 0) {
    # threshold is 0 when the tree has only one split
    new_cp <- model$cptable[which(model$cptable[,4] < threshold)[1], 1]
    model <- rpart::prune(model, cp = new_cp)
  }

  predictions <- as.data.frame(stats::predict(model, bin_mat, type = "prob"))
  rownames(predictions) <- rownames(bin_mat)
  # this will output the real label into a new column
  predictions$observed <- pheno_mat$pheno
  # the result of the function is a dataframe where the misslabelled strains have
  # label XtoY
  predictions$observed[is.element(rownames(predictions), oneTtwo)] <- paste0(phenotype1, "to", phenotype2)
  predictions$observed[is.element(rownames(predictions), twoTone)] <- paste0(phenotype2, "to", phenotype1)
  return(predictions)
}
