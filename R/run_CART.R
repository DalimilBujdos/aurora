#' Run CART models in Outlier calculation phase
#'
#' @noRd
run_CART <- function(phenotype,
                     bin_mat,
                     pheno_mat,
                     dist_mat_or_tree,
                     bagging,
                     bag_size,
                     max_per_bag) {
  # runs CART classificator
  phenotype1 <- substr(phenotype, 1,1)
  phenotype2 <- substr(phenotype, 2,2)
  # get the bin_mat and pheno_mat for this run
  pheno_mat <- pheno_mat[pheno_mat$pheno == phenotype1 |
                           pheno_mat$pheno == phenotype2, ]
  pheno_mat$pheno <- as.character(pheno_mat$pheno)

  proxy_mat <- matrix(NA, nrow = nrow(pheno_mat), ncol = nrow(pheno_mat))
  colnames(proxy_mat) <- rownames(proxy_mat) <- pheno_mat$ids

  bin_mat <- bin_mat[is.element(rownames(bin_mat), pheno_mat$ids),
                     colnames(bin_mat) != "ids"]

  if (bagging == "random_walk") {
    dist_mat_or_tree <- dist_mat_or_tree[is.element(rownames(dist_mat_or_tree), pheno_mat$ids),
                         is.element(colnames(dist_mat_or_tree), pheno_mat$ids)]
  }

  if (bagging == "phylogenetic_walk") {
    other_strains <- dist_mat_or_tree$tip.label[!is.element(dist_mat_or_tree$tip.label, pheno_mat$ids)]
    dist_mat_or_tree <- ape::drop.tip(dist_mat_or_tree, tip = other_strains)
  }

  result_get_bags <- get_bags(pheno_mat,
                              dist_mat_or_tree,
                              bagging,
                              misslabel_mat = FALSE,
                              bag_size,
                              no_rf = 1,
                              max_per_bag)

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
    names(model$where) <- rownames(bin_mat_train)
  }

  predictions <- as.data.frame(predict(model, bin_mat, type = "prob"))
  rownames(predictions) <- rownames(bin_mat)
  # this will output the real label into a new column
  predictions$ids <- pheno_mat$ids

  # now summarise the variable importance
  importance <- 100/sum(model$variable.importance) * model$variable.importance
  index <- as.numeric(substr(names(importance), 2, nchar(names(importance))))
  df_importance <- data.frame(importance = rep(0,ncol(bin_mat)))
  df_importance$importance[index] <- importance

  # calculate how many splits in general was necessary in the model
  no_splits <- model$cptable[nrow(model$cptable),2]
  # calculate the AUCS
  predictions$observed <- pheno_mat$pheno
  curr_auc <- pROC::roc(as.factor(predictions$observed),
                        predictions[,1],
                        quiet = TRUE)

  # sometimes the model is pruned just to its root. If that is the case than just return
  # an empty proxy mat
  if (nrow(model$frame) == 1) {
    list_return <- list(probs = predictions,
                        importance = df_importance,
                        splits = no_splits,
                        curr_auc = curr_auc$auc,
                        heatmat = proxy_mat)
    return(list_return)
  }

  # for some reason some nodes in model$frame are shifted by one. Fix it here
  old <- unique(model$where)
  old <- old[order(old)]
  convers_tbl <- data.frame(new = rownames(model$frame)[model$frame$var == "<leaf>"],
                            old = old)
  new <- convers_tbl$new[match(model$where, convers_tbl$old)]
  names(new) <- names(model$where)
  model$where <- new

  # compute the distance matrix
  # remove the .number
  names(model$where) <- sub("\\.(\\d+)$", "", names(model$where))

  # now create a df that will hold all the node labels and the corresponding distances
  # path.rpart outputs things into command line
  invisible(capture.output(paths <- rpart::path.rpart(model, nodes = unique(model$where))))
  pairs <- as.data.frame(expand.grid(unique(model$where), unique(model$where)))
  pairs <- as.data.frame(sapply(pairs, as.character))
  pairs$dist <- rep(NA, nrow(pairs))
  for (i in 1:nrow(pairs)) {
    if (pairs$Var1[i] == pairs$Var2[i]) {
      pairs$dist[i] <- 0
      next
    }
    path1 <- paths[[which(names(paths) == pairs$Var1[i])]]
    path2 <- paths[[which(names(paths) == pairs$Var2[i])]]

    tip1 <- strsplit(path1[length(path1)], "=")[[1]][1]
    tip2 <- strsplit(path2[length(path2)], "=")[[1]][1]
    if (tip1 == tip2) {
      pairs$dist[i] <- 1
      next
    }

    path1 <- gsub("=miss", "", path1)
    path1 <- gsub("=hit", "", path1)
    path1 <- path1[path1 != "root"]
    path2 <- gsub("=miss", "", path2)
    path2 <- gsub("=hit", "", path2)
    path2 <- path2[path2 != "root"]

    not_in_path1 <- setdiff(path1, path2)
    not_in_path2 <- setdiff(path2, path1)
    total_count <- length(not_in_path1) + length(not_in_path2)
    pairs$dist[i] <- 1+total_count
  }

  # get only the unique vals
  unique_vals <- unique(names(model$where))
  indexes <- sapply(unique_vals, function(x) which(names(model$where) == x)[1])
  strains_nodes <- model$where[indexes]

  for (i in 1:nrow(proxy_mat)) {
    strain <- rownames(proxy_mat)[i]
    if (!any(names(strains_nodes) == strain)) {
      next
    }
    node_strain <- strains_nodes[names(strains_nodes) == strain]
    pairs_curr <- pairs[pairs$Var1 == node_strain, ]
    dists <- rep(NA, ncol(proxy_mat))
    for (j in 1:ncol(proxy_mat)) {
      index <- which(names(strains_nodes) == colnames(proxy_mat)[j])
      if (length(index) > 0) {
        dists[j] <- pairs_curr$dist[pairs_curr$Var1 == node_strain & pairs_curr$Var2 == strains_nodes[index]]
      }
    }
    proxy_mat[i,] <- dists
  }


  list_return <- list(probs = predictions,
                      importance = df_importance,
                      splits = no_splits,
                      curr_auc = curr_auc$auc,
                      heatmat = proxy_mat)
  return(list_return)
}
