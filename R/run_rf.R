run_rf <- function(bin_mat_model,
                   predictors,
                   pheno,
                   sampsize,
                   mtry,
                   ntree,
                   maxnodes,
                   final_proxy_mat_sum) {

  to_dist <- function(x) {
    # function that converts coefficients to distances
    return(abs(x-1))
  }

  model <- randomForest::randomForest(data = bin_mat_model,
                                      pheno ~ .,
                                      xtest = predictors,
                                      ytest = pheno,
                                      replace = FALSE,
                                      sampsize = sampsize,
                                      oob.prox = FALSE,
                                      importance = TRUE,
                                      mtry = mtry,
                                      ntree = ntree,
                                      proximity = TRUE,
                                      maxnodes = maxnodes,
                                      keep.forest = TRUE)

  # record the proxy_mat - median
  dist_mat <- apply(model$test$proximity[1:nrow(predictors),1:nrow(predictors)], 2, to_dist)
  rownames(dist_mat) <- names(pheno)
  colnames(dist_mat) <- names(pheno)

  # record the proxy_mat - sum
  final_proxy_mat_sum <- final_proxy_mat_sum + dist_mat
  rownames(final_proxy_mat_sum) <- colnames(final_proxy_mat_sum) <- names(pheno)

  # run the prediction on the whole dataset
  predictions <- as.data.frame(predict(model, predictors, type = "prob"))

  # this will output the real label
  predictions$observed <- pheno

  return(list(dist_mat = dist_mat,
              predictions = predictions,
              final_proxy_mat_sum = final_proxy_mat_sum,
              df_importance = data.frame(model$importance)))
}
