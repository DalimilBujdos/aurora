run_rf_miss <- function(bin_mat_miss,
                        pheno_mat_miss,
                        predictors_test,
                        pheno_mat_test,
                        sampsize,
                        mtry,
                        ntree,
                        maxnodes) {

  bin_mat_miss$pheno <- as.factor(pheno_mat_miss$pheno)

  # fit random forest
  model_1 <- randomForest::randomForest(data = bin_mat_miss,
                                        pheno ~ .,
                                        replace = FALSE, # sample without replacement (already did the sampling)
                                        sampsize = sampsize, # sampsize has to be smaller -1 because it otherwise runs into error
                                        mtry = mtry,
                                        ntree = ntree,
                                        maxnodes = maxnodes)

  # nodesize = ceiling(node_multip * nrow(predictors_miss))
  predictions <- as.data.frame(stats::predict(model_1, predictors_test, type = "prob"))
  # this will output the real label into a new column
  predictions$observed <- pheno_mat_test
  return(predictions)
}
