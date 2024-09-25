train_gnb <- function(phenotype, phenotypes, df_threshold) {
  # this function trains gaussian naive bayes model on data from CART classificator

  # get the training data for the model
  pheno_search <- paste0(phenotype, "to", phenotypes)
  pheno_search <- c(phenotype, pheno_search[-as.numeric(phenotype)])
  train_data <- df_threshold[is.element(df_threshold$pheno, pheno_search), ]
  train_pheno <- as.factor(as.character(train_data$pheno))# this rerets the original levels
  train_data <- as.matrix(train_data[,-ncol(train_data)])
  train_data <- apply(train_data, 2, as.numeric)
  prior_prob <- rep(1/nlevels(train_pheno), nlevels(train_pheno))
  # fit the model
  model <- gaussian_naive_bayes(train_data, train_pheno, prior = prior_prob)
  return(model)
}
