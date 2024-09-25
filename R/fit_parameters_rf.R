#' Fits hyperparameters for Random Forest models
#'
#' @noRd
fit_parameters_rf <- function(bags,
                              predictors,
                              pheno_mat,
                              phenotypes,
                              repeats) {
  # we could use caret library for grid search here. The only problem is that caret
  # optimises only for mtry hyperparameter which is not enough

  # instead of cross validation (cv) that is normally used in parameter fine tuning
  # we are using validation on multiple bags

  # create a grid that will hold all the parameters to test

  # the 0.3* is here cuz only 30% of the bag is used for training but if the bag size is small (<70) use 60%
  if (length(bags[[1]]) >= 70) {
    smp_size <- c(0.2*floor(0.3*length(bags[[1]])),
                  0.4*floor(0.3*length(bags[[1]])),
                  0.6*floor(0.3*length(bags[[1]])),
                  0.8*floor(0.3*length(bags[[1]])))
    fraction <- 0.3
  } else {
    smp_size <- c(0.2*floor(0.6*length(bags[[1]])),
                  0.4*floor(0.6*length(bags[[1]])),
                  0.6*floor(0.6*length(bags[[1]])),
                  0.8*floor(0.6*length(bags[[1]])))
    fraction <- 0.6
  }

  mtry_arg <- c(10, 50, 200, 500, 1000, 2000)
  mtry_arg <- mtry_arg[mtry_arg < ncol(predictors)]

  ntree_arg <- c(10, 50, 100, 500, 1000)

  maxnodes_arg <- c(4, 8 ,12)
  if (smp_size[1] < 4) {
    maxnodes_arg <- 2
  } else {
    maxnodes_arg <- maxnodes_arg[maxnodes_arg < smp_size[1]]
  }


  Grid <-  expand.grid(sampsize = smp_size,
                       mtry = mtry_arg,
                       ntree = ntree_arg,
                       maxnodes = maxnodes_arg)

  run_rf_param <- function(x, train_data, train_data_pheno, test_data, test_data_pheno, phenotypes) {
    # runs random forest and tries to find the best parameters of the model
    x <- as.numeric(x)
    model <- randomForest::randomForest(x = train_data,
                                        y = train_data_pheno,
                                        sampsize = x[1],
                                        mtry = x[2],
                                        ntree = x[3],
                                        maxnodes = x[4],
                                        replace = FALSE,
                                        oob.prox = FALSE,
                                        importance = FALSE,
                                        proximity = FALSE)

    predictions <- as.data.frame(stats::predict(model,
                                                test_data,
                                                type = "prob"))

    predictions$observed <- test_data_pheno

    # calculate accuracy of a model.
    predict_best <- function(x, phenotypes) {
      return(names(x)[which.max(x[1:length(phenotypes)])])
    }

    predictions$predicted <- apply(predictions, 1, predict_best, phenotypes = phenotypes)
    correct_pred <- ifelse(predictions$predicted == predictions$observed, 1, 0)
    accuracy <- sum(correct_pred)/length(correct_pred)

    return(accuracy)
  }

  df_final_accuracy <- data.frame(index = 1:nrow(Grid))

  for (i in 1:repeats) {
    #  construct pheno_mat for this current run
    index_select <- c()
    for (j in 1:length(bags[[i]])) {
      index_row <- which(pheno_mat$ids == bags[[i]][j])
      index_select <- c(index_select, index_row)
    }
    pheno_mat_model <- pheno_mat[index_select, ]

    # construct bin_mat for this current run
    index_select <- c()
    for (j in 1:nrow(pheno_mat_model)) {
      the_row <- which(rownames(predictors) == pheno_mat_model$ids[j])
      index_select <- c(index_select, the_row[1])
    }
    predictors_train <- predictors[index_select, ]
    predictors_train <- predictors_train[,colnames(predictors_train) != "ids"]
    predictors_train$pheno <- pheno_mat_model$pheno

    predictors_learn <- data.frame()
    predictors_test <- data.frame()
    # split data into test and training set
    for (j in 1:length(phenotypes)) {
      subset <- predictors_train[predictors_train$pheno == phenotypes[j], ]
      # 30% or 60% of the whole dataset is for learning and the rest for testing

      smpl_size <- round(fraction * nrow(subset))
      index <- base::sample(1:nrow(subset), as.integer(smpl_size))
      predictors_learn <- rbind(predictors_learn, subset[index, ])
      predictors_test <- rbind(predictors_test, subset[-index, ])
    }

    accuracy_result <- apply(Grid,
                             1,
                             run_rf_param,
                             train_data = predictors_learn[,!(colnames(predictors_learn) == "pheno")],
                             train_data_pheno = as.factor(predictors_learn$pheno),
                             test_data = predictors_test[,!(colnames(predictors_test) == "pheno")],
                             test_data_pheno = as.factor(predictors_test$pheno),
                             phenotype = phenotypes)

    df_final_accuracy <- cbind(df_final_accuracy, as.data.frame(accuracy_result))

    cat("Random forest grid search finished", i, "out of", repeats, "\n" )
  }
  df_final_accuracy <- df_final_accuracy[,-1]
  medians_acc <- apply(df_final_accuracy, 1, stats::median)
  return(Grid[which.max(medians_acc), ])


}
