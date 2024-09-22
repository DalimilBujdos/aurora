#' Postprocess data from Random Forest and AdaBoost in Outlier calculation phase
#'
#' @noRd
postprocess_fit <- function(predictions,
                            probs,
                            df_importance = NA,
                            phenotypes,
                            rf = FALSE,
                            ada = FALSE) {
  # calculate the area under the curve
  roc <- pROC::multiclass.roc(as.factor(predictions$observed),
                              predictions[,1:length(phenotypes)],
                              quiet = TRUE)
  auc <- roc$auc

  # store all the probabilities into one list
  for (i in 1:nrow(predictions)) {
    index <- which(names(probs) == rownames(predictions)[i])
    new_probs <- as.data.frame(predictions[i,1:length(phenotypes)])
    probs[[index]] <- rbind(probs[[index]], new_probs)
  }

  # store all the feature importances into one list
  if (rf == TRUE) {
    return(list(auc = auc, probs = probs,
                importance_gini = df_importance[,colnames(df_importance) == "MeanDecreaseGini"],
                importance_accuracy = df_importance[,colnames(df_importance) == "MeanDecreaseAccuracy"]))
  }

  if (ada == TRUE) {
    return(list(auc = auc, probs = probs, importance_gini = df_importance[,2]))
  }

  return(list(auc = auc, probs = probs))
}
