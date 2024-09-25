predict_auto_allo <- function(strain_prob, model, phenotypes) {
  # functio that predicts the auto/allo label based on probabilities from
  # outlier calculation
  prob <- as.character(predict(model, as.matrix(strain_prob)))
  final_probs <- setNames(rep(0, length(phenotypes)), phenotypes)
  tbl <- as.data.frame(table(substr(prob, nchar(prob), nchar(prob))))
  final_probs[match(tbl$Var1, names(final_probs))] <- tbl$Freq
  # normalise the result
  final_probs <- (final_probs/sum(final_probs))*100
  return(final_probs)
}
