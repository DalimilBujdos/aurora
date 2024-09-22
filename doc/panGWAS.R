## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = FALSE,
  warning = FALSE,
  message = FALSE,
  cache.lazy = FALSE
)

## ----setup--------------------------------------------------------------------
library(aurora)
library(readr)
library(ape)
library(gplots)

## -----------------------------------------------------------------------------
#save_dir <- "/path/to/your/dir/"

#if (!dir.exists(save_dir)) {dir.create(save_dir)}

data("pheno_mat_reuteri") # loads the phenotype for each strain
data("bin_mat_reuteri") # loads the pangenome matrix
data("tree_reuteri") # loads the core-genome phylogenetic tree

## ----eval = FALSE-------------------------------------------------------------
#  # I want to use all ML algorithmns
#  results <- aurora_pheno(pheno_mat = pheno_mat,
#                          bin_mat = bin_mat,
#                          tree = tree,
#                          condaenv_path = "/path/to/conda/env", # an example of a path to conda environment
#                          save_dir = save_dir)
#  
#  # I don't want to use AdaBoost and log regression
#  results <- aurora_pheno(pheno_mat = pheno_mat,
#                          bin_mat = bin_mat,
#                          tree = tree,
#                          save_dir = save_dir)
#  # lets also save the final list for GWAS analysis later
#  save(results, file = paste0(save_dir, "results_aurora.Rdata"))

## -----------------------------------------------------------------------------
data("aurora_pheno_results_reuteri") # these are only truncated results!

## -----------------------------------------------------------------------------
func <- function(x) {
  # this function will make the results look prettier
  if (is.na(x)) {
    return(NA)
  }
   if (x == 0) {
    return(0)
  }
  y <- 1
  while (x*y < 1 || x*y > 10) {
    y <- y * 10
  }
  x <- x * y
  x_rounded <- round(x, 2)
  x_final <- x_rounded / y
  return(x_final)
}

mat <- apply(results$results$results_random_forest$p_val_mat,
             c(1,2), func)
heatmap.2(mat,
          trace = "none",
          dendrogram = "none",
          cellnote=mat,
          notecol="black",
          density.info="none",
          key = FALSE)

mat <- apply(results$results$results_adaboost$p_val_mat,
             c(1,2), func)
heatmap.2(mat,
          trace = "none",
          dendrogram = "none",
          cellnote=mat,
          notecol="black",
          density.info="none",
          key = FALSE)

mat <- apply(results$results$results_log_reg$p_val_mat,
             c(1,2), func)
heatmap.2(mat,
          trace = "none",
          dendrogram = "none",
          cellnote=mat,
          notecol="black",
          density.info="none",
          key = FALSE)

mat <- apply(results$results$results_CART$p_val_mat,
             c(1,2), func)
heatmap.2(mat,
          trace = "none",
          dendrogram = "none",
          cellnote=mat,
          notecol="black",
          density.info="none",
          key = FALSE)


## -----------------------------------------------------------------------------
plot(results$results$results_random_forest$aucs, ylim = c(0,1), main = "Random Forest AUCs")
plot(results$results$results_adaboost$aucs, ylim = c(0,1), main = "AdaBoost AUCs")
plot(results$results$results_log_reg$aucs, ylim = c(0,1), main = "Log reg AUCs")

## -----------------------------------------------------------------------------
barplot(results$results$results_CART$complexity$AUC, names.arg = results$results$results_CART$complexity$classes, main = "CART AUCs", horiz = T, las = 1)

## ----eval = FALSE-------------------------------------------------------------
#  results$results$results_random_forest$plot
#  results$results$results_adaboost$plot
#  results$results$results_log_reg$plot
#  results$results$results_CART$plot

## -----------------------------------------------------------------------------
# Get all non-typical rodent strains identified by Random Forest
res_rf <- results$results$results_random_forest$auto_allo_results
res_rf_rodent <- res_rf[res_rf$observed == "rodent", ]
ids_rf_not_typ <- res_rf_rodent$strain_ids[grepl("not typical", res_rf_rodent$rodent)]

# Run GWAS analysis
res <- aurora_GWAS(bin_mat = bin_mat,
                   pheno_mat = pheno_mat,
                   tree = tree,
                   remove_strains = ids_rf_not_typ,
                   aurora_results = results,
                   mode = "consensus",
                   write_data = FALSE,
                   use_log = FALSE)

## -----------------------------------------------------------------------------
GWAS_rodent <- res$GWAS_results
# Sort the genes based on rodent standardized residuals and give them a rank
GWAS_rodent <- GWAS_rodent[order(GWAS_rodent$`rodent std_residual`, decreasing = TRUE), ]
GWAS_rodent$rank <- rank(1-GWAS_rodent$`rodent std_residual`, ties.method = "average")
# Print the rank of the known colonization factors
index <- grepl("ure", GWAS_rodent$gene_family) |
  GWAS_rodent$gene_family == "group_2460"
GWAS_rodent[index, c(1,28)]

## -----------------------------------------------------------------------------
GWAS_human <- res$GWAS_results
GWAS_human <- GWAS_human[GWAS_human$`homo sapiens precision` > 0.4, ]
GWAS_human <- GWAS_human[order(GWAS_human$`homo sapiens F1`, decreasing = TRUE), ]
# print the top 10 human colonization factors
head(GWAS_human[,c(1,19)], n = 10)

