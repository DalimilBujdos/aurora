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
library(phangorn)

## -----------------------------------------------------------------------------
data("pheno_mat_reuteri")
data("bin_mat_reuteri")

# Bootstrap the dataset to make it larger
index <- sample(1:nrow(pheno_mat), size = 2000, replace = TRUE)
pheno_mat <- pheno_mat[index, ]
bin_mat_tmp <- bin_mat[,-1:-14]
bin_mat <- cbind(bin_mat[,1:14], bin_mat_tmp[, match(pheno_mat$ids, colnames(bin_mat_tmp))])
# give the strains generic names because the indexes must be unique
pheno_mat$ids <- paste0("strain", 1:nrow(pheno_mat))
colnames(bin_mat)[-1:-14] <- paste0("strain", 1:nrow(pheno_mat))
# Create a random tree
tree <- rtree(2000, tip.label = pheno_mat$ids)

## -----------------------------------------------------------------------------
# run GWAS analysis without results from aurora_pheno()
res <- aurora_GWAS(bin_mat = bin_mat,
                   pheno_mat = pheno_mat,
                   tree = midpoint(tree),
                   write_data = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  # Plot the phyogenetic tree first to see if the dataset is clonal
#  plot(tree, type = "fan", show.tip.label = FALSE)
#  # Thy phylogenetic tree does not appear to contain large clonal lineages that we could collapse into a single strain. Let's thus look at other ways to reduce the computational time
#  
#  # If one class in the dataset is not dominant them we could reduce the argument bag_size
#  table(pheno_mat$pheno)
#  # While rodent and poultry isolates are clearly more abundant the numbers are roughly equal and we can reduce bag_size
#  table_tmp <- table(pheno_mat$pheno)
#  bag_size <- rep(table_tmp[which.min(table_tmp)]*2, length(table_tmp))
#  
#  # Another thing that will drastically reduce the computational time is lowering the number of classes. The results then has to be interpreted accordingly.
#  pheno_mat$pheno[pheno_mat$pheno != "rodent"] <- "other" # make only two classes: rodent and other
#  
#  aurora_pheno(pheno_mat = pheno_mat,
#               bin_mat = bin_mat,
#               tree = tree,
#               bag_size = unname(bag_size),
#               bagging = "random_walk", # use random walk which is faster than phylogenetic walk
#               repeats = 5, # This will lower the time needed to fit hyperparameters to Random Forest
#               ovr_log_reg = FALSE, # Do not use log regression
#               adaboost = FALSE , # Do not use AdaBoost
#               no_rounds = 50, # From literature we know that reuteri colonization factors have high effect size. Therefore we can reduce this argument
#               write_data = FALSE)

