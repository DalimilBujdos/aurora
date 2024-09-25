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
library(mlr)
library(dplyr)
library(randomForest)
library(xgboost)
library(ggplot2)
library(proxy)
library(FactoMineR)

## -----------------------------------------------------------------------------
data("pheno_mat_reuteri") # loads the phenotype for each strain
data("bin_mat_reuteri") # loads the pangenome matrix
data("tree_reuteri") # loads the core-genome phylogenetic tree

# convert the pangenome matrix to binary matrix
bin_mat <- as.data.frame(bin_mat)
fut_row_names <- bin_mat$Gene
bin_mat <- bin_mat[, -1:-14]
bin_mat <- bin_mat %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate_all(~replace(., .!=0, 1))
bin_mat <- sapply(bin_mat, as.numeric)
bin_mat <- t(bin_mat)
colnames(bin_mat) <- fut_row_names

bin_mat <- bin_mat[match(pheno_mat$ids, rownames(bin_mat)), ]
no_gf <- nrow(bin_mat)
cutoff <- no_gf * (95/100) # remove abundat features
select_ID_up <- which(colSums(bin_mat) < cutoff)

cutoff <- no_gf * (5/100) # remove rare features
non_select_ID_low <- which(colSums(bin_mat) > cutoff)
bin_mat <- bin_mat[, intersect(select_ID_up, non_select_ID_low)]

## -----------------------------------------------------------------------------
# Get the train dataset
bag <- get_bags(pheno_mat,
                dist_mat_or_tree = tree, 
                bag_size = c(200, 200, 200, 200, 200) # number of strains in each class
                )
# The number of strains in each class should be balanced for Random Forest
# Get bags also outputs $count which shows the number of strain repetitions

# Construct new pheno_mat and bin_mat
boot_pheno_mat <- data.frame(ids = bag$bags[[1]],
                             pheno = pheno_mat$pheno[match(bag$bags[[1]], pheno_mat$ids)])
boot_bin_mat<- as.data.frame(bin_mat[match(bag$bags[[1]], rownames(bin_mat)), ])
boot_bin_mat$host <- boot_pheno_mat$pheno

# gene names from panaroo are not valid column names
gene_names <- colnames(boot_bin_mat)[1:ncol(boot_bin_mat)-1]
colnames(boot_bin_mat)[1:ncol(boot_bin_mat)-1] <- paste0(rep("gene", ncol(boot_bin_mat)-1), 1:(ncol(boot_bin_mat)-1))

## -----------------------------------------------------------------------------
task <- makeClassifTask(data = boot_bin_mat, target = "host")
# define hyperparameter search space
param_set <- makeParamSet(
  makeIntegerParam("sampsize", lower = 50, upper = 200),
  makeIntegerParam("ntree", lower = 50, upper = 500),
  makeIntegerParam("mtry", lower = 10, upper = 1000), 
  makeIntegerParam("maxnodes", lower = 4, upper = 20)
  )

ctrl <- makeTuneControlRandom(maxit = 10) # set search method as random

result <- tuneParams(
  makeLearner("classif.randomForest"),
  task,
  resampling = makeResampleDesc("CV", iters = 4),  # 4-fold cross-validation
  measures = list(acc),  # Use accuracy as the performance measure
  par.set = param_set,
  control = ctrl
)

best_params <- result$x

## ----eval = FALSE-------------------------------------------------------------
#  data("pheno_mat_reuteri")
#  data("bin_mat_reuteri")
#  data("tree_reuteri")
#  results_aurora <- aurora_pheno(pheno_mat = pheno_mat,
#                                 bin_mat = bin_mat,
#                                 type_bin_mat = "panaroo",
#                                 tree = tree,
#                                 fit_parameters = FALSE,
#                                 sampsize = best_params$sampsize,
#                                 mtry = best_params$mtry,
#                                 ntree = best_params$ntree,
#                                 maxnodes = best_params$maxnodes,
#                                 ovr_log_reg = FALSE,
#                                 adaboost = FALSE,
#                                 CART = FALSE,
#                                 write_data = FALSE)
#  

## -----------------------------------------------------------------------------
# Get the train dataset
bag <- get_bags(pheno_mat,
                dist_mat_or_tree = tree, 
                bag_size = c(20, 20, 20, 20, 20), # number of strains in each class
                no_rf = 10, # specify the number of datasets
                bagging = "phylogenetic_walk" # specify the bagging algorithm
                )


## -----------------------------------------------------------------------------
# make a list that will hold all the feature importances
importances <- vector(mode = "list", 10)
# make a list that will hold all the classification probabilities
probs <- vector(mode = "list", 10)
# gene names from panaroo are not valid column names
gene_names <- colnames(bin_mat)[1:ncol(bin_mat)]
bin_mat_adjust <- bin_mat
colnames(bin_mat_adjust)[1:ncol(bin_mat_adjust)] <- paste0(rep("gene", ncol(bin_mat_adjust)), 1:(ncol(bin_mat_adjust)))
names(gene_names) <- colnames(bin_mat_adjust)

for (i in 1:10) {
  # Construct new pheno_mat and bin_mat
  boot_pheno_mat <- data.frame(ids = bag$bags[[i]],
                             pheno = pheno_mat$pheno[match(bag$bags[[i]], pheno_mat$ids)])
  boot_pheno_mat$pheno <- ifelse(boot_pheno_mat$pheno == "rodent", 1, 0) # 1 = rodent, 0 = other host
  boot_bin_mat<- as.data.frame(bin_mat_adjust[match(bag$bags[[i]], rownames(bin_mat_adjust)), ])
  
  # train XGBoost model
  model <- xgboost(data = as.matrix(boot_bin_mat),
                   label = boot_pheno_mat$pheno,
                   nrounds = 10,
                   objective = "binary:logistic",
                   verbose = 0)
  importances[[i]] <- as.matrix(xgb.importance(model = model)[,c(1,2)])
  probs[[i]] <- predict(model, bin_mat_adjust)
}

## -----------------------------------------------------------------------------
# analyze the classification probabilities
medians <- apply(do.call(cbind, probs), 1, median)
df_probs <- data.frame(ids = rownames(bin_mat_adjust),
                       probs_med = medians,
                       host = pheno_mat$pheno)

# analyze the feature importances
df_importances <- as.data.frame(do.call(rbind, importances))
df_importances$Feature <- gene_names[match(df_importances$Feature, names(gene_names))]
features <- as.data.frame(table(df_importances$Feature))
features <- features[features$Freq > 1, ]
features$importance <- rep(NA, nrow(features))
for (i in 1:nrow(features)) {
  features$importance[i] <- median(as.numeric(df_importances$Gain[df_importances$Feature == features$Var1[i]]))
}
features <- features[order(features$Freq, features$importance, decreasing = TRUE), ]
colnames(features) <- c("Gene", "Frequence", "Importance")
# print the top 15 features along with their frequences in the 10 models and the median of their importances
knitr::kable(features[1:15,], format = "html", table.attr = "class='table table-striped'")
head(features, 15)

## -----------------------------------------------------------------------------
df_probs <- df_probs[order(df_probs$probs_med, decreasing = TRUE), ]

ggplot(df_probs, aes(x = 1:nrow(df_probs), y = probs_med, color = host)) +
  geom_point() +
  labs(title = "XGBoost results",
       x = "analysed strains",
       y = "Classification probability",
       color = "Host")


## ----echo=FALSE, figures-side, fig.show="hold", out.width="50%"---------------
# plot the original dataset
jaccard_dist <- 1 - proxy::dist(as.matrix(bin_mat), method = "Jaccard")
pca_result <- FactoMineR::PCA(as.matrix(jaccard_dist), graph = FALSE)
pca_coordinates <- as.data.frame(pca_result$ind$coord)
pca_coordinates$host <- pheno_mat$pheno
ggplot(pca_coordinates, aes(x = Dim.1, y = Dim.2, color = host)) +
  geom_point(size = 3) +
  labs(title = "PCA unadjusted", x = "PC1", y = "PC2") +
  theme_minimal()

# plot the adjusted dataset
bag <- get_bags(pheno_mat,
                dist_mat_or_tree = tree, 
                bag_size = c(100, 100, 100, 100, 100)
                )

boot_pheno_mat <- data.frame(ids = bag$bags[[1]],
                             pheno = pheno_mat$pheno[match(bag$bags[[1]], pheno_mat$ids)])
boot_bin_mat<- as.data.frame(bin_mat[match(bag$bags[[1]], rownames(bin_mat)), ])
colnames(boot_bin_mat) <- rep("col", ncol(boot_bin_mat))
jaccard_dist <- 1 - proxy::dist(as.matrix(boot_bin_mat), method = "Jaccard")
pca_result <- FactoMineR::PCA(as.matrix(jaccard_dist), graph = FALSE)
pca_coordinates <- as.data.frame(pca_result$ind$coord)
pca_coordinates$host <- boot_pheno_mat$pheno
ggplot(pca_coordinates, aes(x = Dim.1, y = Dim.2, color = host)) +
  geom_point(size = 3) +
  labs(title = "PCA adjusted", x = "PC1", y = "PC2") +
  theme_minimal()

