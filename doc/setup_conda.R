## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = FALSE,
  warning = FALSE,
  message = FALSE,
  cache.lazy = FALSE
)

## ----eval = FALSE-------------------------------------------------------------
#  library(aurora)
#  # set up conda environment using packaged .yaml file
#  conda_path <- setup_conda() # make sure that the path is in this format: /path/to/dir
#  
#  # now you can use the path in aurora_pheno() and log regression and AdaBoost will be run
#  results <- aurora_pheno(pheno_mat = pheno_mat,
#                          bin_mat = bin_mat,
#                          tree = tree,
#                          condaenv_path = conda_path,
#                          save_dir = "/path/to/save/dir")
#  

