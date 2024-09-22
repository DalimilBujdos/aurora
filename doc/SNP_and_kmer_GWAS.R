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

## -----------------------------------------------------------------------------
data("tree_neisseria")
data("pheno_mat_neisseria")
bin_mat_snps <- paste(system.file(package="aurora"), ".inst/extdata/neisseria_snps.vcf", sep="/")

## ----eval = FALSE-------------------------------------------------------------
#  results <- aurora_pheno(pheno_mat = pheno_mat,
#                          bin_mat_snps = bin_mat_snps,
#                          which_snps = "all_alleles",
#                          tree = tree,
#                          low_perc_cutoff = 5,
#                          upp_perc_cutoff = 95,
#                          jaccard_filter = TRUE,
#                          run_chisq = TRUE,
#                          cutoff_chisq = 0.1,
#                          condaenv_path = "/path/to/conda/env",
#                          write_data = FALSE)
#  

## ----eval = FALSE-------------------------------------------------------------
#  data("aurora_pheno_results_neisseria")
#  res <- aurora_GWAS(bin_mat_snps = bin_mat_snps,
#                     pheno_mat = pheno_mat,
#                     tree = tree,
#                     aurora_results = results,
#                     mode = "consensus",
#                     write_data = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  res$GWAS_results <- res$GWAS_results[order(res$GWAS_results$`resistant std_residual`, decreasing = TRUE),]
#  
#  head(res$GWAS_results[,c(9:12)], 5)

