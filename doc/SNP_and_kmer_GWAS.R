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
#                          bin_mat = bin_mat_snps,
#                          type_bin_mat = "SNPs",
#                          which_snps = "all_alleles",
#                          tree = tree,
#                          low_perc_cutoff = 5,
#                          upp_perc_cutoff = 95,
#                          jaccard_filter = TRUE,
#                          run_chisq = TRUE,
#                          cutoff_chisq = 0.1,
#                          write_data = FALSE)
#  

## -----------------------------------------------------------------------------
data("aurora_pheno_results_neisseria") # load the precalculated results
res <- aurora_GWAS(bin_mat = bin_mat_snps,
                   type_bin_mat = "SNPs",
                   pheno_mat = pheno_mat,
                   tree = tree,
                   aurora_results = results,
                   mode = "consensus",
                   write_data = FALSE)

## -----------------------------------------------------------------------------
res$GWAS_results <- res$GWAS_results[order(res$GWAS_results$`resistant std_residual`, decreasing = TRUE),]
knitr::kable(res$GWAS_results[1:5,c(1:4,9)], format = "html", table.attr = "class='table table-striped'")

