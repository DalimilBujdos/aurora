#' Function that run GWAS analysis
#'
#' @param pheno_mat Data frame that contains unique indexes in the first column
#'                  and the phenotype classes in the second column. The unique indexes
#'                  should contain only letters, numbers and special signs "_", ".".
#' @param bin_mat Binary matrix containing the genomic variants. See \code{type_bin_mat} to check what can be supplied as a binary matrix.
#' @param type_bin_mat Specifies the type of binary matrix. Options are: "panaroo"|"roary"; Exact system path to a csv file containing a pangenome matrix
#'                in Roary format (gene_presence_absence_roary.csv) or the pangenome matrix loaded as data frame. The pangenome in this format is produced by both Roary and Panaroo.
#'                "DRAM"; Exact system path to a file containing the summary of metabolism produced by DRAM (metabolism_summary.xlsx). "SCARAP"; Exact system path to a file pangenome.tsv
#'                or a dataframe containing the results of a pangenome tool SCARAP. "custom"; Data frame or matrix containing custom binary matrix. The features should be in columns and the strains in rows.
#'                "k-mers"|"unitigs"; Exact system path to a .gz file containing unitigs called by unitig-counter or k-mers by fsm-lite. "PIRATE"; Exact system path to a file containing the pangenome
#'                produced by PIRATE (PIRATE.gene_families.tsv) or the pangenome loaded as a dataframe. "SNPs"; Exact system path to a VCF file. aurora expect these columns:
#'                #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT and a column for each analysed strain.
#' @param which_snps There are two options: "biallelic" and  "all_alleles". Setting which_snps = "biallelic" will remove all SNPs that
#'                   have more than one alternative allele. Setting which_snps = "all_alleles" will create a new column for every alternative allele.
#' @param bagging Bagging algorithm applied to capture the population structure.
#'                Select from "random_walk" or "phylogenetic_walk". Default: "phylogenetic_walk"
#' @param tree  Phylogenetic tree loaded as an object of class \code{Phylo}. The tree needs to contain
#'              edge lengths and the tips should be the same as the indexes in pheno_mat.
#' @param alternative_dist_mat  Distance matrix that contains pairwise phylogenetic
#'                              distances. Row names and col names needs to be the same
#'                              as indexes in \code{pheno_mat}.
#' @param reduce_outlier  If TRUE the phylogenetic distances of a very distant strains
#'                        to the rest of the population will be reduced. Default: TRUE
#' @param cutoff_outlier  The number of standard deviations to the right of the mean
#'                        of the distance distribution. Distances with z-score higher than the cutoff
#'                        will be shrunk to the cutoff. Default: 3.
#' @param remove_strains  Character vector that specifies which strains to explicitly remove.
#' @param aurora_results  Output of a function \code{aurora_pheno}.
#' @param mode  Specifies how to handle results from \code{aurora_pheno}. Value can be either "consensus"
#'              or "strict". Consensus mode removes only strains that were found to be allochthonous
#'              by all ML tools used in \code{aurora_pheno}. Strict mode removes a strain if it was identified
#'              as allochthonous by at least one ML tool. Default: "consensus".
#' @param rm_non_typical  Specifies if strains that were labelled as non-typical should be treated as allochthonous.
#'                        Default: FALSE. WARNING: setting to TRUE may remove many strains.
#' @param use_rf  If set to TRUE, than results from Random Forest are considered. Default: TRUE.
#' @param use_ada If set to TRUE, than results from AdaBoost are considered. Default: TRUE.
#' @param use_log If set to TRUE, than results from log regression are considered. Default: TRUE.
#' @param use_CART  If set to TRUE, than results from CART are considered. Default: TRUE.
#' @param bag_size  The size of the bag for each class. Default: NA. If NA than the
#'                  bag_size is 1000 for each class of the phenotype. Provide the size
#'                  as a number for each class i.e., c(50, 50, 50) for a phenotype with three classes.
#' @param max_per_bag Maximum number of times a strain can be repeated in the bag. Default: Infinite.
#' @param get_bagging_counts  If set to TRUE then the number of times each stain appears in the
#'                            bags will be outputted.
#' @param write_data  Indicates whether aurora should write the data in a directory specified
#'                    by an argument save_dir. Default: TRUE.
#' @param save_dir An exact system path to the directory where the result will be written.
#' @param run_hogwash If set to TRUE then mGWAS tool hogwash will be run with the same input as \code{aurora_GWAS}. Default: FALSE.
#' @param run_treeWAS If set to TRUE then mGWAS tool TreeWAS will be run with the same input as \code{aurora_GWAS}. Default: FALSE.
#' @param get_scoary If set to TRUE then an input into mGWAS tool Scoary will be constructed. Default: FALSE.
#' @param get_pyseer If set to TRUE then an input into mGWAS tool Pyseer will be constructed. Default: FALSE.
#'
#' @returns The output of this function is a data frame that contains: the relative abundance of the feature in the class,
#'          standardized residuals, precision, recall and F1 values for each feature and class. Additionally the output contains a
#'          data frame that shows which strains were removed. If get_bagging_counts was set to TRUE than the output will also contain
#'          a data frame that shows how many times was each strain repeated in the bag. If the user also requested to run hogwash or
#'          TreeWAS then the result of these two tools will also be part of the output.
#' @export
#'
#' @examples
#' \dontrun{
#'   data(tree_reuteri)
#'   data(pheno_mat_reuteri)
#'   data(bin_mat_reuteri)
#'   data(aurora_pheno_results_reuteri)
#'
#'   aurora_GWAS(bin_mat = bin_mat,
#'               pheno_mat = pheno_mat,
#'               tree = tree,
#'               aurora_results = results,
#'               save_dir = "/path/to/my/dir/")
#' }
aurora_GWAS <- function(pheno_mat,
                        bin_mat = NA,
                        type_bin_mat = "roary",
                        which_snps = "all_alleles",
                        bagging = "phylogenetic_walk",
                        tree = NA,
                        alternative_dist_mat = NA,
                        reduce_outlier = TRUE,
                        cutoff_outlier = 3,
                        remove_strains = NA,
                        aurora_results = NA,
                        mode = "consensus",
                        rm_non_typical = FALSE,
                        use_rf = TRUE,
                        use_ada = TRUE,
                        use_log = TRUE,
                        use_CART = TRUE,
                        bag_size = NA,
                        max_per_bag = 100000,
                        get_bagging_counts = FALSE,
                        write_data = TRUE,
                        save_dir = NA,
                        run_hogwash = FALSE,
                        run_treeWAS = FALSE,
                        get_scoary = FALSE,
                        get_pyseer = FALSE) {

  if (is.na(save_dir) & write_data == TRUE) {
    stop("Specify the 'save_dir' or set 'write_data = FALSE'")
  }

  # if the user disables all algos then remove the AURORA list
  if (length(which(c(use_rf, use_ada, use_log, use_CART))) == 0) {
    aurora_results <- NA
  }
  # prepare pheno_mat
  pheno_mat <- as.data.frame(pheno_mat)
  names(pheno_mat) <- c("ids", "pheno")
  # remove the strains that the user requested now if the aurora results were not provided
  if (!is.list(aurora_results) & is.character(remove_strains)) {
    if (write_data == TRUE) {
      dir <- paste(save_dir, "aurora_results_removed_strains_", Sys.Date(), ".csv", sep = "")
      df_removed_strains <- pheno_mat
      colnames(df_removed_strains) <- c("strain_ids", "removed")
      df_removed_strains$removed <- "stay"
      df_removed_strains$removed[is.element(df_removed_strains$strain_ids, remove_strains)] <- "remove"
      write_csv(df_removed_strains, dir)
    }
    pheno_mat <- pheno_mat[!is.element(pheno_mat$ids, remove_strains), ]
  }
  lookup_tbl <- setNames(unique(pheno_mat$pheno),
                         1:length(unique(pheno_mat$pheno)))
  pheno_mat$pheno <- as.factor(names(lookup_tbl)[match(pheno_mat$pheno, lookup_tbl)])
  phenotypes <- as.factor(1:length(unique(pheno_mat$pheno)))
  # prepare bin_mat based on input source
  pangenome <- FALSE
  DRAM <-  FALSE
  custom <- FALSE
  SNPs <- FALSE
  scarap <- FALSE
  unitigs <- FALSE
  PIRATE <- FALSE
  if (type_bin_mat == "roary" | type_bin_mat == "panaroo") {
    bin_mat <- get_pangenome(bin_mat, pheno_mat)
    pangenome <- TRUE
    annotations <- bin_mat$annotations
    bin_mat <- bin_mat$bin_mat
  }

  if (type_bin_mat == "PIRATE") {
    bin_mat <- get_PIRATE(bin_mat, pheno_mat)
    PIRATE <- TRUE
    annotations <- bin_mat$annotations
    bin_mat <- bin_mat$bin_mat
  }
  if (type_bin_mat == "DRAM") {
    bin_mat <- get_DRAM(bin_mat, pheno_mat)
    DRAM <- TRUE
  }

  if (type_bin_mat == "SNPs") {
    bin_mat <- vcf_to_bin(bin_mat, pheno_mat, which_snps)
    SNPs <- TRUE
  }

  if (type_bin_mat == "unitigs" | type_bin_mat == "k-mers") {
    bin_mat <- unitigs_to_bin(bin_mat, pheno_mat)
    unitigs <- TRUE
  }

  if (type_bin_mat == "SCARAP") {
    bin_mat <- scarap_to_bin(bin_mat, pheno_mat)
    scarap <- TRUE
  }

  if (type_bin_mat == "custom") {
    bin_mat <- as.data.frame(bin_mat)
    custom <- TRUE
    # order it the same way as pheno_mat
    bin_mat <- bin_mat[match(pheno_mat$ids, rownames(bin_mat)), ]
  }

  if (is.list(tree)) {
    # check that phylo tree has the correct leafs. If not remove them and print warning
    drop_tips <- setdiff(tree$tip.label, pheno_mat$ids)
    if (length(drop_tips) > 0) {
      tree <- ape::drop.tip(tree, tip = drop_tips)
      warning("The tree contains strains that are not in pheno_mat. ",
              "These strains (", length(drop_tips), ") were removed. ",
              "However, for optimal performance the tree should be recalculated without these strains.\n")
    }
    # root the tree if not rooted
    if (!ape::is.rooted(tree)) {tree <- phangorn::midpoint(tree)}

    pheno_tips <- setdiff(pheno_mat$ids, tree$tip.label)
    if (length(pheno_tips) > 0) {
      messa <- paste("pheno_mat contains strains:", paste(pheno_tips, collapse = " "), "that are not in phylogenetic tree.", sep = " ")
      stop(messa, "\n")
    }
  }


  result <- distances(tree,
                      alternative_dist_mat,
                      bagging,
                      reduce_outlier,
                      cutoff_outlier)

  dist_mat_or_tree <- result$dist_mat_or_tree

  # remove strains based on AURORA results. There are two modes. Strict removes
  # ALL strains where predicted does not match observed and addtionaly all strains with
  # "atypical" label.
  # Consensus removes only strains where all provided results say that predicted does not observed
  if (is.list(aurora_results)) {
    strains_rmv_rf <- c()
    strains_rmv_ada <- c()
    strains_rmv_log <- c()
    strains_rmv_CART <- c()

    mat_to_df <- function(pval_mat) {
      # this function transformes p-val mats into a better format
      combos <- as.data.frame(t(combn(rownames(pval_mat), 2)))
      combos$dist1 <- rep(NA, nrow(combos))
      combos$dist2 <- rep(NA, nrow(combos))
      for (i in 1:nrow(combos)) {
        combos$dist1[i] <- pval_mat[rownames(pval_mat) == combos$V1[i],
                                    colnames(pval_mat) == combos$V2[i]]

        combos$dist2[i] <- pval_mat[rownames(pval_mat) == combos$V2[i],
                                    colnames(pval_mat) == combos$V1[i]]
      }
      return(combos)
    }

    if (is.list(aurora_results$results$results_random_forest) & use_rf) {
      pval_mat <- aurora_results$results$results_random_forest$p_val_mat
      df_pval_mat <- mat_to_df(pval_mat)
      decide <- any(df_pval_mat$dist1 > 0.05 & df_pval_mat$dist2 > 0.05)
      if (decide) {
        messa <- paste("Random forest cannot distinguish between some classes.",
                       "Consider removing random forest results by setting 'use_rf = FALSE'")
        warning(messa)
      }
      res_classif <- aurora_results$results$results_random_forest$auto_allo_results
      index <- c()
      for (i in 1:nrow(res_classif)) {
        if (res_classif$predicted[i] != res_classif$observed[i] & !grepl("not typical", res_classif$predicted[i])) {
          index <- c(index, i)
        }
      }
      # if rm_non_typical is TRUE remove also strains that are not typical
      if (rm_non_typical == TRUE) {
        res_classif1 <- res_classif[, c(-1, -(ncol(res_classif)-1), -ncol(res_classif))]
        vect <- apply(res_classif1, 1, paste, collapse = " ")
        index <- c(index, which(grepl("not typical", vect)))
      }

      cat("Random forest removes", length(index), "strains \n")
      strains_rmv_rf <- res_classif$strain_ids[index]
    }

    if (is.list(aurora_results$results$results_adaboost) & use_ada) {
      pval_mat <- aurora_results$results$results_adaboost$p_val_mat
      df_pval_mat <- mat_to_df(pval_mat)
      decide <- any(df_pval_mat$dist1 > 0.05 & df_pval_mat$dist2 > 0.05)
      if (decide) {
        messa <- paste("Adaboost cannot distinguish between some classes.",
                       "Consider removing adaboost results by setting 'use_ada = FALSE'")
        warning(messa)
      }
      res_classif <- aurora_results$results$results_adaboost$auto_allo_results
      index <- c()
      for (i in 1:nrow(res_classif)) {
        if (res_classif$predicted[i] != res_classif$observed[i] & !grepl("not typical", res_classif$predicted[i])) {
          index <- c(index, i)
        }
      }
      if (rm_non_typical == TRUE) {
        res_classif1 <- res_classif[, c(-1, -(ncol(res_classif)-1), -ncol(res_classif))]
        vect <- apply(res_classif1, 1, paste, collapse = " ")
        index <- c(index, which(grepl("not typical", vect)))
      }

      cat("Adaboost removes", length(index), "strains \n")
      strains_rmv_ada <- res_classif$strain_ids[index]
    }

    if (is.list(aurora_results$results$results_log_reg) & use_log) {
      pval_mat <- aurora_results$results$results_log_reg$p_val_mat
      df_pval_mat <- mat_to_df(pval_mat)
      decide <- any(df_pval_mat$dist1 > 0.05 & df_pval_mat$dist2 > 0.05)
      if (decide) {
        messa <- paste("Log regression cannot distinguish between some classes.",
                       "Consider removing log regression results by setting 'use_log = FALSE'")
        warning(messa)
      }
      res_classif <- aurora_results$results$results_log_reg$auto_allo_results
      index <- c()
      for (i in 1:nrow(res_classif)) {
        if (res_classif$predicted[i] != res_classif$observed[i] & !grepl("not typical", res_classif$predicted[i])) {
          index <- c(index, i)
        }
      }

      if (rm_non_typical == TRUE) {
        res_classif1 <- res_classif[, c(-1, -(ncol(res_classif)-1), -ncol(res_classif))]
        vect <- apply(res_classif1, 1, paste, collapse = " ")
        index <- c(index, which(grepl("not typical", vect)))
      }

      cat("Log regression removes", length(index), "strains \n")
      strains_rmv_log <- res_classif$strain_ids[index]
    }

    if (is.list(aurora_results$results$results_CART) & use_CART) {
      pval_mat <- aurora_results$results$results_CART$p_val_mat
      df_pval_mat <- mat_to_df(pval_mat)
      decide <- any(df_pval_mat$dist1 > 0.05 & df_pval_mat$dist2 > 0.05)
      if (decide) {
        messa <- paste("CART cannot distinguish between some classes.",
                       "Consider removing CART results by setting 'use_CART = FALSE'")
        warning(messa)
      }
      res_classif <- aurora_results$results$results_CART$auto_allo_results
      index <- c()
      for (i in 1:nrow(res_classif)) {
        if (res_classif$predicted[i] != res_classif$observed[i] & !grepl("not typical", res_classif$predicted[i])) {
          index <- c(index, i)
        }
      }
      if (rm_non_typical == TRUE) {
        res_classif1 <- res_classif[, c(-1, -(ncol(res_classif)-1), -ncol(res_classif))]
        vect <- apply(res_classif1, 1, paste, collapse = " ")
        index <- c(index, which(grepl("not typical", vect)))
      }

      cat("CART removes", length(index), "strains \n")
      strains_rmv_CART <- res_classif$strain_ids[index]
    }

    if (mode == "strict") {
      strains_rmv <- c()
      for (i in 1:nrow(pheno_mat)) {
        rmv <- FALSE
        if (any(strains_rmv_rf == pheno_mat$ids[i]) & use_rf) {rmv <- TRUE}
        if (any(strains_rmv_ada == pheno_mat$ids[i]) & use_ada) {rmv <- TRUE}
        if (any(strains_rmv_log == pheno_mat$ids[i]) & use_log) {rmv <- TRUE}
        if (any(strains_rmv_CART == pheno_mat$ids[i]) & use_CART) {rmv <- TRUE}
        strains_rmv <- c(strains_rmv, rmv)
      }
      # if the user want other strains to be removed remove them here
      strains_rmv2 <- is.element(pheno_mat$ids, remove_strains)
      strains_rmv <- strains_rmv | strains_rmv2

      x <- rep("stay", nrow(pheno_mat))
      x[strains_rmv] <- "remove"
      df_removed_strains <- data.frame(strain_ids = pheno_mat$ids, removed = x)
      if (write_data == TRUE) {
        dir <- paste(save_dir, "aurora_results_removed_strains_", Sys.Date(), ".csv", sep = "")
        write_csv(df_removed_strains, dir)
      }

      cat("Strict mode removed in total", length(which(strains_rmv)),
          "strains. There are", nrow(pheno_mat) - length(which(strains_rmv)), "strains left. \n")
      if ((nrow(pheno_mat) - length(which(strains_rmv))) == 0) {
        stop("Strict mode removed all strains. Try to set 'mode = consensus'")
      }

      if (bagging == "phylogenetic_walk") {
        dist_mat_or_tree <- drop.tip(dist_mat_or_tree, tip = pheno_mat$ids[strains_rmv])
      }

      pheno_mat <- pheno_mat[!strains_rmv, ]
      bin_mat <- bin_mat[!strains_rmv, ]
      if (bagging == "random_walk") {
        dist_mat <- dist_mat_or_tree[is.element(pheno_mat$ids, colnames(dist_mat_or_tree)),
                                     is.element(pheno_mat$ids, rownames(dist_mat_or_tree))]
      }
    }

    if (mode == "consensus") {
      strains_rmv <- c()
      no_classif <- length(which(c(is.list(aurora_results$results$results_random_forest) & use_rf,
                                   is.list(aurora_results$results$results_adaboost) & use_ada,
                                   is.list(aurora_results$results$results_log_reg) & use_log,
                                   is.list(aurora_results$results$results_CART) & use_CART)))
      for (i in 1:nrow(pheno_mat)) {
        count <- 0
        if (any(strains_rmv_rf == pheno_mat$ids[i]) & use_rf) {count <- count+1}
        if (any(strains_rmv_ada == pheno_mat$ids[i]) & use_ada) {count <- count+1}
        if (any(strains_rmv_log == pheno_mat$ids[i]) & use_log) {count <- count+1}
        if (any(strains_rmv_CART == pheno_mat$ids[i]) & use_CART) {count <- count+1}

        if (count == no_classif) {
          strains_rmv <- c(strains_rmv, TRUE)
        } else {
          strains_rmv <- c(strains_rmv, FALSE)
        }
      }
      # if the user want other strains to be removed remove them here
      strains_rmv2 <- is.element(pheno_mat$ids, remove_strains)
      strains_rmv <- strains_rmv | strains_rmv2

      x <- rep("stay", nrow(pheno_mat))
      x[strains_rmv] <- "remove"
      df_removed_strains <- data.frame(strain_ids = pheno_mat$ids, removed = x)
      if (write_data == TRUE) {
        dir <- paste(save_dir, "aurora_results_removed_strains_", Sys.Date(), ".csv", sep = "")
        write_csv(df_removed_strains, dir)
      }

      cat("Consensus mode removed in total", length(which(strains_rmv)),
          "strains. There are", nrow(pheno_mat) - length(which(strains_rmv)),"strains left. \n")
      if ((nrow(pheno_mat) - length(which(strains_rmv))) == 0) {
        stop("Consensus mode removed all strains. It looks like there is no adaptation")
      }

      if (bagging == "phylogenetic_walk") {
        dist_mat_or_tree <- drop.tip(dist_mat_or_tree, tip = pheno_mat$ids[strains_rmv])
      }

      pheno_mat <- pheno_mat[!strains_rmv, ]
      bin_mat <- bin_mat[!strains_rmv, ]
      if (bagging == "random_walk") {
        dist_mat_or_tree <- dist_mat_or_tree[is.element(pheno_mat$ids, colnames(dist_mat_or_tree)),
                                             is.element(pheno_mat$ids, rownames(dist_mat_or_tree))]
      }
    }
  }

  # first, calculate the frequency of each gf in each phenotype
  calcul_ratio <- function(gf, pheno_mat, phenotypes) {
    calcul_ratio1 <- function(phenotype, gf, pheno_mat) {
      jmenovatel <- length(which(pheno_mat$pheno == phenotype))
      citatel <- pheno_mat$ids[pheno_mat$pheno == phenotype]
      citatel <- sum(gf[is.element(names(gf), citatel)])
      return(citatel/jmenovatel)
    }
    result <- sapply(phenotypes, calcul_ratio1, gf = gf, pheno_mat = pheno_mat)
    return(result)
  }

  df_ratios <- apply(as.matrix(bin_mat), 2, calcul_ratio,
                      pheno_mat = pheno_mat,
                      phenotypes = phenotypes)

  df_ratios <- as.data.frame(t(df_ratios))
  colnames(df_ratios) <- paste(lookup_tbl, "ratio")

  # the maximum bag size is 1000
  if (!is.numeric(bag_size)) {
    bag_size <- rep(1000, length(phenotypes))
  }

  if (bagging == "random_walk") {
    result_get_bags <- get_bags(pheno_mat = pheno_mat,
                                dist_mat_or_tree = dist_mat_or_tree,
                                bagging = bagging,
                                bag_size = bag_size,
                                no_rf = 1,
                                max_per_bag = max_per_bag)
  }

  if (bagging == "phylogenetic_walk") {
    result_get_bags <- get_bags(pheno_mat = pheno_mat,
                                dist_mat_or_tree = dist_mat_or_tree,
                                bagging = bagging,
                                bag_size = bag_size,
                                no_rf = 1,
                                max_per_bag = max_per_bag)
  }

  bags <- result_get_bags$bags[[1]]
  bags_counts <- result_get_bags$count
  names(bags_counts) <- lookup_tbl

  if (get_bagging_counts == TRUE) {
    # calculate a df that will show how many times each strain was used to construct the bags
    # this can then be examines in iTOL
    result_iTOL <- data.frame()
    for (i in 1:length(bags_counts)) {
      no_strains <- length(bags_counts[[i]])
      df_tmp <- data.frame(species = rep(lookup_tbl[i], no_strains),
                           strain_ids = names(bags_counts[[i]]),
                           counts = bags_counts[[i]])
      result_iTOL <- rbind(result_iTOL, df_tmp)
    }

    if (write_data == TRUE) {
      dir <- paste(save_dir, "bagging_results_", Sys.Date(), ".csv", sep = "")
      write_csv(result_iTOL, dir)
    }
  }

  # this is needed to calculate the chisq accueately
  host_counts_original <- as.data.frame(table(pheno_mat$pheno))
  # calculate sensitivity specificity, F1, and chisq test for each gf in each pheno
  calcul_stats <- function(gf, pheno_mat, phenotypes, host_counts_original) {
    # calculate the chisq and the residuals
    contig_mat <- matrix(NA, nrow = 2,
                         ncol = length(phenotypes),
                         dimnames = list(c("present", "absent"), as.character(phenotypes)))
    for (i in 1:length(phenotypes)) {
      gf_pheno <- gf[which(pheno_mat$pheno == phenotypes[i])]
      present_ratio <- length(which(gf_pheno == 1))/length(gf_pheno)
      absent_ratio <- length(which(gf_pheno == 0))/length(gf_pheno)
      host_no <- host_counts_original$Freq[host_counts_original$Var1 == phenotypes[i]]
      contig_mat[1,i] <- round(host_no*present_ratio)
      contig_mat[2,i] <- round(host_no*absent_ratio)
    }
    test_res <- suppressWarnings(chisq.test(contig_mat, simulate.p.value = TRUE))
    p_val <- test_res$p.value
    stand_res <- test_res$stdres # maybe just non-std residuals

    calcul_stats1 <- function(phenotype, gf, pheno_mat) {
      index <- which(pheno_mat$pheno == phenotype)
      TP <- sum(gf[index])
      FP <- sum(gf[-index])
      precision <- TP/(TP + FP)
      FN <- length(which(gf[index] == 0))
      recall <- TP/(TP + FN)
      F1 <- 2*((precision*recall)/(precision+recall))
      return(c(precision, recall, F1))
    }
    result <- sapply(phenotypes, calcul_stats1, gf = gf, pheno_mat = pheno_mat)
    result <- t(result)
    result <- cbind(stand_res[1,], result)

    colnames(result) <- c("std_residual", "precision", "recall", "F1")
    rownames(result) <- as.character(phenotypes)
    result_one <- c()
    for (i in 1:nrow(result)) {
      line <- result[i,]
      names(line) <- c(paste0(rownames(result)[i],"_",names(line)[1]),
                       paste0(rownames(result)[i],"_",names(line)[2]),
                       paste0(rownames(result)[i],"_",names(line)[3]),
                       paste0(rownames(result)[i],"_",names(line)[4]))
      result_one <- c(result_one, line)
    }
    result_one <- c(p_val, result_one)
    names(result_one)[1] <- "chisq_pval"
    return(result_one)
  }

  # get the right bin_mat and pheno_mat
  bin_mat_bags <- bin_mat[match(bags, rownames(bin_mat)), ]
  pheno_mat_bags <-  pheno_mat[match(bags, pheno_mat$ids), ]

  mat_results <- apply(as.matrix(bin_mat_bags),2, calcul_stats,
                       pheno_mat = pheno_mat_bags,
                       phenotypes = phenotypes,
                       host_counts_original = host_counts_original)

  df_vals <- as.data.frame(t(mat_results))

  rep4 <- function(x) {
    rep(x, 4)
  }
  cols1 <- as.character(sapply(lookup_tbl, rep4))
  colnames(df_vals) <- c("chisq p.val", paste(cols1, c("std_residual", "precision", "recall", "F1")))
  df_final <- cbind(data.frame(Variant = rownames(df_ratios)), df_ratios, df_vals)

  # if pangenome was used then add the annotation as a second column
  if (pangenome == TRUE) {
    df_final <- cbind(df_final[,1], annotations, df_final[,2:ncol(df_final)])
    colnames(df_final)[1] <- "Variant"
  }

  # if PIRATE was used then add the annotation as a second column
  if (PIRATE == TRUE) {
    df_final <- cbind(df_final[,1], annotations, df_final[,2:ncol(df_final)])
    colnames(df_final)[1] <- "Variant"
  }

  if (write_data == TRUE) {
    dir <- paste(save_dir, "aurora_results_GWAS_", Sys.Date(), ".csv", sep = "")
    write_csv(df_final, dir)
  }

  results <- list(GWAS_results = df_final)
  if (is.list(aurora_results)) {
    results <- rlist::list.append(results, removed_strains = df_removed_strains)
  }
  if (get_bagging_counts == TRUE) {
    results <- rlist::list.append(results, bagging_results = result_iTOL)
  }
  # run Hogwash if requested
  if (run_hogwash == TRUE) {
    if (!require("hogwash", character.only = TRUE)) {
      warning("Hogwash was not run because the library is not installed or attached. Please see https://github.com/katiesaund/hogwash.\n")
    } else if (!is.na(alternative_dist_mat)) {
      warning("Hogwash was not run because a phylogenetic tree needs to be supplied.\n")
    } else if (nlevels(phenotypes) != 2){
      warning("Hogwash was not run because the number of phenotype classes is not two.\n")
    } else {
      # Matrix must be formatted with samples in matrix in the same order as dist_mat_or_tree$tip.label
      select_bin_mat <- c()
      select_pheno_mat <- c()

      for (i in 1:length(dist_mat_or_tree$tip.label)) {
        select_bin_mat <- c(select_bin_mat, which(rownames(bin_mat) == dist_mat_or_tree$tip.label[i]))
        select_pheno_mat <- c(select_pheno_mat, which(pheno_mat$ids == dist_mat_or_tree$tip.label[i]))
      }

      bin_mat_hog <- bin_mat[select_bin_mat,]
      pheno_mat_hog <- pheno_mat[select_pheno_mat, ]
      # convert pheno_mat to matrix
      pheno_mat_mat <- matrix(data = as.numeric(pheno_mat_hog$pheno), nrow = nrow(pheno_mat_hog))
      rownames(pheno_mat_mat) <- pheno_mat_hog$ids
      colnames(pheno_mat_mat) <- "pheno"

      # tree edges cannot have value == 0 and NA
      tree_hog <- dist_mat_or_tree
      tree_hog$edge.length[is.na(tree_hog$edge.length)] <- min(tree_hog$edge.length[tree_hog$edge.length > 0 &
                                                                                      !is.na(tree_hog$edge.length)])
      tree_hog$edge.length[tree_hog$edge.length == 0] <- min(tree_hog$edge.length[tree_hog$edge.length > 0])

      # not all trees have bootstrap values. Set all bootstrap vals to 100 to avoid errors
      tree_hog$node.label <- rep("100",  length(tree_hog$node.label))

      cat("Running Hogwash\n")
      res_hogwash <- hogwash::hogwash(pheno =  pheno_mat_mat,
                                      geno = bin_mat_hog,
                                      tree = tree_hog,
                                      test = "both",
                                      dir = save_dir,
                                      fdr = 0.99,
                                      perm = 5000)

      results <- rlist::list.append(results, results_hogwash = res_hogwash)
    }
  }

  # run TreeWAS if requested
  if (run_treeWAS == TRUE) {
    if (!is.na(alternative_dist_mat)) {
      warning("TreeWAS was not run because a phylogenetic tree needs to be supplied.\n")
    } else {
      # prepare the pheno_mat to a format for TreeWAS
      pheno_mat_WAS <- pheno_mat$pheno
      pheno_mat_WAS <- unname(lookup_tbl[match(pheno_mat_WAS,names(lookup_tbl))])
      names(pheno_mat_WAS) <- pheno_mat$ids

      # just a check that names in bin_mat equals names in pheno_mat and also to remove no observation variants
      bin_mat_WAS <- bin_mat[is.element(rownames(bin_mat), names(pheno_mat_WAS)), ]
      bin_mat_WAS <- bin_mat_WAS[, which(colSums(bin_mat_WAS) != 0 & colSums(bin_mat_WAS) != nrow(bin_mat_WAS))]
      # NA means that treeWAS will treat it as binary
      phen_type <- ifelse(length(phenotypes) == 2, NA,"categorical")
      cat("Running TreeWAS\n")
      result_treeWAS <- treeWAS::treeWAS(bin_mat_WAS,
                                         pheno_mat_WAS,
                                         phen.type = phen_type,
                                         tree = dist_mat_or_tree,
                                         n.snps.sim = 100*ncol(bin_mat),
                                         p.value = 0.99,
                                         p.value.correct = "fdr",
                                         plot.tree = FALSE, # makes the output too large
                                         plot.manhattan = TRUE,
                                         plot.null.dist = TRUE,
                                         plot.dist = TRUE,
                                         chunk.size = 1000,
                                         filename.plot = paste0(save_dir, "TreeWAS_result_", Sys.Date(),".pdf"))

      results <- rlist::list.append(results, result_treeWAS = result_treeWAS)

      terminal <- data.frame(ID_terminal = names(result_treeWAS$terminal$corr.dat),
                             score_terminal = result_treeWAS$terminal$corr.dat,
                             p_val_terminal = result_treeWAS$terminal$p.vals)
      terminal <- terminal[order(terminal$score_terminal), ]

      simul <- data.frame(ID_simultaneous = names(result_treeWAS$simultaneous$corr.dat),
                          score_simultaneous = result_treeWAS$simultaneous$corr.dat,
                          p_val_simultaneous = result_treeWAS$simultaneous$p.vals)
      simul <- simul[order(simul$score_simultaneous), ]

      subseq <- data.frame(ID_subsequent = names(result_treeWAS$subsequent$corr.dat),
                           score_subsequent = result_treeWAS$subsequent$corr.dat,
                           p_val_subsequent = result_treeWAS$subsequent$p.vals)
      subseq <- subseq[order(subseq$score_subsequent), ]

      df_save <- cbind(terminal, simul, subseq)

      results <- rlist::list.append(results, result_treeWAS_ordered = df_save)
      if (is.character(save_dir)) {
        readr::write_csv(df_save, paste0(save_dir, "/treeWAS_results_", Sys.Date(), ".csv"))
      }
    }
  }

  # make scoary data if requested
  if (get_scoary == TRUE & is.character(save_dir)) {
    if (!is.na(alternative_dist_mat)) {
      warning("Scoary input data could not be constructed because a phylogenetic tree needs to be supplied.\n")
    } else {
      pheno_mat_sc <- pheno_mat
      pheno_mat_sc$pheno <- unname(lookup_tbl[match(pheno_mat_sc$pheno,names(lookup_tbl))])
      colnames(pheno_mat_sc) <- c("", "ids")
      readr::write_csv(pheno_mat_sc, paste0(save_dir, "/pheno_mat_scoary.csv"))

      ape::write.tree(dist_mat_or_tree, paste0(save_dir, "/tree_scoary.treefile"))

      bin_mat_scoary <- bin_mat[, which(colSums(bin_mat) != 0)]
      bin_mat_scoary <- t(bin_mat_scoary)
      buffer_matrix <- as.data.frame(matrix(NA, nrow = nrow(bin_mat_scoary), ncol = 13))
      colnames(buffer_matrix) <- paste0(rep("buffer", 13), 1:13)
      bin_mat_scoary <- cbind(data.frame(Gene = rownames(bin_mat_scoary)), buffer_matrix, bin_mat_scoary)
      readr::write_csv(bin_mat_scoary, paste0(save_dir, "/bin_mat_scoary.csv"))
    }
  }

  # make Pyseer data if requested
  if (get_pyseer == TRUE & is.character(save_dir)) {
    if (!is.na(alternative_dist_mat)) {
      warning("Pyseer input data could not be constructed because a phylogenetic tree needs to be supplied.\n")
    } else {
      pheno_mat_py <- pheno_mat
      pheno_mat_py$pheno <- unname(lookup_tbl[match(pheno_mat_py$pheno,names(lookup_tbl))])
      colnames(pheno_mat_py) <- c("ids", "ids")
      readr::write_delim(pheno_mat_py, paste0(save_dir, "/pheno_mat_pyseer.tsv"), delim = "\t")

      ape::write.tree(dist_mat_or_tree, paste0(save_dir, "/tree_scoary.treefile"))

      bin_mat_pyseer <- bin_mat[, which(colSums(bin_mat) != 0)]
      bin_mat_pyseer <- t(bin_mat_pyseer)
      bin_mat_pyseer <- cbind(data.frame(Gene = rownames(bin_mat_pyseer)), bin_mat_pyseer)

      readr::write_delim(bin_mat_pyseer, paste0(save_dir, "/bin_mat_pyseer.tsv"))
    }
  }
  cat("Finished \n")
  return(results)
}

