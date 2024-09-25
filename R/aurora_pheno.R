#' Main function of the aurora package. This function determines if the analysed species is adapted to the investigated phenotype and it also identifies
#' mislabeled strains. The output of this function can be used as an input into \code{aurora_GWAS} function.
#'
#' @param pheno_mat Data frame that contains unique indexes in the first column
#'                  and the phenotype classes in the second column. The unique indexes
#'                  should contain only letters, numbers and special signs "_", ".". The maximum number of unique classes is 9.
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
#'                Select from "random_walk" or "phylogenetic_walk". Default: "phylogenetic_walk".
#' @param tree  Phylogenetic tree loaded as an object of class \code{Phylo}. The tree needs to contain
#'              edge lengths and the tips should be the same as the indexes in pheno_mat.
#' @param alternative_dist_mat  Distance matrix that contains pairwise phylogenetic
#'                              distances. Row names and column names needs to be the same
#'                              as indexes in \code{pheno_mat}.
#' @param reduce_outlier  If TRUE the phylogenetic distances of very distant strains
#'                        to the rest of the population will be reduced. Default: TRUE.
#' @param cutoff_outlier  The number of standard deviations to the right of the mean
#'                        of the distance distribution. Distances with z-score higher than the cutoff
#'                        will be shrunk to the cutoff. Default: 3.
#' @param fit_parameters  Indicates if parameters to AdaBoost and Random Forest should be fitted.
#'                        Default: TRUE.
#' @param repeats Number of times the parameters are fitted. Default: 10.
#' @param random_forest Indicates whether to use Random Forest. Default: TRUE.
#' @param plot_random_forest  Indicates if the distance matrix from Random Forest
#'                            should be plotted. Default: TRUE.
#' @param sampsize  Determines the size of the random sample used to build each tree in Random Forest.
#'                  Default: 50. Ignored if \code{fit_parameters} = TRUE.
#' @param mtry  Controls the number of features considered at each split during
#'              tree building in Random Forest. Default: 200. Ignored if \code{fit_parameters} = TRUE.
#' @param ntree Determines the number of trees in Random Forest. Default: 100. Ignored if \code{fit_parameters} = TRUE.
#' @param maxnodes  Limits the maximum number of terminal nodes in each Random Forest tree. Default: 12. Ignored if \code{fit_parameters} = TRUE.
#' @param ovr_log_reg Indicates whether to use log regression with one vs rest strategy.
#'                    Default: TRUE.
#' @param C_val Regularization term for l1 penalty in log regression. Default: 0.5.
#' @param adaboost Indicates whether to use AdaBoost. Default: TRUE.
#' @param max_depth Indicates the maximum depth of the decision tree in AdaBoost. Default: 1.
#'                  Changing the value will not effect anything. It is only here because
#'                  in future version the user will be able to experiment with it.
#' @param n_estimators Determines the number of weak learners to combine in AdaBoost. Default: 500. Ignored if \code{fit_parameters} = TRUE.
#' @param learning_rate Controls the contribution of each weak learner in the final AdaBoost model. Default: 0.3. Ignored if \code{fit_parameters} = TRUE.
#' @param CART  Indicates whether to use classification and regression tree (CART) model. Default: TRUE.
#' @param CART_plot Indicates whether the proximities from CART models
#'                  should be plotted. Default: TRUE.
#' @param condaenv_path An exact system path to a conda environment that will be used for log regression and
#'                      AdaBoost. If no path is provided or if the path is set to NA then AdaBoost and log
#'                      regression will not be run.
#' @param low_perc_cutoff The lower cutoff for filtering the supplied binary matrix. Default = 3.
#'                        This means that all features present in less than 3% of all strains will be removed.
#' @param upp_perc_cutoff The upper cutoff for filtering the supplied binary matrix. Default = 99.
#'                        This means that all features present in more than 99% of all strains will be removed.
#' @param run_chisq Indicates whether chi-square filter should be run. Use only if the number of features after
#'                  the initial filtering is still large (> 10,000). Default: FALSE.
#' @param cutoff_chisq  A cutoff p-value for the chi-square test. Features with p-value
#'                      higher that the cutoff will be removed. Default: 0.1.
#' @param jaccard_filter  Grouping of correlated features based on Jaccard distance matrix and DBSCAN.
#'                        Default: FALSE. Only one grouping method (\code{jaccard_filter} or \code{hamming_filter})
#'                        can be used.
#' @param minPts_val  \code{minPts} parameter in DBSCAN. Specifies the minimum number of neighboring data
#'                    points for a point to be considered a core point. Default: 3.
#' @param eps_val \code{eps} parameter in DBSCAN. Determines the maximum distance between two points
#'                for them to be considered neighbors. Default: 0.01.
#' @param hamming_filter  Grouping of correlated features based on Hamming distance.
#'                        Default: TRUE. Only one grouping method (\code{jaccard_filter} or \code{hamming_filter})
#'                        can be used.
#' @param hamming_cutoff  Maximum Hamming distance. If features have intrarcluster
#'                        distances lower or equal to hamming_cutoff then they are grouped
#'                        into one feature. Default: 3.
#' @param ancest_rec_filter Indicates if ancestral reconstruction filter should be used.
#'                          The filter removes features broadly distributed along the
#'                          phylogenetic tree. These features are often common plasmid genes
#'                          and IS elements. Default: TRUE.
#' @param cutoff_asr  The number of standard deviations to the right of the mean
#'                    of the distribution produced by ancest_rec_filter. Features
#'                    with z-score higher than the cutoff will be removed. Default: 2.
#' @param bag_size  The size of the bag for each class. Default: NA. If NA than the
#'                  bag_size is calculated as 5* the number of strains in the class
#'                  with the fewest strains. Provide the size as a number for each class
#'                  i.e., c(50, 50, 50) for a phenotype with three classes.
#' @param no_rounds Number of times the aurora algorithm iterates. Default: 100. You can decrease
#'                  this number to lower the computational time
#' @param max_per_bag Maximum number of times a strain can be repeated in the bag. Default: NA
#'                    If NA than the max_per_bag is calculated so that none of the strains
#'                    exceeds 20% of the bag of each class.
#' @param misslabel_no  The number of mislabeled strains in the threshold calculation
#'                      phase. Default: 1. Do not modify this argument. In the future
#'                      version the user will be able to experiment with this.
#' @param write_data  Indicates whether aurora should write the data in a directory specified
#'                    by an argument \code{save_dir}. Default: TRUE.
#' @param save_dir An exact system path to the directory where the result will be written.
#'
#' @return  The output is a nested list that mainly contains a table for each ML tool.
#'          The tables show which strains were identified as autochthonous and allochthonous.
#'          The list also contains p-value matrices that show if the entire species is adapted
#'          to the investigated phenotype. If \code{write_csv} was set to \code{TRUE} and \code{save_dir}
#'          was provided then distance matrices are written to the folder. For more information about the
#'          output check https://dalimilbujdos.github.io/aurora/articles/outputs.html.
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   data(pheno_mat_reuteri)
#'   data(bin_mat_reuteri)
#'   data(tree_reuteri)
#'
#'   aurora_pheno(bin_mat = bin_mat,
#'                type_bin_mat = "panaroo",
#'                pheno_mat = pheno_mat,
#'                tree = tree,
#'                save_dir = "/path/to/my/dir/")
#' }

aurora_pheno <- function(
    pheno_mat,
    bin_mat = NA,
    type_bin_mat = "roary",
    which_snps = "all_alleles",
    bagging = "phylogenetic_walk",
    tree = NA,
    alternative_dist_mat = NA,
    reduce_outlier = TRUE,
    cutoff_outlier = 3,
    fit_parameters = TRUE,
    repeats = 10,
    random_forest = TRUE,
    plot_random_forest = TRUE,
    sampsize = 50,
    mtry = 200,
    ntree = 100,
    maxnodes = 12,
    ovr_log_reg = TRUE,
    C_val = 0.5,
    adaboost = TRUE,
    max_depth = 1,
    n_estimators = 500,
    learning_rate = 0.3,
    CART = TRUE,
    CART_plot = TRUE,
    condaenv_path = NA,
    low_perc_cutoff = 3,
    upp_perc_cutoff = 99,
    run_chisq = FALSE,
    cutoff_chisq = 0.1,
    jaccard_filter = FALSE,
    minPts_val = 3,
    eps_val = 0.01,
    hamming_filter = TRUE,
    hamming_cutoff = 3,
    ancest_rec_filter = TRUE,
    cutoff_asr = 2,
    bag_size = NA,
    no_rounds = 100,
    max_per_bag = NA,
    misslabel_no = 1,
    write_data = TRUE,
    save_dir = NA) {

  if (is.na(condaenv_path)) {
    # if the user did not provide path to the conda env then disable ada and log reg
    adaboost <- FALSE
    ovr_log_reg <- FALSE
  }

  if (adaboost == TRUE) {
    # source all functions for ada_boost
    reticulate::use_condaenv(condaenv_path)
    path <- paste(system.file(package="aurora"), ".inst/run_ada_miss.py", sep="/")
    reticulate::source_python(path)
    path <- paste(system.file(package="aurora"), ".inst/run_ada.py", sep="/")
    reticulate::source_python(path)
    if (fit_parameters == TRUE) {
      path <- paste(system.file(package="aurora"), ".inst/fit_parameters_ada.py", sep="/")
      reticulate::source_python(path)
    }
  }

  if (ovr_log_reg == TRUE) {
    # source all functions for ovr log regression
    reticulate::use_condaenv(condaenv_path)
    path <- paste(system.file(package="aurora"), ".inst/run_log_reg_miss.py", sep="/")
    reticulate::source_python(path)
    path <- paste(system.file(package="aurora"), ".inst/run_log_reg.py", sep="/")
    reticulate::source_python(path)
  }


  if (is.na(save_dir) & write_data == TRUE) {
    stop("Provide absolute path to a directory where the data will be saved or set write_data = FALSE")
  }

  # prepare pheno_mat
  pheno_mat <- as.data.frame(pheno_mat)
  names(pheno_mat) <- c("ids", "pheno")
  lookup_tbl <- setNames(unique(pheno_mat$pheno),
                         1:length(unique(pheno_mat$pheno)))
  pheno_mat$pheno <- as.factor(names(lookup_tbl)[match(pheno_mat$pheno, lookup_tbl)])
  phenotypes <- as.factor(levels(pheno_mat$pheno))
  # check that the classes are roughly balanced and that we have enough strains to run the analysis
  pheno_tbl <- as.data.frame(table(pheno_mat$pheno))
  if (any(pheno_tbl$Freq < 20)) {
    warning("Some phenotype categories have less than 20 strains. In such case results from aurora_pheno() may not be reliable.")
  }

  maxi <- pheno_tbl$Freq[which.max(pheno_tbl$Freq)]
  pheno_tbl$Freq <- maxi/pheno_tbl$Freq

  if (any(pheno_tbl$Freq < 20)) {
    warning("The number of strains in phenotype categories is not balanced. In such case results from aurora_pheno() may not be reliable.")
  }

  # prepare bin_mat based on input source (panaroo/roary, DRAM, SNPs, unitigs, kmers, scarap, custom)
  pangenome <- FALSE
  DRAM <-  FALSE
  custom <- FALSE
  SNPs <- FALSE
  unitigs <- FALSE
  scarap <- FALSE
  PIRATE <- FALSE
  if (type_bin_mat == "roary" | type_bin_mat == "panaroo") {
    bin_mat <- get_pangenome(bin_mat, pheno_mat)
    pangenome <- TRUE
    bin_mat <- bin_mat$bin_mat
  }

  if (type_bin_mat == "PIRATE") {
    bin_mat <- get_PIRATE(bin_mat, pheno_mat)
    PIRATE <- TRUE
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

  if (is.character(type_bin_mat == "scarap")) {
    bin_mat <- scarap_to_bin(bin_mat, pheno_mat)
    scarap <- TRUE
  }

  if (type_bin_mat == "custom") {
    bin_mat <- as.data.frame(bin_mat)
    custom <- TRUE
    # order it the same way as pheno_mat
    bin_mat <- bin_mat[match(pheno_mat$ids, rownames(bin_mat)), ]
    colnames(bin_mat) <- gsub("___", "_", colnames(bin_mat))
  }

  if (!is.list(tree)) {
    # remove acr filter if tree was not provided
    ancest_rec_filter <- FALSE
  } else {
    # check that phylo tree has the correct leafs. If not remove them and print warning
    drop_tips <- setdiff(tree$tip.label, pheno_mat$ids)
    if (length(drop_tips) > 0) {
      tree <- ape::drop.tip(tree, tip = drop_tips)
      warning("The tree contains strains that are not in the pheno_mat. ",
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

  # copy for the final csv calculation
  bin_mat_copy <- bin_mat

  to_bin_result <- to_bin(bin_mat,
                          pheno_mat,
                          low_perc_cutoff,
                          upp_perc_cutoff,
                          run_chisq,
                          cutoff_chisq)

  sample_names <- rownames(to_bin_result$bin_mat)

  post_proc_data <- preprocess_data(to_bin_result$bin_mat,
                                    pheno_mat,
                                    tree,
                                    ancest_rec_filter,
                                    cutoff_asr,
                                    jaccard_filter,
                                    eps_val,
                                    minPts_val,
                                    hamming_filter,
                                    hamming_cutoff)

  bin_mat <- post_proc_data$bin_mat

  # convert bin_mat to factor variable

  for (i in 1:ncol(bin_mat)) {
    index_hit <- which(bin_mat[,i] == 1)
    bin_mat[index_hit,i] <- "hit"
    index_miss <- which(bin_mat[,i] == "0")
    bin_mat[index_miss,i] <- "miss"
  }
  # will convert all strings into factors
  bin_mat <- as.data.frame(unclass(bin_mat), stringsAsFactors =TRUE)
  row.names(bin_mat) <- sample_names
  bonferroni_threshold <- 0.05/ncol(bin_mat)


  result <- distances(tree,
                      alternative_dist_mat,
                      bagging,
                      reduce_outlier,
                      cutoff_outlier)

  ################################################################################
  ################################################################################
  ################################################################################
  # calculate the allochthony distributions

  # prepare the features input
  predictors <- bin_mat
  predictors$ids <- rownames(predictors)
  # the problem with randomForest is that the colnames has to be changed to generic names
  # otherwise it run into an error
  key_gfsnames <- colnames(predictors)[colnames(predictors) != "ids"] # the original gfs name
  names(predictors)[1:(ncol(predictors) - 1)] <- paste0("V", 1:(ncol(predictors) - 1))
  # prepare the test set
  predictors_test <- predictors
  pheno_mat_test <- pheno_mat$pheno[match(predictors_test$ids, pheno_mat$ids)]
  predictors_test <- predictors[, colnames(predictors_test) != "ids"]

  # initialize a list that will hold the pheno_mats for simulations
  pheno_mats_list <- vector("list", no_rounds)
  # initialize a list that will hold the bin_mats for simulations
  bin_mats_list <- vector("list", no_rounds)
  # initialize a list that will hold the misslabel_mats for simulations
  misslabel_mats_list <- vector("list", no_rounds)

  if (CART == TRUE) {
    # initilize the objects that will hold results for CART classification
    grid_inp <- expand.grid(phenotype1 = phenotypes,
                            phenotype2 = phenotypes)
    phenotypes1 <- c()
    phenotypes2 <- c()
    existing_pairs <- c()
    for (i in 1:nrow(grid_inp)) {
      if (grid_inp[i,1] == grid_inp[i,2]) {next}
      pair1 <- paste(grid_inp[i,1],grid_inp[i,2], sep = "")
      pair2 <- paste(grid_inp[i,2],grid_inp[i,1], sep = "")
      if (!any(c(existing_pairs == pair1, existing_pairs == pair2))) {
        phenotypes1 <- c(phenotypes1, grid_inp[i,2])
        phenotypes2 <- c(phenotypes2, grid_inp[i,1])
        existing_pairs <- c(existing_pairs, pair1, pair2)
      }
    }
    phenotypes_vect <- paste(phenotypes1, phenotypes2, sep = "")
    CART_threshold <- vector("list", length(phenotypes_vect))
    names(CART_threshold) <- phenotypes_vect
  }


  # calculate bags for fit_parameters
  if (fit_parameters == TRUE) {
    result_get_bags <- get_bags(pheno_mat = pheno_mat,
                                dist_mat_or_tree = result$dist_mat_or_tree,
                                bagging = bagging,
                                misslabel_mat = FALSE,
                                bag_size = bag_size,
                                no_rf = repeats,
                                max_per_bag = max_per_bag)

    bags <- result_get_bags$bags
  }

  if (fit_parameters == TRUE & random_forest == TRUE) {
    parameters_rf <- fit_parameters_rf(bags,
                                       predictors,
                                       pheno_mat,
                                       phenotypes,
                                       repeats)

    sampsize <- parameters_rf[[1]]
    mtry <- parameters_rf[[2]]
    ntree <- parameters_rf[[3]]
    maxnodes <- parameters_rf[[4]]
  }

  if (fit_parameters == TRUE & adaboost == TRUE) {
    bags_ada <- as.data.frame(bags)
    names(bags_ada) <- paste0("ada", 1:ncol(bags_ada))
    pheno_mat$pheno <- as.character(pheno_mat$pheno)

    parameters_ada <- fit_parameters_ada(bags_ada,
                                         predictors,
                                         pheno_mat,
                                         as.character(phenotypes),
                                         repeats)

    max_depth <- parameters_ada[1]
    n_estimators <- parameters_ada[2]
    learning_rate <- parameters_ada[3]
  }

  bags <- vector("list", no_rounds)
  for (i in 1:no_rounds) {
    # here misslabel the samples
    if (bagging == "random_walk") {
      result_misslabel <- misslabel(pheno_mat, phenotypes, misslabel_no)
    }
    if (bagging == "phylogenetic_walk") {
      result_misslabel <- misslabel_phylo_walk(pheno_mat, phenotypes, misslabel_no, result$dist_mat_or_tree)
    }
    pheno_mat_miss <- result_misslabel$pheno_mat_miss
    misslabel_mats_list[[i]] <- result_misslabel$misslabel_mat
    # generate bag for this misslabelled pheno_mat
    result_get_bags <- get_bags(pheno_mat = pheno_mat_miss,
                                dist_mat_or_tree = result$dist_mat_or_tree,
                                bagging = bagging,
                                misslabel_mat = result_misslabel$misslabel_mat,
                                bag_size = bag_size,
                                no_rf = 1,
                                max_per_bag = max_per_bag)

    bags[[i]] <- result_get_bags$bags[[1]]
    curr_bag_order <- result_get_bags$bags[[1]]
    # construct bin_mat and pheno_mat for this current run
    index_select <- c()
    for (j in 1:length(curr_bag_order)) {
      the_row <- which(rownames(predictors) == curr_bag_order[j])
      index_select <- c(index_select, the_row)
    }

    # save bin_mat
    predictors_miss <- predictors[index_select, ]
    predictors_miss <- predictors_miss[,colnames(predictors_miss) != "ids"]
    bin_mats_list[[i]] <- predictors_miss

    # save pheno_mat
    pheno_mat_miss_bag <- pheno_mat_miss[index_select, ]
    pheno_mats_list[[i]] <- pheno_mat_miss_bag
  }

  df_threshold_rf <- data.frame()
  df_threshold_ada <- data.frame()
  df_threshold_log_reg <- data.frame()

  # initiate time vec
  time_diffs <- c()
  for (i in 1:no_rounds) {
    # make the initial timestamp
    initial_t <- Sys.time()

    if (CART == TRUE) {
      # CART builds own pheno_mat
      CART_result <- lapply(phenotypes_vect, run_CART_miss,
                            bin_mat = predictors,
                            pheno_mat = pheno_mat,
                            dist_mat_or_tree = result$dist_mat_or_tree,
                            bagging = bagging,
                            misslabel_no = misslabel_no,
                            bag_size = bag_size,
                            max_per_bag = max_per_bag)

      for (j in 1:length(CART_threshold)) {
        CART_threshold[[j]] <- rbind(CART_threshold[[j]], CART_result[[j]])
      }
      # calculate expected time to finish

      cat("Threshold Calculation Phase:", i, "CART model build out of", no_rounds, "\n")
    }

    if (random_forest == TRUE) {
      rf_result <- run_rf_miss(bin_mats_list[[i]],
                               pheno_mats_list[[i]],
                               predictors_test,
                               pheno_mat_test,
                               sampsize,
                               mtry,
                               ntree,
                               maxnodes)

      cat("Treshold calculation:", i, "random forest model build out of", no_rounds, "\n")
      rf_result1 <- postprocess_fit_miss(rf_result,
                                         misslabel_mats_list[[i]],
                                         phenotypes)

      df_threshold_rf <- rbind(df_threshold_rf, rf_result1)
    }

    if (adaboost == TRUE) {
      bin_mat_ada <- as.data.frame(bin_mats_list[[i]])
      pheno_mat_ada <- as.data.frame(pheno_mats_list[[i]])
      bin_mat_test_ada <- as.data.frame(predictors_test)
      pheno_mat_test_ada <- data.frame(pheno = pheno_mat_test)

      ada_result <- run_ada_miss(bin_mat_ada,
                                 pheno_mat_ada,
                                 bin_mat_test_ada,
                                 pheno_mat_test_ada,
                                 max_depth,
                                 n_estimators,
                                 learning_rate)

      rownames(ada_result) <- rownames(predictors_test)
      colnames(ada_result) <-  phenotypes
      ada_result$observed <- as.factor(pheno_mat_test_ada$pheno)
      cat("Threshold Calculation Phase:", i, "adaboost model build out of", no_rounds, "\n")

      ada_result1 <- postprocess_fit_miss(ada_result,
                                          misslabel_mats_list[[i]],
                                          phenotypes)

      df_threshold_ada <- rbind(df_threshold_ada, ada_result1)
    }

    if (ovr_log_reg == TRUE) {
      bin_mat_log_reg <- as.data.frame(bin_mats_list[[i]])
      pheno_mat_log_reg <- as.data.frame(pheno_mats_list[[i]])
      bin_mat_test_log_reg <- as.data.frame(predictors_test)
      pheno_mat_test_log_reg <- data.frame(pheno = pheno_mat_test)

      log_reg_result <- run_log_reg_miss(bin_mat_log_reg,
                                         pheno_mat_log_reg,
                                         bin_mat_test_log_reg,
                                         C_val)

      rownames(log_reg_result) <- rownames(predictors_test)
      colnames(log_reg_result) <-  phenotypes
      log_reg_result$observed <- as.factor(pheno_mat_test_log_reg$pheno)
      cat("Threshold Calculation Phase:", i, "log regression model build out of", no_rounds, "\n")

      log_reg_result1 <- postprocess_fit_miss(log_reg_result,
                                              misslabel_mats_list[[i]],
                                              phenotypes)

      df_threshold_log_reg <- rbind(df_threshold_log_reg, log_reg_result1)
    }
    elaps_time <- Sys.time()
    time_diffs <- c(as.numeric(elaps_time) - as.numeric(initial_t), time_diffs)
    time_secs <- median(time_diffs) * (no_rounds-i)
    hours <- floor(time_secs / 3600)  # Calculate hours
    minutes <- floor((time_secs %% 3600) / 60)
    cat("The estimated time to finish the Threshold Calculation Phase is", hours, "hours and", minutes, "minutes.\n")
  }

  # here you have to make sure that a strain misslabelled from pheno1 to pheno2
  # is distinguishible from strains is pheno2. The loops below use kolmogorov-smirnov
  # test to find out if the distribution pheno2 probs of non-misslabelled strain is
  # greater then distribution of pheno2 probs of misslabelled strains
  if (random_forest == TRUE) {
    colnames(df_threshold_rf) <- c(as.character(phenotypes), "pheno")
    ks_mat_rf <- matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
    for (i in 1:length(phenotypes)) {
      for (j in 1:length(phenotypes)) {
        if (i==j) {next}
        # i is row and j is column. P-val 1, 2 is a p val of strains in pheno1 vs
        # strains of pheno2 misslabelled into pheno1
        # get the origin_probs
        orig_probs <- as.numeric(df_threshold_rf[df_threshold_rf$pheno == i, i])
        misslabel_label <- paste0(j,"to",i)
        misslabel_probs <- as.numeric(df_threshold_rf[df_threshold_rf$pheno == misslabel_label, i])
        # calculate kolmogorov-smirnov test
        # null hypothesis of this test is that misslabel_probs is not greater then
        # orig_probs. If the p-value is small it means we can distinguish smaller
        # distribution misslabel_probs from greater distribution orig_probs
        ks_result <- suppressWarnings(stats::ks.test(misslabel_probs,
                                                     orig_probs,
                                                     alternative = "greater"))
        ks_mat_rf[i,j] <- ks_result$p.value
      }
    }
    rownames(ks_mat_rf) <- colnames(ks_mat_rf) <- lookup_tbl
  }

  if (adaboost == TRUE) {
    colnames(df_threshold_ada) <- c(as.character(phenotypes), "pheno")
    ks_mat_ada <- matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
    for (i in 1:length(phenotypes)) {
      for (j in 1:length(phenotypes)) {
        if (i==j) {next}
        orig_probs <- as.numeric(df_threshold_ada[df_threshold_ada$pheno == i, i])
        misslabel_label <- paste0(j,"to",i)
        misslabel_probs <- as.numeric(df_threshold_ada[df_threshold_ada$pheno == misslabel_label, i])
        ks_result <- suppressWarnings(stats::ks.test(misslabel_probs,
                                                     orig_probs,
                                                     alternative = "greater"))
        ks_mat_ada[i,j] <- ks_result$p.value
      }
    }
    rownames(ks_mat_ada) <- colnames(ks_mat_ada) <- lookup_tbl
  }


  if (ovr_log_reg == TRUE) {
    colnames(df_threshold_log_reg) <- c(as.character(phenotypes), "pheno")
    ks_mat_log_reg <- matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
    for (i in 1:length(phenotypes)) {
      for (j in 1:length(phenotypes)) {
        if (i==j) {next}
        orig_probs <- as.numeric(df_threshold_log_reg[df_threshold_log_reg$pheno == i, i])
        misslabel_label <- paste0(j,"to",i)
        misslabel_probs <- as.numeric(df_threshold_log_reg[df_threshold_log_reg$pheno == misslabel_label, i])
        ks_result <- suppressWarnings(stats::ks.test(misslabel_probs,
                                                     orig_probs,
                                                     alternative = "greater"))
        ks_mat_log_reg[i,j] <- ks_result$p.value
      }
    }
    rownames(ks_mat_log_reg) <- colnames(ks_mat_log_reg) <- lookup_tbl
  }

  if (CART == TRUE) {
    ks_mat_CART <- matrix(NA, nrow = length(phenotypes), ncol = length(phenotypes))
    for (i in 1:length(CART_threshold)) {
      current_mat <- CART_threshold[[i]]
      pheno1 <- colnames(current_mat)[1]
      pheno2 <- colnames(current_mat)[2]
      orig_probs <- as.numeric(current_mat[current_mat$observed == pheno1, pheno1])
      misslabel_label <- paste0(pheno2, "to", pheno1)
      misslabel_probs <- as.numeric(current_mat[current_mat$observed == misslabel_label, pheno1])
      ks_result <- suppressWarnings(stats::ks.test(misslabel_probs,
                                                   orig_probs,
                                                   alternative = "greater"))
      ks_mat_CART[as.numeric(pheno1), as.numeric(pheno2)] <- ks_result$p.value

      orig_probs <- as.numeric(current_mat[current_mat$observed == pheno2, pheno2])
      misslabel_label <- paste0(pheno1, "to", pheno2)
      misslabel_probs <- as.numeric(current_mat[current_mat$observed == misslabel_label, pheno2])
      ks_result <- suppressWarnings(stats::ks.test(misslabel_probs,
                                                   orig_probs,
                                                   alternative = "greater"))
      ks_mat_CART[as.numeric(pheno2), as.numeric(pheno1)] <- ks_result$p.value
    }
    rownames(ks_mat_CART) <- colnames(ks_mat_CART) <- lookup_tbl
  }

  ################################################################################
  # outlier calculation

  # get new bags
  if (random_forest | adaboost | ovr_log_reg) {
    result_get_bags <- get_bags(pheno_mat,
                                result$dist_mat_or_tree,
                                bagging,
                                misslabel_mat = FALSE,
                                bag_size,
                                no_rounds,
                                max_per_bag)

    bags <- result_get_bags$bags
  }

  # it is necessary to give the features this generic names because otherwise
  # the randomForest function is not gonna work
  names(bin_mat)[1:ncol(bin_mat)] <- paste0("V", 1:ncol(bin_mat))

  # convert pheno to factor so RF does not run regression
  pheno_mat$pheno <- as.factor(pheno_mat$pheno)

  # now convert the pheno_mat and bin_mat so they can be used as a evaluation of the models
  predictors <- bin_mat
  pheno <- pheno_mat$pheno[match(rownames(predictors), pheno_mat$ids)]
  names(pheno) <- rownames(predictors)


  # prepare objects that will hold the results
  if (random_forest == TRUE) {
    # will hold the AUC for each iteration
    aucs_rf <- rep(NA, no_rounds)

    # also create a list that will store all the prediction prob for the strains
    # the list contains all the probs for all phenotypes in format
    # probs[[strain_ID]][df with all probs]
    probs_rf <- setNames(replicate(nrow(bin_mat),data.frame()), rownames(bin_mat))

    # two dataframes that hold feature accuracy and decrease in gini
    df_gini_rf <- data.frame(gfs = key_gfsnames)
    df_accuracy_rf <- data.frame(gfs = key_gfsnames)

    # list that will hold all the mats for outlier visualisation
    list_proxy_mats_median <- vector("list", no_rounds)
    # mat that will do the same but for the sum not median
    final_proxy_mat_sum <- matrix(0, nrow = nrow(predictors), ncol = nrow(predictors))

  }

  if (adaboost == TRUE) {
    aucs_ada <- rep(NA, no_rounds)
    probs_ada <- setNames(replicate(nrow(bin_mat),data.frame()), rownames(bin_mat))
    df_gini_ada <- data.frame(gfs = key_gfsnames)
  }

  if (ovr_log_reg == TRUE) {
    aucs_log_reg <- rep(NA, no_rounds)
    probs_log_reg <- setNames(replicate(nrow(bin_mat),data.frame()), rownames(bin_mat))
    df_coef_log_reg <- vector("list", no_rounds)
  }

  if (CART == TRUE) {
    aucs_CART  <- matrix(NA, ncol = length(phenotypes_vect), nrow = no_rounds)
    colnames(aucs_CART) <- phenotypes_vect

    splits_CART  <- matrix(NA, ncol = length(phenotypes_vect), nrow = no_rounds)
    colnames(splits_CART) <- phenotypes_vect

    probs_CART <- vector("list", length(phenotypes_vect))
    names(probs_CART) <- phenotypes_vect
    for (i in 1:length(probs_CART)) {
      phenotype1 <- substr(phenotypes_vect[i], 1,1)
      phenotype2 <- substr(phenotypes_vect[i], 2,2)
      pheno_mat_tmp <- pheno_mat[pheno_mat$pheno == phenotype1 |
                                   pheno_mat$pheno == phenotype2, ]
      probs_CART[[i]] <- vector("list", length(pheno_mat_tmp$ids))
      names(probs_CART[[i]]) <- pheno_mat_tmp$ids
      for (j in 1:length(pheno_mat_tmp$ids)) {
        probs_CART[[i]][[j]] <- data.frame()
      }
    }

    importances_CART <- vector("list", length(phenotypes_vect))
    names(importances_CART) <- phenotypes_vect
    for (i in 1:length(importances_CART)) {
      importances_CART[[i]] <- data.frame(gfs = key_gfsnames)
    }
    # list that will hold all the proximity mats for CART
    list_proxy_mats_CART <- vector("list", no_rounds)
  }

  # initiate vec that will hold the times
  time_diffs <- c()

  for (i in 1:no_rounds) {
    # Make the initial timestamp
    initial_t <- Sys.time()

    # first build a bin_mat for the models
    index_select <- c()
    for (j in 1:length(bags[[i]])) {
      the_row <- which(rownames(bin_mat) == bags[[i]][j])
      index_select <- c(index_select, the_row)
    }
    bin_mat_model <- bin_mat[index_select, ]

    # now the same with pheno_mat
    index_select <- c()
    for (j in 1:length(bags[[i]])) {
      index_row <- which(pheno_mat$ids == bags[[i]][j])
      index_select <- c(index_select, index_row)
    }
    pheno_mat_model <- pheno_mat[index_select, ]

    if (adaboost == TRUE) {
      # adaboost is run by python script

      result_ada <- run_ada(bin_mat_model,
                            pheno_mat_model,
                            predictors,
                            data.frame(pheno = pheno),
                            max_depth,
                            n_estimators,
                            learning_rate)
      # sometimes it can happen that a Runtime warning is returned from function
      # run_ada. This is because a divide-by-zero error is encountered during
      # feature importances calculation. As a results there might be a few features
      # that have importance == NaN. Set their value to 0
      result_ada[[2]]$MeanDecreaseGini[
        which(is.na(result_ada[[2]]$MeanDecreaseGini))] <- 0

      cat("Outlier Calculation Phase:", i, "adaboost model build out of", no_rounds, "\n")
      predictions <- result_ada[[1]]
      predictions$observed <- as.factor(predictions$observed)
      rownames(predictions) <- rownames(predictors)
      colnames(predictions) <-  c(phenotypes, "observed")
      result_ada1 <- postprocess_fit(predictions,
                                     probs_ada,
                                     result_ada[[2]],
                                     phenotypes = phenotypes,
                                     ada = TRUE)

      df_gini_ada <- cbind(df_gini_ada, result_ada1$importance_gini)
      aucs_ada[i] <- result_ada1$auc
      probs_ada <- result_ada1$probs

    }

    if (ovr_log_reg == TRUE) {
      # log regression is run by python script
      result_log_reg <- run_log_reg(bin_mat_model,
                                    pheno_mat_model,
                                    predictors,
                                    data.frame(pheno = pheno),
                                    C_val)

      df_coef_log_reg[[i]] <- result_log_reg[[2]]

      cat("Outlier Calculation Phase:", i, "log regression model build out of", no_rounds, "\n")
      predictions <- result_log_reg[[1]]
      predictions$observed <- as.factor(predictions$observed)
      rownames(predictions) <- rownames(predictors)
      colnames(predictions) <-  c(phenotypes, "observed")
      result_log_reg1 <- postprocess_fit(predictions,
                                        probs_log_reg,
                                        phenotypes = phenotypes)

      aucs_log_reg[i] <- result_log_reg1$auc
      probs_log_reg <- result_log_reg1$probs
    }

    if (random_forest == TRUE) {
      bin_mat_model$ids <- rownames(bin_mat_model)
      bin_mat_model <- merge(bin_mat_model, pheno_mat_model, by = "ids")

      key_ids <- bin_mat_model$ids
      bin_mat_model <- bin_mat_model[,-1]

      result_rf <- run_rf(bin_mat_model,
                          predictors,
                          pheno,
                          sampsize,
                          mtry,
                          ntree,
                          maxnodes,
                          final_proxy_mat_sum)

      list_proxy_mats_median[[i]] <- result_rf$dist_mat
      final_proxy_mat_sum <- result_rf$final_proxy_mat_sum

      cat("Outlier Calculation Phase:", i, "random forest model build out of", no_rounds, "\n")

      result_rf <- postprocess_fit(result_rf$predictions,
                                   probs_rf,
                                   result_rf$df_importance,
                                   phenotypes = phenotypes,
                                   rf = TRUE)
      df_accuracy_rf <- cbind(df_accuracy_rf, result_rf$importance_accuracy)
      df_gini_rf <- cbind(df_gini_rf, result_rf$importance_gini)
      aucs_rf[i] <- result_rf$auc
      probs_rf <- result_rf$probs
    }

    if (CART == TRUE) {
      # CART again uses their own bags
      CART_result <- lapply(phenotypes_vect, run_CART,
                            bin_mat = bin_mat,
                            pheno_mat = pheno_mat,
                            dist_mat_or_tree = result$dist_mat_or_tree,
                            bagging = bagging,
                            bag_size = bag_size,
                            max_per_bag = max_per_bag)
      # name a mat that will hold the final proxies
      final_proxy_mat <- matrix(NA, nrow = nrow(pheno_mat), ncol = nrow(pheno_mat))
      colnames(final_proxy_mat) <- rownames(final_proxy_mat) <- pheno_mat$ids

      for (j in 1:length(phenotypes_vect)) {
        res_run <- CART_result[[j]]

        df_tmp <- probs_CART[[which(names(probs_CART) == phenotypes_vect[j])]]
        for (k in 1:nrow(res_run$probs)) {
          # adds probs from CART tree to list probs_CART
          index <- which(names(df_tmp) == res_run$probs$ids[k])
          probs_CART[[which(names(probs_CART) ==  phenotypes_vect[j])]][[index]] <-
            rbind(df_tmp[[index]], res_run$probs[k, c(1,2)])
        }
        # add importances to the list
        importances_CART[[which(names(importances_CART) == phenotypes_vect[j])]] <-
          cbind(importances_CART[[which(names(importances_CART) == phenotypes_vect[j])]], res_run$importance)
        # add the number of splits
        index <- which(colnames(splits_CART) == phenotypes_vect[j])
        splits_CART[which(is.na(splits_CART[, index]))[1], index] <- res_run$splits
        # add the AUCs from this run
        index <- which(colnames(aucs_CART) == phenotypes_vect[j])
        aucs_CART[which(is.na(aucs_CART[, index]))[1], index] <- res_run$curr_auc
        # construct the final proxy mat
        index <- which(is.element(rownames(final_proxy_mat), rownames(res_run$heatmat)))
        final_proxy_mat[index,index] <- res_run$heatmat
      }
      list_proxy_mats_CART[[i]] <- final_proxy_mat
      cat("Outlier Calculation Phase:", i, "CART model build out of", no_rounds, "\n")
    }

    elaps_time <- Sys.time()
    time_diffs <- c(as.numeric(elaps_time) - as.numeric(initial_t), time_diffs)
    time_secs <- median(time_diffs) * (no_rounds-i)
    hours <- floor(time_secs / 3600)
    minutes <- floor((time_secs %% 3600) / 60)
    cat("The estimated time to finish the Outlier Calculation Phase is", hours, "hours and", minutes, "minutes.\n")
  }
  ################################################################################
  # Analyse results from the model part
  # the original bin_mat is necessary for the calculation
  bin_mat_original <- as.data.frame(sapply(as.data.frame(bin_mat_copy),
                                           as.numeric))
  rownames(bin_mat_original) <- rownames(bin_mat_copy)


  if (random_forest == TRUE) {
    # use kolmogorov smirnov test to find out if some strains were allochthonous
    result_probs <- lapply(probs_rf,
                           kol_smir_test,
                           df_threshold = df_threshold_rf,
                           pheno_mat = pheno_mat,
                           phenotypes = phenotypes,
                           lookup_tbl = lookup_tbl)
    result_probs <- t(as.data.frame(result_probs))
    result_probs <- cbind(data.frame(strain_ids = rownames(result_probs)), result_probs)
    colnames(result_probs) <- c("strain_ids", lookup_tbl, "predicted", "observed")
    result_probs_rf <- result_probs
    if (write_data == TRUE) {
      dir <- paste(save_dir, "aurora_results_rf_", Sys.Date(), ".csv", sep = "")
      readr::write_csv(result_probs_rf, dir)
    }

    # create a plot for the results
    for (i in 1:length(lookup_tbl)) {
      result_probs_tmp <- result_probs[,colnames(result_probs) == lookup_tbl[i]]
      index <- which(grepl("not typical", result_probs_tmp))
      result_probs$predicted[grepl("not typical", result_probs_tmp)] <- paste0(lookup_tbl[i], " not typical")
    }

    grid <- as.data.frame(expand.grid(unique(result_probs$predicted), lookup_tbl))
    df_plot_long <- cbind(grid, rep(0, nrow(grid)))
    for (i in 1:length(lookup_tbl)) {
      tbl_predicted <- table(result_probs$predicted[result_probs$observed == lookup_tbl[i]])
      tbl_predicted <- as.data.frame(tbl_predicted)
      indexes <- match(tbl_predicted$Var1, df_plot_long[df_plot_long[,2] == lookup_tbl[i], 1])
      df_plot_long[which(df_plot_long[,2] == lookup_tbl[i])[indexes], 3] <- tbl_predicted$Freq
    }
    colnames(df_plot_long) <- c("col2", "col1", "col3")
    df_plot_long$col2 <- as.character(df_plot_long$col2)
    df_plot_long$col1 <- as.character(df_plot_long$col1)
    for (i in 1:nrow(df_plot_long)) {
      if (df_plot_long$col2[i] == df_plot_long$col1[i]) {
        df_plot_long$col2[i] <- "autochthonous"
      }

      if (grepl("not typical", df_plot_long$col2[i])) {
        df_plot_long$col2[i] <- "weak autochthonous"
      }
    }

    df_plot_long <- df_plot_long[df_plot_long$col3 > 0, ]

    plot_rf <- ggplot2::ggplot(df_plot_long, aes(fill = col2, y = col3, x=col1)) +
      geom_bar(position="fill", stat="identity") +
      ggtitle("aurora: random forest result") +
      xlab("") + ylab("") +
      theme(text = element_text(size = 20)) +
      guides(fill=guide_legend(title=""))

    # now analyze the features
    df_gini_rf <- split_gfs(df_gini_rf, pangenome, DRAM)
    df_accuracy_rf <- split_gfs(df_accuracy_rf, pangenome, DRAM)

    top_features <- function(x, no_gfs) {
      x <- x[x > 0]
      if (length(x) > 0) {
        return(stats::median(x))
      } else {return(0)}
    }

    df_gini_median <- apply(df_gini_rf[,2:ncol(df_gini_rf)], 1, top_features, no_gfs = no_gfs)
    df_accuracy_median <- apply(df_accuracy_rf[,2:ncol(df_accuracy_rf)], 1, top_features, no_gfs = no_gfs)

    df_final <- data.frame(gf_ID = df_gini_rf[, 1],
                           MeanDecreaseGini_median = df_gini_median,
                           MeanDecreaseAccuracy_median = df_accuracy_median)
    rownames(df_final) <- df_gini_rf[, 1]
    # now create columns that show how the gfs are distributed in classes
    # Here we have to work with the original bin_mat
    # Calculate fraction of strains with the GF and without for each class
    calc_ratio <- function(x) {
      paste0(sum(x), "|", length(x))
    }

    df_ratios <- as.data.frame(matrix(NA, nrow = nrow(df_final),
                                      ncol = length(phenotypes)))
    for (i in 1:length(phenotypes)) {
      indexes <- pheno_mat$ids[pheno_mat$pheno == phenotypes[i]]
      bin_mat_pheno <- bin_mat_original[is.element(rownames(bin_mat_original), indexes), ]
      fin_vect <- apply(bin_mat_pheno ,2 ,calc_ratio)
      df_ratios[,i] <- fin_vect[match(rownames(df_final), names(fin_vect))]
      colnames(df_ratios)[i] <- paste0(as.character(lookup_tbl[i]), " ratio")
    }
    df_final <- cbind(df_final, df_ratios)
    df_final <- df_final[order(df_final$MeanDecreaseGini_median, decreasing = TRUE), ]
    df_final_gfs_rf <- df_final
    if (write_data == TRUE) {
      dir <- paste(save_dir, "classification_variants_random_forest_", Sys.Date(), ".csv", sep = "")
      readr::write_csv(df_final_gfs_rf, dir)
    }
  }

  if (adaboost == TRUE) {
    result_probs <- lapply(probs_ada,
                           kol_smir_test,
                           df_threshold = df_threshold_ada,
                           pheno_mat = pheno_mat,
                           phenotypes = phenotypes,
                           lookup_tbl = lookup_tbl)
    result_probs <- t(as.data.frame(result_probs))
    result_probs <- cbind(data.frame(strain_ids = rownames(result_probs)), result_probs)
    colnames(result_probs) <- c("strain_ids", lookup_tbl, "predicted", "observed")
    result_probs_ada <- result_probs
    if (write_data == TRUE) {
      dir <- paste(save_dir, "aurora_results_ada_", Sys.Date(), ".csv", sep = "")
      readr::write_csv(result_probs_ada, dir)
    }

    # create a plot for the results
    for (i in 1:length(lookup_tbl)) {
      result_probs_tmp <- result_probs[,colnames(result_probs) == lookup_tbl[i]]
      index <- which(grepl("not typical", result_probs_tmp))
      result_probs$predicted[grepl("not typical", result_probs_tmp)] <- paste0(lookup_tbl[i], " not typical")
    }

    grid <- as.data.frame(expand.grid(unique(result_probs$predicted), lookup_tbl))
    df_plot_long <- cbind(grid, rep(0, nrow(grid)))
    for (i in 1:length(lookup_tbl)) {
      tbl_predicted <- table(result_probs$predicted[result_probs$observed == lookup_tbl[i]])
      tbl_predicted <- as.data.frame(tbl_predicted)
      indexes <- match(tbl_predicted$Var1, df_plot_long[df_plot_long[,2] == lookup_tbl[i], 1])
      df_plot_long[which(df_plot_long[,2] == lookup_tbl[i])[indexes], 3] <- tbl_predicted$Freq
    }
    colnames(df_plot_long) <- c("col2", "col1", "col3")
    df_plot_long$col2 <- as.character(df_plot_long$col2)
    df_plot_long$col1 <- as.character(df_plot_long$col1)
    for (i in 1:nrow(df_plot_long)) {
      if (df_plot_long$col2[i] == df_plot_long$col1[i]) {
        df_plot_long$col2[i] <- "autochthonous"
      }

      if (grepl("not typical", df_plot_long$col2[i])) {
        df_plot_long$col2[i] <- "weak autochthonous"
      }
    }

    df_plot_long <- df_plot_long[df_plot_long$col3 > 0, ]

    plot_ada <- ggplot2::ggplot(df_plot_long, aes(fill = col2, y = col3, x=col1)) +
      geom_bar(position="fill", stat="identity") +
      ggtitle("aurora: adaboost result") +
      xlab("") + ylab("") +
      theme(text = element_text(size = 20)) +
      guides(fill=guide_legend(title=""))

    # now analyze the features
    df_gini_ada <- split_gfs(df_gini_ada, pangenome, DRAM)

    top_features <- function(x, no_gfs) {
      x <- x[x > 0]
      if (length(x) > 0) {
        return(stats::median(x))
      } else {return(0)}
    }

    df_gini_median <- apply(df_gini_ada[,2:ncol(df_gini_ada)], 1, top_features, no_gfs = no_gfs)
    df_final <- data.frame(gf_ID = df_gini_ada[, 1],
                           MeanDecreaseGini_median = df_gini_median)
    rownames(df_final) <- df_gini_ada[, 1]

    calc_ratio <- function(x) {
      paste0(sum(x), "|", length(x))
    }

    df_ratios <- as.data.frame(matrix(NA, nrow = nrow(df_final),
                                      ncol = length(phenotypes)))
    for (i in 1:length(phenotypes)) {
      indexes <- pheno_mat$ids[pheno_mat$pheno == phenotypes[i]]
      bin_mat_pheno <- bin_mat_original[is.element(rownames(bin_mat_original), indexes), ]
      fin_vect <- apply(bin_mat_pheno ,2 ,calc_ratio)
      df_ratios[,i] <- fin_vect[match(rownames(df_final), names(fin_vect))]
      colnames(df_ratios)[i] <- paste0(as.character(lookup_tbl[i]), " ratio")
    }
    df_final <- cbind(df_final, df_ratios)
    df_final <- df_final[order(df_final$MeanDecreaseGini_median, decreasing = TRUE), ]
    df_final_gfs_ada <- df_final

    if (write_data == TRUE) {
      dir <- paste(save_dir, "classification_variants_adaboost_", Sys.Date(), ".csv", sep = "")
      readr::write_csv(df_final_gfs_ada, dir)
    }
  }

  # analyse results from log regression
  if (ovr_log_reg == TRUE) {
    result_probs <- lapply(probs_log_reg,
                           kol_smir_test,
                           df_threshold = df_threshold_log_reg,
                           pheno_mat = pheno_mat,
                           phenotypes = phenotypes,
                           lookup_tbl = lookup_tbl)
    result_probs <- t(as.data.frame(result_probs))
    result_probs_log_reg <- cbind(data.frame(strain_ids = rownames(result_probs)), result_probs)
    colnames(result_probs_log_reg) <- c("strain_ids", lookup_tbl, "predicted", "observed")

    if (write_data == TRUE) {
      dir <- paste(save_dir, "aurora_results_log_regression_", Sys.Date(), ".csv", sep = "")
      readr::write_csv(result_probs_log_reg, dir)
    }

    # create a plot for the results
    for (i in 1:length(lookup_tbl)) {
      result_probs_tmp <- result_probs_log_reg[,colnames(result_probs_log_reg) == lookup_tbl[i]]
      index <- which(grepl("not typical", result_probs_tmp))
      result_probs_log_reg$predicted[grepl("not typical", result_probs_tmp)] <- paste0(lookup_tbl[i], " not typical")
    }

    grid <- as.data.frame(expand.grid(unique(result_probs_log_reg$predicted), lookup_tbl))
    df_plot_long <- cbind(grid, rep(0, nrow(grid)))
    for (i in 1:length(lookup_tbl)) {
      tbl_predicted <- table(result_probs_log_reg$predicted[result_probs_log_reg$observed == lookup_tbl[i]])
      tbl_predicted <- as.data.frame(tbl_predicted)
      indexes <- match(tbl_predicted$Var1, df_plot_long[df_plot_long[,2] == lookup_tbl[i], 1])
      df_plot_long[which(df_plot_long[,2] == lookup_tbl[i])[indexes], 3] <- tbl_predicted$Freq
    }
    colnames(df_plot_long) <- c("col2", "col1", "col3")
    df_plot_long$col2 <- as.character(df_plot_long$col2)
    df_plot_long$col1 <- as.character(df_plot_long$col1)
    for (i in 1:nrow(df_plot_long)) {
      if (df_plot_long$col2[i] == df_plot_long$col1[i]) {
        df_plot_long$col2[i] <- "autochthonous"
      }

      if (grepl("not typical", df_plot_long$col2[i])) {
        df_plot_long$col2[i] <- "weak autochthonous"
      }
    }

    df_plot_long <- df_plot_long[df_plot_long$col3 > 0, ]

    plot_log_reg <- ggplot2::ggplot(df_plot_long, aes(fill = col2, y = col3, x=col1)) +
      geom_bar(position="fill", stat="identity") +
      ggtitle("aurora: log regression result") +
      xlab("") + ylab("") +
      theme(text = element_text(size = 20)) +
      guides(fill=guide_legend(title=""))

    # now analyze the features
    get_med_log <- function(x) {
      x <- x[x != 0]
      if (length(x) == 0) {return(0)}
      return(stats::median(x))
    }
    # matrix that will hold the final medians of the log reg coeffs
    coefs_feature_final <- matrix(NA, nrow = nrow(df_coef_log_reg[[1]]), ncol = length(phenotypes))
    for (i in 1:nrow(df_coef_log_reg[[1]])) {
      coefs_feature <- matrix(NA, nrow = no_rounds, ncol = length(phenotypes))
      for (j in 1:length(df_coef_log_reg)) {
        coefs_feature[j,] <- as.numeric(df_coef_log_reg[[j]][i,])
      }
      coefs_feature_final[i, ] <- apply(coefs_feature, 2, get_med_log)
    }
    coefs_feature_final <- cbind(data.frame(gfs = key_gfsnames),
                                 as.data.frame(coefs_feature_final))

    coefs_feature_final <- split_gfs(coefs_feature_final, pangenome, DRAM)
    colnames(coefs_feature_final) <- c("gfs", lookup_tbl)
    # now create columns that show how the gfs are distributed in classes
    # Here we have to work with the original bin_mat
    # Calculate fraction of strains with the GF and without for each class
    calc_ratio <- function(x) {
      paste0(sum(x), "|", length(x))
    }

    df_ratios <- as.data.frame(matrix(NA, nrow = nrow(coefs_feature_final),
                                      ncol = length(phenotypes)))
    for (i in 1:length(phenotypes)) {
      indexes <- pheno_mat$ids[pheno_mat$pheno == phenotypes[i]]
      bin_mat_pheno <- bin_mat_original[is.element(rownames(bin_mat_original), indexes), ]
      fin_vect <- apply(bin_mat_pheno ,2 ,calc_ratio)
      df_ratios[,i] <- fin_vect[match(coefs_feature_final$gfs, names(fin_vect))]
      colnames(df_ratios)[i] <- paste0(as.character(lookup_tbl[i]), " ratio")
    }
    df_final <- cbind(coefs_feature_final, df_ratios)
    col1 <- 2:(ncol(coefs_feature_final))
    col2 <- (ncol(coefs_feature_final)+1):(ncol(df_ratios)+ncol(coefs_feature_final))
    index <- c()
    for (i in 1:length(col1)) {
      index <- c(index, col2[i], col1[i])
    }
    df_final_log_reg <- df_final[,c(1,index)]
    if (write_data == TRUE) {
      dir <- paste(save_dir, "classification_variants_log_reg_", Sys.Date(), ".csv", sep = "")
      readr::write_csv(df_final_log_reg, dir)
    }
  }

  if (CART == TRUE) {
    # CART needs a special way to deal with probs and importances
    vect_CART <- vector("list", length(phenotypes_vect))
    names(vect_CART) <- phenotypes_vect
    for (i in 1:length(phenotypes_vect)) {
      phenotype1 <- substr(phenotypes_vect[i], 1,1)
      phenotype2 <- substr(phenotypes_vect[i], 2,2)
      probs_curr <- probs_CART[[i]]
      threshold_curr <- CART_threshold[[i]]
      colnames(threshold_curr)[3] <- "pheno"
      pheno_mat_curr <- pheno_mat[pheno_mat$pheno == phenotype1 |
                                    pheno_mat$pheno == phenotype2, ]
      lookup_tbl_curr <- lookup_tbl[as.numeric(c(phenotype1, phenotype2))]
      result_probs <- lapply(probs_curr,
                             kol_smir_test,
                             df_threshold = threshold_curr,
                             pheno_mat = pheno_mat_curr,
                             phenotypes = as.factor(c(phenotype1, phenotype2)),
                             lookup_tbl = lookup_tbl_curr)
      vect_CART[[i]] <- cbind(data.frame(strain_ids = rownames(t(as.data.frame(result_probs)))),
                              t(as.data.frame(result_probs)))
    }

    df_final <- as.data.frame(matrix(NA, ncol = length(phenotypes)+3,
                                     nrow = nrow(pheno_mat)))
    df_final[,1] <- pheno_mat$ids
    df_final[,ncol(df_final)] <- lookup_tbl[match(pheno_mat$pheno, names(lookup_tbl))]

    for (i in 1:nrow(pheno_mat)) {
      obsrv_pheno <- pheno_mat$pheno[i]
      # get indexes of all mats where the phenotype is present
      index <- which(grepl(obsrv_pheno, phenotypes_vect))
      result_label <- c()
      result_vect <- setNames(rep(NA, length(phenotypes)), phenotypes)
      for (j in 1:length(index)) {
        df_tmp <- vect_CART[[index[j]]]
        phenos <- str_split(names(vect_CART)[index[j]],"")[[1]]
        othr_pheno <- which(phenos != obsrv_pheno)
        othr_pheno_probs <- df_tmp[df_tmp$strain_ids == pheno_mat$ids[i],
                                   c(2,3)[othr_pheno]]
        result_vect[names(result_vect) == phenos[phenos != obsrv_pheno]] <- othr_pheno_probs
        result_label <- c(result_label, df_tmp[df_tmp$strain_ids == pheno_mat$ids[i], 4])
      }
      pheno <- unique(str_split(paste(result_label, collapse = " or "), pattern = " or ")[[1]])
      df_final[i,ncol(df_final)-1] <- paste(pheno, collapse = " or ")
      obsrv_pheno1 <- lookup_tbl[obsrv_pheno]
      if(grepl(obsrv_pheno1, paste(pheno, collapse = " or "))) {
        result_vect[obsrv_pheno] <- TRUE
      } else {
        result_vect[obsrv_pheno] <- FALSE
      }
      df_final[i, 2:(length(phenotypes)+1)] <- result_vect
    }
    colnames(df_final) <- c("strain_ids", lookup_tbl, "predicted", "observed")
    df_final_CART <- df_final

    if (write_data == TRUE) {
      dir <- paste(save_dir, "aurora_results_CART_", Sys.Date(), ".csv", sep = "")
      readr::write_csv(df_final_CART, dir)
    }


    # create a plot for the results
    for (i in 1:length(lookup_tbl)) {
      result_probs_tmp <- df_final_CART[,colnames(df_final_CART) == lookup_tbl[i]]
      index <- which(grepl("not typical", result_probs_tmp))
      df_final_CART$predicted[grepl("not typical", result_probs_tmp)] <- paste0(lookup_tbl[i], " not typical")
    }

    grid <- as.data.frame(expand.grid(unique(df_final_CART$predicted), lookup_tbl))
    df_plot_long <- cbind(grid, rep(0, nrow(grid)))
    for (i in 1:length(lookup_tbl)) {
      tbl_predicted <- table(df_final_CART$predicted[df_final_CART$observed == lookup_tbl[i]])
      tbl_predicted <- as.data.frame(tbl_predicted)
      indexes <- match(tbl_predicted$Var1, df_plot_long[df_plot_long[,2] == lookup_tbl[i], 1])
      df_plot_long[which(df_plot_long[,2] == lookup_tbl[i])[indexes], 3] <- tbl_predicted$Freq
    }
    colnames(df_plot_long) <- c("col2", "col1", "col3")
    df_plot_long$col2 <- as.character(df_plot_long$col2)
    df_plot_long$col1 <- as.character(df_plot_long$col1)
    for (i in 1:nrow(df_plot_long)) {
      if (df_plot_long$col2[i] == df_plot_long$col1[i]) {
        df_plot_long$col2[i] <- "autochthonous"
      }

      if (grepl("not typical", df_plot_long$col2[i])) {
        df_plot_long$col2[i] <- "weak autochthonous"
      }
    }
    df_plot_long <- df_plot_long[df_plot_long$col3 > 0, ]

    plot_CART <- ggplot2::ggplot(df_plot_long, aes(fill = col2, y = col3, x=col1)) +
      geom_bar(position="fill", stat="identity") +
      ggtitle("aurora: CART result") +
      xlab("") + ylab("") +
      theme(text = element_text(size = 20)) +
      guides(fill=guide_legend(title=""))

    # now analyse the features
    # gets median of all non-zero elements
    get_med <- function(x) {
      x <- as.numeric(x[-1])
      x <- x[x > 0]
      if (length(x) == 0) {
        return (0)
      } else {return (stats::median(x))}
    }

    for (i in 1:length(phenotypes_vect)) {
      importances_CART_curr <- importances_CART[[i]]
      medians <- apply(importances_CART_curr, 1 ,get_med)
      df_tmp <- data.frame(gfs = importances_CART_curr$gfs, y = medians)

      if (i == 1) {
        df_importances_CART <- data.frame(gfs = df_tmp$gfs)
      }
      df_importances_CART <- cbind(df_importances_CART, as.data.frame(df_tmp$y))

      phenotype1 <- substr(phenotypes_vect[i], 1,1)
      phenotype2 <- substr(phenotypes_vect[i], 2,2)
      col_nm <- paste(lookup_tbl[names(lookup_tbl) == phenotype1],"vs",
                      lookup_tbl[names(lookup_tbl) == phenotype2])
      colnames(df_importances_CART)[ncol(df_importances_CART)] <- col_nm
    }
    df_importances_CART <- split_gfs(df_importances_CART, pangenome, DRAM)
    # calcualate the presence/absence ratio for each gf
    calc_ratio <- function(x) {
      paste0(sum(x), "|", length(x))
    }

    df_ratios <- as.data.frame(matrix(NA, nrow = nrow(df_importances_CART),
                                      ncol = length(phenotypes)))
    for (i in 1:length(phenotypes)) {
      indexes <- pheno_mat$ids[pheno_mat$pheno == phenotypes[i]]
      bin_mat_pheno <- bin_mat_original[is.element(rownames(bin_mat_original), indexes), ]
      fin_vect <- apply(bin_mat_pheno ,2 ,calc_ratio)
      df_ratios[,i] <- fin_vect[match(df_importances_CART$gfs, names(fin_vect))]
      colnames(df_ratios)[i] <- paste0(as.character(lookup_tbl[i]), " ratio")
    }
    df_importances_CART <- cbind(df_importances_CART, df_ratios)
    if (write_data == TRUE) {
      dir <- paste(save_dir, "classification_variants_CART_", Sys.Date(), ".csv", sep = "")
      readr::write_csv(df_importances_CART, dir)
    }
    # get df with average number of splits and the aucs to show how hard it is
    # to differentiate classes
    df_one_vs_one <- as.data.frame(matrix(NA, nrow = length(phenotypes_vect), ncol = 5))
    colnames(df_one_vs_one) <- c("classes", "AUC", "number_of_splits", "min_splits", "max_splits")
    df_one_vs_one$classes <- colnames(df_importances_CART)[2:(length(phenotypes_vect)+1)]
    for (i in 1:length(phenotypes_vect)) {
      df_one_vs_one$AUC[i] <- stats::median(aucs_CART[,i])
      uniqv <- unique(splits_CART[,i])
      df_one_vs_one$number_of_splits[i] <- uniqv[which.max(tabulate(match(splits_CART[,i], uniqv)))]
      df_one_vs_one$min_splits[i] <- splits_CART[,i][which.min(splits_CART[,i])]
      df_one_vs_one$max_splits[i] <- splits_CART[,i][which.max(splits_CART[,i])]
    }

    if (write_data == TRUE) {
      dir <- paste(save_dir, "aurora_CART_paramteres_", Sys.Date(), ".csv", sep = "")
      readr::write_csv(df_one_vs_one, dir)
    }

  }

  ###############################################################################

  if (random_forest == TRUE) {
    # analyse the proximities from the random forest
    no_samples <- nrow(list_proxy_mats_median[[1]])
    final_proxy_mat <- matrix(NA, nrow = no_samples, ncol = no_samples )
    colnames(final_proxy_mat) <- colnames(list_proxy_mats_median[[1]])
    rownames(final_proxy_mat) <- colnames(list_proxy_mats_median[[1]])

    for (i in 1:no_samples) {
      for (j in 1:no_samples) {
        vect_dist <- c()
        for (k in 1:length(list_proxy_mats_median)) {
          vect_dist <- c(vect_dist, list_proxy_mats_median[[k]][j,i])
        }
        final_proxy_mat[j,i] <- stats::median(vect_dist)
      }
    }
    # plot the results
    if (plot_random_forest == TRUE & length(phenotypes) > 11 & write_data == TRUE) {
      plot_random_forest <- FALSE
      warning("The number of analysed phenotypes is higher then 11. Results of random forest will not be plotted. \n")
    }

    if (plot_random_forest == TRUE & write_data == TRUE) {
      # the maximum pheno is 11 because then we run out of colors
      auto_colors = c("darkred",
                      "blue4",
                      "green4",
                      "black",
                      "darkgoldenrod4",
                      "lightpink3",
                      "darkorange4",
                      "cyan3",
                      "deeppink4",
                      "antiquewhite4",
                      "aquamarine4")

      allo_colors = c("firebrick1",
                      "blue1",
                      "green1",
                      "gray20",
                      "darkgoldenrod1",
                      "lightpink1",
                      "darkorange1",
                      "cyan1",
                      "deeppink1",
                      "antiquewhite1",
                      "aquamarine1")
      col_color <- c()
      row_color <- c()
      for (i in 1:nrow(result_probs_rf)) {
        index <- which(lookup_tbl == result_probs_rf$observed[i])
        row_color <- c(row_color, auto_colors[index])
        if (result_probs_rf$predicted[i] == result_probs_rf$observed[i]) {
          if (result_probs_rf[i, 1+index] == "TRUE not typical") {
            col_color <- c(col_color, allo_colors[index])
          } else {
            col_color <- c(col_color, auto_colors[index])
          }
        } else {
          col_color <- c(col_color, "white")
        }
      }

      dir <- paste(save_dir, "aurora_heatmap_sum_random_forest_",
                   Sys.Date(), ".png", sep = "")

      dist_matrix <- as.dist(final_proxy_mat_sum)
      hc_result <- stats::hclust(dist_matrix, method = "ward.D")

      png(filename=dir)
      gplots::heatmap.2(final_proxy_mat_sum,
                ColSideColors = col_color, # predicted
                RowSideColors = row_color, # observed
                Rowv = as.dendrogram(hc_result),
                Colv = as.dendrogram(hc_result),
                trace = "none",
                key = FALSE,
                labRow = NA,
                labCol = NA)
      dev.off()

      dir <- paste(save_dir, "aurora_heatmap_median_random_forest_",
                   Sys.Date(), ".png", sep = "")

      dist_matrix <- as.dist(final_proxy_mat)
      hc_result <- stats::hclust(dist_matrix, method = "ward.D")

      png(filename=dir)
      gplots::heatmap.2(final_proxy_mat,
                ColSideColors = col_color, # predicted
                RowSideColors = row_color, # observed
                Rowv = as.dendrogram(hc_result),
                Colv = as.dendrogram(hc_result),
                trace = "none",
                key = FALSE,
                labRow = NA,
                labCol = NA,
                na.rm = TRUE)
      dev.off()



      colors = c("red",
                 "blue",
                 "green",
                 "black",
                 "darkgoldenrod",
                 "lightpink",
                 "darkorange",
                 "cyan",
                 "deeppink",
                 "antiquewhite",
                 "aquamarine")
      heatmap_colors <- lookup_tbl
      # heatmap() does not have a legend for the colors so return this object
      names(heatmap_colors) <- colors[1:length(lookup_tbl)]
    }
  }

  if (CART == TRUE) {
    # analyse the proximities from CART
    no_samples <- nrow(list_proxy_mats_CART[[1]])
    final_proxy_mat_CART <- matrix(NA, nrow = no_samples, ncol = no_samples )
    colnames(final_proxy_mat_CART) <- colnames(list_proxy_mats_CART[[1]])
    rownames(final_proxy_mat_CART) <- colnames(list_proxy_mats_CART[[1]])

    for (i in 1:no_samples) {
      for (j in 1:no_samples) {
        vect_dist <- c()
        for (k in 1:length(list_proxy_mats_CART)) {
          vect_dist <- c(vect_dist, list_proxy_mats_CART[[k]][j,i])
        }
        vect_dist <- as.numeric(vect_dist)
        vect_dist <- vect_dist[!is.na(vect_dist)]
        final_proxy_mat_CART[j,i] <- stats::median(as.numeric(vect_dist))
      }
    }

    # plot the results
    if (CART_plot == TRUE & length(phenotypes) > 11 & write_data == TRUE) {
      CART_plot <- FALSE
      warning("The number of analysed phenotypes is higher then 11. Results of CART will not be plotted. \n")
    }

    if (CART_plot == TRUE & write_data == TRUE) {
      # the maximum pheno is 11 because then we run out of colors
      auto_colors = c("darkred",
                      "blue4",
                      "green4",
                      "black",
                      "darkgoldenrod4",
                      "lightpink3",
                      "darkorange4",
                      "cyan3",
                      "deeppink4",
                      "antiquewhite4",
                      "aquamarine4")

      allo_colors = c("firebrick1",
                      "blue1",
                      "green1",
                      "gray20",
                      "darkgoldenrod1",
                      "lightpink1",
                      "darkorange1",
                      "cyan1",
                      "deeppink1",
                      "antiquewhite1",
                      "aquamarine1")
      col_color <- c()
      row_color <- c()
      for (i in 1:nrow(df_final_CART)) {
        index <- which(lookup_tbl == df_final_CART$observed[i])
        row_color <- c(row_color, auto_colors[index])
        if (df_final_CART$predicted[i] == df_final_CART$observed[i]) {
          if (df_final_CART[i, 1+index] == "TRUE not typical") {
            col_color <- c(col_color, allo_colors[index])
          } else {
            col_color <- c(col_color, auto_colors[index])
          }
        } else {
          col_color <- c(col_color, "white")
        }
      }

      dir <- paste(save_dir, "aurora_heatmap_CART_",
                   Sys.Date(), ".png", sep = "")

      if (any(is.na(final_proxy_mat_CART)|is.nan(final_proxy_mat_CART)|is.infinite(final_proxy_mat_CART))) {
        # if there is an NA in the matrix then impute the missing values
        final_proxy_mat_CART[is.infinite(final_proxy_mat_CART)] <- NA
        final_proxy_mat_CART[is.nan(final_proxy_mat_CART)] <- NA
        for (i in 1:ncol(final_proxy_mat_CART)) {
          med <- median(final_proxy_mat_CART[, i][!is.na(final_proxy_mat_CART[, i])])
          if (is.na(med)) {med <- 0}
          final_proxy_mat_CART[is.na(final_proxy_mat_CART[, i]), i] <- med
        }
      }

      dist_matrix <- as.dist(final_proxy_mat_CART)
      hc_result <- stats::hclust(dist_matrix, method = "ward.D", )

      png(filename=dir)
      gplots::heatmap.2(final_proxy_mat_CART,
                ColSideColors = col_color, # predicted
                RowSideColors = row_color, # observed
                Rowv = as.dendrogram(hc_result),
                Colv = as.dendrogram(hc_result),
                trace = "none",
                key = FALSE,
                labRow = NA,
                labCol = NA,
                na.rm = TRUE)
      dev.off()

      colors = c("red",
                 "blue",
                 "green",
                 "black",
                 "darkgoldenrod",
                 "lightpink",
                 "darkorange",
                 "cyan",
                 "deeppink",
                 "antiquewhite",
                 "aquamarine")

      heatmap_colors <- lookup_tbl
      # heatmap() does not have a legend for the colors so return this object
      names(heatmap_colors) <- colors[1:length(lookup_tbl)]
    }
  }
  # create empty list that will store all the results
  AURORA_results <- vector("list")
  # return list with all important objects
  AURORA_results <- rlist::list.append(AURORA_results, call = match.call()) # AURORA call with all args
  AURORA_results <- rlist::list.append(AURORA_results, bin_mat = bin_mat_original) # original bin_mat before any filtering
  AURORA_results <- rlist::list.append(AURORA_results, phenotypes = unname(lookup_tbl)) # all analysed phenotypes
  # put results from all tools here
  results_nested <- vector("list")
  if (random_forest == TRUE) {
    results_random_forest <- vector("list")
    results_random_forest <- rlist::list.append(results_random_forest, p_val_mat = ks_mat_rf)
    results_random_forest <- rlist::list.append(results_random_forest, auto_allo_results = result_probs_rf)
    results_random_forest <- rlist::list.append(results_random_forest, plot = plot_rf)
    results_random_forest <- rlist::list.append(results_random_forest, distance_heatmap_median = final_proxy_mat)
    results_random_forest <- rlist::list.append(results_random_forest, distance_heatmap_sum = final_proxy_mat_sum)
    # add plots for these heatmaps
    if (plot_random_forest == TRUE & write_data == TRUE) {
      results_random_forest <- rlist::list.append(results_random_forest, legend_heatmap = heatmap_colors)
    }
    results_random_forest <- rlist::list.append(results_random_forest, gfs_importance = df_final_gfs_rf)
    results_random_forest <- rlist::list.append(results_random_forest, aucs = aucs_rf)

    results_nested <- rlist::list.append(results_nested, results_random_forest = results_random_forest)
  }
  if (adaboost == TRUE) {
    results_adaboost <- vector("list")
    results_adaboost <- rlist::list.append(results_adaboost, p_val_mat = ks_mat_ada)
    results_adaboost <- rlist::list.append(results_adaboost, auto_allo_results = result_probs_ada)
    results_adaboost <- rlist::list.append(results_adaboost, plot = plot_ada)
    results_adaboost <- rlist::list.append(results_adaboost, gfs_importance = df_final_gfs_ada)
    results_adaboost <- rlist::list.append(results_adaboost, aucs = aucs_ada)

    results_nested <- rlist::list.append(results_nested, results_adaboost = results_adaboost)

  }

  if (ovr_log_reg == TRUE) {
    results_log_reg <- vector("list")
    results_log_reg <- rlist::list.append(results_log_reg, p_val_mat = ks_mat_log_reg)
    results_log_reg <- rlist::list.append(results_log_reg, auto_allo_results = result_probs_log_reg)
    results_log_reg <- rlist::list.append(results_log_reg, plot = plot_log_reg)
    results_log_reg <- rlist::list.append(results_log_reg, gfs_importance = df_final_log_reg)
    results_log_reg <- rlist::list.append(results_log_reg, aucs = aucs_log_reg)

    results_nested <- rlist::list.append(results_nested, results_log_reg = results_log_reg)
  }

  if (CART == TRUE) {
    results_CART <- vector("list")
    results_CART <- rlist::list.append(results_CART, p_val_mat = ks_mat_CART)
    results_CART <- rlist::list.append(results_CART, auto_allo_results = df_final_CART)
    results_CART <- rlist::list.append(results_CART, plot = plot_CART)
    results_CART <- rlist::list.append(results_CART, complexity = df_one_vs_one)
    results_CART <- rlist::list.append(results_CART, gfs_importance = df_importances_CART)
    results_CART <- rlist::list.append(results_CART, distance_heatmap = final_proxy_mat_CART)

    if (CART_plot == TRUE & write_data == TRUE) {
      results_CART <- rlist::list.append(results_CART, legend_heatmap = heatmap_colors)
    }
    results_nested <- rlist::list.append(results_nested, results_CART = results_CART)
  }

  AURORA_results <- rlist::list.append(AURORA_results, results = results_nested)
  # no of gfs removed in the initial filter (low_perc_cutoff, upp_perc_cutoff)
  AURORA_results <- rlist::list.append(AURORA_results, abundance_filter_removed = to_bin_result$no_removed)
  if (run_chisq == TRUE) {
    # no of gfs removed in chi-square test filter
    AURORA_results <- rlist::list.append(AURORA_results, chisq_filter_removed = to_bin_result$chisq_removed)
  }
  if (ancest_rec_filter == TRUE) {
    AURORA_results <- rlist::list.append(AURORA_results, acr_filter_removed = post_proc_data$ancent_filt_col)
    AURORA_results <- rlist::list.append(AURORA_results, plot_acr_cutoff = post_proc_data$p_acr)
  }

  if (jaccard_filter == TRUE) {
    # no of gfs removed is sum of all collapsed gfs - no of new gfs
    AURORA_results <- rlist::list.append(AURORA_results, jaccard_filter_removed = post_proc_data$jacc_filt_rem)
  }

  if (hamming_filter == TRUE) {
    # no of gfs removed is sum of all collapsed gfs - no of new gfs
    AURORA_results <- rlist::list.append(AURORA_results, hamming_filter_removed = post_proc_data$hamm_filt_rem)
  }

  if (reduce_outlier == TRUE) {
    # output a histogram with pairwise distances and the cutoff value (blue line)
    AURORA_results <- rlist::list.append(AURORA_results, pairwise_dist_outlier = result$p_outlier)
  }

  if (fit_parameters == TRUE) {
    # return the fitted parameters
    fitted_parameters <- vector("list")
    if (random_forest == TRUE) {
      parameters_random_forest <- as.list(setNames(c(sampsize, mtry, ntree, maxnodes),
                                                   c("sampsize", "mtry", "ntree", "maxnodes")))
      fitted_parameters <- rlist::list.append(fitted_parameters,
                                              parameters_random_forest = parameters_random_forest)
    }

    if (adaboost == TRUE) {
      parameters_adaboost <- as.list(setNames(c(max_depth, n_estimators, learning_rate),
                                              c("max_depth", "n_estimators", "learning_rate")))

      fitted_parameters <- rlist::list.append(fitted_parameters,
                                              parameters_adaboost = parameters_adaboost)
    }
    AURORA_results <- rlist::list.append(AURORA_results, fitted_parameters = fitted_parameters)
  }

  AURORA_results <- rlist::list.append(AURORA_results, bonferroni_threshold = bonferroni_threshold)

  cat("Finished \n")
  return(AURORA_results)
}
