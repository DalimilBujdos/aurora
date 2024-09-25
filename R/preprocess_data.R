#' Preprocess the binary matrix
#'
#' @noRd
preprocess_data <- function(bin_mat,
                            pheno_mat,
                            tree = NA,
                            ancest_rec_filter = TRUE,
                            cutoff_asr = 2,
                            jaccard_filter = FALSE,
                            eps_val = 0.05,
                            minPts_val = 3,
                            hamming_filter = TRUE,
                            hamming_cutoff = 3) {


  if (is.list(tree)) {
    if (!is.null(tree$node.label)) {
      # print warning messages if the tree quality is bad
      bootstrap <- length(which(as.numeric(tree$node.label) < 75))
      if (bootstrap > 1) {
        warning("The tree has ", bootstrap ," nodes with bootstrap values below 75\n")
      }
    }

    if (!is.null(tree$edge.length)) {
      # test if there are not branches that are unusually long
      total_len <- sum(tree$edge.length)
      tree_len <- length(which(any(tree$edge.length > total_len*0.05)))
      if (tree_len > 1) {
        warning("The tree has ", tree_len ," edges that are longer than 5% of the total tree edge length. These edges might point to outliers.\n")
      }
    }
  }

  messa <- ncol(bin_mat)

  if (ancest_rec_filter == TRUE) {
    # right now node labels in tree are just bootstrap values
    tree$node.label <- paste0("node", seq(1:tree$Nnode))

    # initiate a df that will hold distribution of the final values
    distrib_gfs <- data.frame(gfs_id = colnames(bin_mat),
                             orignis_no = rep(NA, ncol(bin_mat)),
                             area = rep(NA, ncol(bin_mat)))

    for (i in 1:ncol(bin_mat)) {
      df_bin_mat <- data.frame(ids = rownames(bin_mat),
                               pheno = bin_mat[,i])

      pheno <- data.frame(ids = tree$tip.label,
                          pheno = df_bin_mat$pheno[match(tree$tip.label, df_bin_mat$ids)])


      pheno$pheno <- as.factor(ifelse(pheno$pheno == 1, "B", "A"))
      pheno_input <- pheno$pheno
      names(pheno_input) <- pheno$ids
      # use treeWAS wrapper around ape function ace to calculate ancestral reconstruction
      output_was <- suppressWarnings(treeWAS::asr(pheno_input, tree))

      select.tip.or.node <- function(element, tree) {
        ifelse(element < ape::Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-ape::Ntip(tree)])
      } # everything that has num in tree$edge below Ntip(tree)+1 is a tip

      edge_table <- data.frame(
        parent = tree$edge[,1],
        parent_name = sapply(tree$edge[,1], select.tip.or.node, tree = tree),
        parent_positive = rep(FALSE, length(tree$edge[,1])),
        child = tree$edge[,2],
        child_name = sapply(tree$edge[,2], select.tip.or.node, tree = tree),
        child_positive = rep(FALSE, length(tree$edge[,1])),
        edge_length = tree$edge.length
      )
      edge_table$parent_positive[is.element(edge_table$parent_name,
                                            names(output_was)[output_was == 1])] <- TRUE
      edge_table$child_positive[is.element(edge_table$child_name,
                                           names(output_was)[output_was == 1])] <- TRUE

      # find how many times a negative parent node had a child that is positive
      edge_table_origins <- edge_table[edge_table$parent_positive == FALSE &
                                         edge_table$child_positive == TRUE, ]
      distrib_gfs[i,2] <- length(unique(edge_table_origins$parent))

      # now calculate the area of the tree covered by the ancestral state reconstruction
      # of the gf
      edge_table_area <- edge_table[edge_table$parent_positive == TRUE &
                                       edge_table$child_positive == TRUE, ]

      distrib_gfs[i,3] <- sum(edge_table_area$edge_length)

    }
    # calculate ratio: number of times a trait originated in evolution/total area
    # the acr covers. This ratio will be very large for genes that are transferred by
    # horizontal gene transfer
    # the ratio is an exponential value so we will log it

    # if orignis_no it mans that the root of the tree is origin. Covert all 0 to 1
    distrib_gfs$orignis_no[distrib_gfs$orignis_no == 0] <- 1
    distrib_gfs$ratio <- log(distrib_gfs$orignis_no/distrib_gfs$area)
    # some elements are infinite because the area is 0
    #distrib_gfs$ratio[!is.finite(distrib_gfs$ratio)] <- max(distrib_gfs$ratio[is.finite(distrib_gfs$ratio)])
    distrib_gfs$zscore <- (distrib_gfs$ratio - mean(distrib_gfs$ratio[is.finite(distrib_gfs$ratio)]))/
      stats::sd(distrib_gfs$ratio[is.finite(distrib_gfs$ratio)])


    remove_gfs <- distrib_gfs$gfs_id[distrib_gfs$zscore > cutoff_asr]

    ancent_filt_col <- length(remove_gfs)
    cat("Ancestral reconstruction filter eliminated", ancent_filt_col
        ,"variants out of", messa, "\n")
    remove_col <- which(is.element(colnames(bin_mat), remove_gfs))
    bin_mat <- bin_mat[, -remove_col]

    # generate the histogram with cutoff
    distrib_gfs <- distrib_gfs[is.finite(distrib_gfs$zscore), ]
    # Basic histogram
    p_acr <- ggplot2::ggplot(distrib_gfs, aes(x=zscore)) +
                             geom_histogram(binwidth=0.3, color="black", fill="white") +
                             geom_vline(aes(xintercept=cutoff_asr),
                             color="blue", linetype="dashed", linewidth=1) +
                             theme_classic() +
                             ggtitle("Remove variants that are spread around the phylogenetic tree")

  }
  messa <- ncol(bin_mat)


  if (jaccard_filter == TRUE) {
    # calculate a distance matrix for gene families
    dist <- proxy::dist(t(bin_mat),
                        method = "Jaccard")

    cat("Jaccard calculation complete \n")

    # convert object "dist"  to matrix for later calculation

    dist_mat <- as.matrix(dist)

    # run dbscan
    # the purpose of dbscan is to find gfs that are inherrited together
    # these gfs will be clustered into one representative point. Thus, we
    # deal with colinearity probelm in rf and adaboost
    dbscan_res <- dbscan::dbscan(dist,
                                 eps = eps_val,
                                 minPts = minPts_val)


    cat("DBSCAN complete \n")

    # now choose the most representative sequence from the idenfidied cluster

    clusters <- unique(dbscan_res$cluster)
    clusters <- clusters[!(clusters == 0)]# here zero indicates noise

    cluster_vect <- rep(0, dim(dist_mat)[1])

    no_new_gfs <- 0

    # calculate the representative points
    for (i in 1:length(clusters)) {
      indexes <- which(dbscan_res$cluster == clusters[i])

      distances <- data.frame(gf_names = rownames(dist_mat)[indexes],
                              distances = rep(NA, length(indexes)))

      # calculate which gf has the smallest distance from other gfs
      for (j in 1:nrow(distances)) {
        set <- dist_mat[rownames(dist_mat) == distances$gf_names[j],
                             is.element(colnames(dist_mat), distances$gf_names)]
        distances$distances[j] <- sum(set)
      }


      # remove the other_gf
      index <- which(is.element(colnames(bin_mat),
                                distances$gf_names[-which.max(distances$distances)]))

      bin_mat <- bin_mat[,-index]

      # the name of the new aggregate gene family is separated by ___
      new_gf_name1 <- distances$gf_names[-which.max(distances$distances)]
      new_gf_name2 <- distances$gf_names[which.max(distances$distances)]
      new_gf_name  <- paste(c(new_gf_name2, new_gf_name1), sep = "___", collapse = "___")

      no_new_gfs  <- no_new_gfs + 1

      # change the name of the reprez_seq
      index <- which(colnames(bin_mat) == distances$gf_names[which.max(distances$distances)])
      colnames(bin_mat)[index] <- new_gf_name

      # create also a vector that will assign a cluster number to gene families in
      # dist_mat_gif - for figure

      indexes <- which(is.element(colnames(dist_mat), distances$gf_names))

      cluster_vect[indexes] <- i

    }

    cat("DBSCAN with parameters minPts:", minPts_val, "eps:", eps_val , "collapsed",
        messa - ncol(bin_mat), "variants. The total number of variants left is", ncol(bin_mat), "\n")

    jacc_filt_rem <- (messa - ncol(bin_mat)) - no_new_gfs

    # now run the mds plot
    mds <- stats::cmdscale(dist_mat, eig = TRUE)
    mds_df <- as.data.frame(mds$points)
    mds_df$keys <- rownames(mds_df)
    colnames(mds_df)[c(1,2)] <- c("dim_1", "dim_2")

    colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
                "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ff33cc", "#99ccff",
                "#66ff66", "#ff6666", "#ffb3e6")

    color_df <- setNames(colors[(1:max(dbscan_res$cluster) - 1) %% length(colors) + 1], 1:max(dbscan_res$cluster))

    p_jacc <- ggplot2::ggplot() +
                       geom_point(data = mds_df, aes(x = dim_1, y = dim_2, color = as.character(cluster_vect))) +
                       scale_colour_manual(values=color_df) +
                       labs(y= "MDS2", x = "MDS1") +
                       theme_classic() +
                       theme(legend.position = "none") +
                       theme(axis.text = element_text(size = 12))

  }

  if (hamming_filter == TRUE) {
    # hamming filter is a filter that clusters only gene families that have
    # less or equal number of missmatching positions then "hamming_cutoff"

    no_collapsed_gfs <- c()
    no_new_gfs <- 0

    hamm_count <- function(num, bin_mat) {
      vect <- bin_mat[, num]
      k <- apply(bin_mat, 2, function(x) length(which(vect != x)))
      return(k)
    }
    # elements of mat are the number of 1/0 pairs between gene families
    hamming_mat <- sapply(1:ncol(bin_mat), hamm_count, bin_mat)
    colnames(hamming_mat) <- colnames(bin_mat)

    list_clusters <- vector(mode = "list")

    bigger_then_cutoff <- function(vect) {
      bigger <- length(which(vect>hamming_cutoff))
      return(bigger)
    }

    # calculates the clusters
    decide <- TRUE
    while(decide) {
      if(!is.matrix(hamming_mat)) {
        break
      } else if (ncol(hamming_mat) <= 1) {
        break
      }

      index <- which(hamming_mat[,1] <= hamming_cutoff)
      if (length(index) == 1) {
        hamming_mat <- hamming_mat[-index,  -index]

        if(!is.matrix(hamming_mat)) {
          break
        } else if (ncol(hamming_mat) <= 1) {
          break
        } else {
          next
        }
      }
      hamming_mat_cluster <- hamming_mat[index,  index]
      bigger <- 1
      # chooses only cluster where all gfs are have only "<= hamming_cutoff" 1/0 pairs
      while(sum(bigger) > 0) {
        bigger <- apply(hamming_mat_cluster,2 , bigger_then_cutoff)
        if (sum(bigger) > 0) {
          hamming_mat_cluster <- hamming_mat_cluster[-which.max(bigger), -which.max(bigger)]
        }
      }

      # choose the gf that has the lowest number of 1/0 pairs and this will be the
      # representative gf
      index_name <- which.min(rowSums(hamming_mat_cluster))
      new_gf_name <-colnames(hamming_mat_cluster)[index_name]
      new_gf_name_2 <- colnames(hamming_mat_cluster)[-index_name]
      new_gf_name <- paste(c(new_gf_name, new_gf_name_2), sep = "___", collapse = "___")
      # find the indexes of the new cluster in bin_mat
      # the index of the representative gfs should be the last
      index_gfs <- which(is.element(colnames(bin_mat), colnames(hamming_mat_cluster)))
      index_gfs <- index_gfs[-index_name]
      index_gfs <- c(index_gfs, which(is.element(colnames(bin_mat),  names(index_name))))

      list_clusters[[length(list_clusters)+1]] <- index_gfs
      names(list_clusters)[length(list_clusters)] <- new_gf_name

      # remove the clustered gfs from the hamming_mat
      index <- which(is.element(colnames(hamming_mat), colnames(hamming_mat_cluster)))
      hamming_mat <- hamming_mat[-index,  -index]

      no_collapsed_gfs <- c(no_collapsed_gfs, nrow(hamming_mat_cluster))
      no_new_gfs <- no_new_gfs+1
    }
    # remove the clustered gfs and add the new representative gf

    cols_to_remove <- c()
    cols_to_add <- data.frame(row.names = row.names(bin_mat))
    for (i in 1:length(list_clusters)) {
      if (length(list_clusters) == 0) {break}
      repres_gf <- bin_mat[, list_clusters[[i]][length(list_clusters[[i]])] ]
      cols_to_remove <- c(cols_to_remove, list_clusters[[i]])
      cols_to_add <- cbind(cols_to_add, as.data.frame(repres_gf))
      names(cols_to_add)[length(cols_to_add)] <- names(list_clusters[i])
    }
    if (length(list_clusters) > 0) {
      bin_mat <- bin_mat[,-cols_to_remove]
      bin_mat <- cbind(bin_mat, cols_to_add)
    }

    cat("Hamming filter collapsed", sum(no_collapsed_gfs), "variants into",
        no_new_gfs, "new variants. The number of variants left is:",
        ncol(bin_mat), "\n")
    hamm_filt_rem <- sum(no_collapsed_gfs) - no_new_gfs
  }

  return_obj <- list(bin_mat = bin_mat,
                     ancent_filt_col =  if (exists("ancent_filt_col")) ancent_filt_col else FALSE,
                     p_acr = if (exists("p_acr")) p_acr else FALSE,
                     jacc_filt_rem = if (exists("jacc_filt_rem")) jacc_filt_rem else FALSE,
                     p_jacc = if (exists("p_jacc")) p_jacc else FALSE,
                     hamm_filt_rem = if (exists("hamm_filt_rem")) hamm_filt_rem else FALSE)

  return(return_obj)
}
