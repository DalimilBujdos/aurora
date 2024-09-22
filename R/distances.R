#' Process the phylogenetic object and reduce distances to outliers
#'
#' @noRd
distances <- function(tree = NA,
                      alternative_dist_mat = NA,
                      bagging,
                      reduce_outlier = TRUE,
                      cutoff_outlier = 3) {

  if (is.matrix(alternative_dist_mat)) {
    # alternative mat
    dist_mat <- apply(alternative_dist_mat, 2, as.numeric)
    # rownames are lost after apply functional is run
    rownames(dist_mat) <- colnames(dist_mat)
    if (bagging == "phylogenetic_walk") {
      tree <- ape::nj(dist_mat)
    }
  }  else if (bagging == "random_walk") {
    dist_mat <- ape::cophenetic.phylo(tree)
  }

  if (bagging == "random_walk") {
    # if bagging algo is random walk then we reduce outliers in dist_mat

    dist_mat_plot <- dist_mat

    # convert the df to triangle
    for (i in 1:ncol(dist_mat_plot)) {
      for (j in 1:nrow(dist_mat_plot)) {
        if (j>=i) {
          dist_mat_plot[i,j] <- NA
        }
      }
    }

    # reduce the outliers by calculating which pairwise distances are large
    # these are the reduced
    if (reduce_outlier == TRUE) {
      # turn the dist_mat_plot into a dataframe

      df_dist_mat <- data.frame(column_id = as.numeric(),
                                row_id = as.numeric(),
                                distance = as.numeric())
      # this loop turns distance matrix into dataframe so we can calculate the z-score
      for (i in 1:ncol(dist_mat_plot)) {
        row_indexes <- which(!is.na(dist_mat_plot[,i]))
        # get values for the column
        df_temp <- data.frame(column_id = rep(i, length(row_indexes)),
                              row_id = row_indexes,
                              distance = dist_mat_plot[row_indexes,i])

        df_dist_mat <- rbind(df_dist_mat, df_temp)
      }

      # compute z-score for the distance values
      df_dist_mat <- df_dist_mat %>%
        dplyr::mutate(zscore = ((df_dist_mat[,3] - mean(df_dist_mat[,3]))/stats::sd(df_dist_mat[,3])))

      # this hist tells me that the distances are not distrib equally and normaliz
      # will be applied

      p_outlier <- ggplot2::ggplot(df_dist_mat, aes(x=zscore)) +
                                   geom_histogram(binwidth=0.3, color="black", fill="white") +
                                   geom_vline(aes(xintercept=cutoff_outlier),
                                    color="blue", linetype="dashed", linewidth=1) +
                                   ggtitle("Reduce pairwise distancetes to the cutoff value")


      # now find points that are above the cutoff_outlier*sd (99.7% for the def value)
      # and lower their value
      outliers <- which(df_dist_mat$zscore > cutoff_outlier)
      # There is a chance that there are not gonna be any outliers - in that case
      # just return the dist_mat

      # calculate the new values in case there are outliers. The newvalue is simply
      # mean+cutoff_outlier*SD
      new_val <- mean(df_dist_mat[,3])+(cutoff_outlier*stats::sd(df_dist_mat[,3]))
      if (length(outliers) > 0) {
        strains <- rownames(df_dist_mat)[which(df_dist_mat$zscore > cutoff_outlier)]
        cat(length(strains), "distances were identified as outliers. These distance were reduced \n")
        df_dist_mat_out <- df_dist_mat[outliers,]

        # remove the original values of a point and add a new_val
        for (i in 1:nrow(df_dist_mat_out)) {
          dist_mat[df_dist_mat_out[i,1], df_dist_mat_out[i,2]] <- new_val
          dist_mat[df_dist_mat_out[i,2], df_dist_mat_out[i,1]] <- new_val
        }

        # some values in the dist_mat are 0. This creates a problems in the random walk
        # algorithm. Instead a 0 give them the smallest non-zero value from dist_mat
        # this also means that it is possible to walk from one point to the same point

        dist_mat[dist_mat == 0] <- min(df_dist_mat$distance[df_dist_mat$distance > 0])

      }
    }
    # Exclude 0. Same as above
    if (reduce_outlier == FALSE) {
      num_dist_mat <- as.numeric(dist_mat)
      dist_mat[dist_mat == 0] <- min(num_dist_mat[num_dist_mat > 0])
    }

    return(list(dist_mat_or_tree = dist_mat,
                p_outlier = if (exists("p_outlier")) p_outlier else FALSE))
  } else {
    # if bagging algo is phylogenetic_walk then we reduce the length of the phylo tree

    if (reduce_outlier == TRUE) {
      # the edges have to be loged cuz most trees have many short branch lengths and
      # just a few long
      # sometimes some edge lenghts are negative which is a problem when we want to apply log
      tree$edge.length[tree$edge.length <=0 ] <- 0
      log_branch <- data.frame(log_length = log(tree$edge.length, base = 10),
                               selected = rep(TRUE, length(tree$edge.length)))
      # is possible that some vals are -inf. Replace it with the lowest value in the tree
      lowest_val <- log_branch$log_length[log_branch$log_length != -Inf]
      lowest_val <- lowest_val[which.min(lowest_val)]
      log_branch$log_length[log_branch$log_length == -Inf] <- lowest_val
      # next we need to remove extremly short branch lengths. We can do this by
      # calculating the percentage of tree area a branch represents and then remove
      # branches that represent less than 100/length(tree$edge.length) of the total tree area
      log_branch_perc <- (log_branch$log_length/sum(log_branch$log_length)) * 100
      log_branch$selected[log_branch_perc < 100/length(tree$edge.length)] <- FALSE
      # calculate IQR. outliers are branches that are longer or equal to median + cutoff_outlier*IQR
      iqr <- stats::IQR(log_branch$log_length[log_branch$selected == TRUE])
      upper_cutoff <- stats::median(log_branch$log_length[log_branch$selected == TRUE]) + cutoff_outlier * iqr
      outliers <- which(log_branch$log_length >= upper_cutoff)
      # transform the cutoff into a real branch length
      branch_upper_cutoff <- 10^upper_cutoff
      df_plot <- data.frame(branch_lengths = tree$edge.length)
      p_outlier <- ggplot2::ggplot(df_plot, aes(x=branch_lengths)) +
                            geom_histogram(bins=30, color="black", fill="white") +
                            geom_vline(aes(xintercept=branch_upper_cutoff),
                              color="blue", linetype="dashed", linewidth=1) +
                            ggtitle("Reduce tree branch lengths to the cutoff value")


      # calculate the new values in case there are outliers. The newvalue is simply
      # the branch_upper_cutoff
      if (length(outliers) > 0) {
        tree$edge.length[outliers] <- branch_upper_cutoff
        cat(length(outliers), "branches were identified as outliers. These branches were reduced \n")
      }
    }
    return(list(dist_mat_or_tree = tree,
                p_outlier = if (exists("p_outlier")) p_outlier else FALSE))
  }
}
