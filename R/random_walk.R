#' Random walk algorithm that obtains training datasets
#'
#' @noRd
random_walk <- function(bags,
                        dist_mat_class,
                        misslabeled_str,
                        bag_size,
                        max_per_bag,
                        max_misslablel = FALSE,
                        no_rounds) {
  # this will measure how many times a certain strain is repeated in all bags
  count <- setNames(rep(0, nrow(dist_mat_class)), row.names(dist_mat_class))

  # the loop that makes the bags
  for (i in 1:no_rounds) {
    # dist_mat is reneved for every run
    dist_mat_class_tmp <- dist_mat_class
    # initialize the first sample
    sample_id <- base::sample(1:nrow(dist_mat_class_tmp), size = 1L)
    # initialize new count for this bag
    count_bag <- setNames(rep(0, nrow(dist_mat_class_tmp)),
                          row.names(dist_mat_class_tmp))
    # this is here to ensure that the each genomes appears only "max_per_bag"
    # amount of times in the bag. it also ensures that misslabelled strains
    # does not overtake the non-misslabelled

    # add that to the bag
    bags[[i]][length(bags[[i]]) + 1] <- rownames(dist_mat_class_tmp)[sample_id]
    # count the sample into the two counters
    count[sample_id] <- count[sample_id] + 1
    count_bag[sample_id] <- 1

    vect_prob <- dist_mat_class_tmp[sample_id, ]

    for (j in 2:bag_size) {
      # now find a new sample
      sample_id <- base::sample(1:length(vect_prob), size = 1L, prob = vect_prob)
      bags[[i]][length(bags[[i]]) + 1] <- names(vect_prob)[sample_id]

      vect_prob <- dist_mat_class_tmp[sample_id, ]

      count_bag_index <- which(names(count_bag) == bags[[i]][length(bags[[i]])])
      count_bag[count_bag_index] <- count_bag[count_bag_index] + 1

      # remove those samples that are in bag more than max_per_bag times
      if (count_bag[count_bag_index] == max_per_bag) {
        name <- names(count_bag)[count_bag_index]
        # update vect_prob
        index <- which(names(vect_prob) == name)
        vect_prob <- vect_prob[-index]
        # update dist_mat_class_tmp
        index <- which(colnames(dist_mat_class_tmp) == name)
        dist_mat_class_tmp <- dist_mat_class_tmp[-index, -index]
      }

      # remove misslabelled strains that are in bag more than max_misslablel times
      if (is.numeric(max_misslablel)) {
        name <- names(count_bag)[count_bag_index]
        if (count_bag[count_bag_index] == max_misslablel &
            any(misslabeled_str == name)) {
          index <- which(names(vect_prob) == name)
          vect_prob <- vect_prob[-index]
          index <- which(colnames(dist_mat_class_tmp) == name)
          dist_mat_class_tmp <- dist_mat_class_tmp[-index, -index]
        }
      }
      # sample_id and sample_id2 are different because of the two criteria count_bag and count
      sample_id2 <- which(names(count) == bags[[i]][length(bags[[i]])])
      count[sample_id2] <- count[sample_id2] + 1
    }
    # its possible that random walk did not select the misslabelled strians. If thats
    # the case add them here
    if (is.character(misslabeled_str)) {
      index <- (length(bags[[i]]) - bag_size + 1):length(bags[[i]])
      for (j in 1:length(misslabeled_str)) {
        if (count_bag[names(count_bag) == misslabeled_str[j]] == 0) {
          replace_str <- base::sample(index, 1)
          index <- index[-replace_str]
          count_bag[names(count_bag) == misslabeled_str[j]] <- 1
          rmv_strain <- bags[[i]][replace_str]
          count_bag[names(count_bag) == rmv_strain] <- count_bag[names(count_bag) == rmv_strain] -1
          bags[[i]][replace_str] <- misslabeled_str[j]
        }
      }
    }
  }
  return(list(bags = bags, count = count))

}
