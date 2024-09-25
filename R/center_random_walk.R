center_random_walk <- function(bags,
                               dist_mat_class,
                               misslabeled_str,
                               bag_size,
                               max_per_bag,
                               max_misslablel = FALSE,
                               no_rf) {
  # this will measure how many times a certain strain is repeated in all bags
  count <- setNames(rep(0, nrow(dist_mat_class)), row.names(dist_mat_class))

  # Function to calculate stress for a given number of dimensions
  calc_stress <- function(k) {
    mds <- suppressWarnings(cmdscale(dist_mat_class, k = k))
    stress <- sqrt(sum((as.dist(dist_mat_class) - dist(mds))^2)) / sqrt(sum(as.dist(dist_mat_class)^2))
    return(stress)
  }

  # Calculate stress for different numbers of dimensions
  dims <- c(nrow(dist_mat_class) - 1, 50)
  max_dim <- dims[which.min(dims)]
  stress_values <- sapply(1:max_dim, calc_stress)

  # take the first no of dim where the stress value is below 0.005 or the lowest
  # stress value
  if (any(stress_values <= 0.01)) {
    dim_mds <- which(stress_values <= 0.01)[1]
  } else {
    dim_mds <- which.min(stress_values)
  }

  # Apply MDS to the distance matrix with the new no of dimensions
  mds <- suppressWarnings(cmdscale(dist_mat_class, k = dim_mds))

  # Calculate the centroid of the MDS coordinates
  centroid <- apply(mds, 2, mean)
  distances <- rbind(t(data.frame(centroid = centroid)),
                     as.data.frame(mds))
  # this calculates euclidean distance between the centroid and the point
  distances1 <- as.matrix(dist(as.matrix(distances)))

  # Extract distances to all strains from the centroid
  distances_to_points <- as.numeric(distances1[1, 2:ncol(distances1)])
  names(distances_to_points) <- colnames(distances1)[-1]
  distances_to_points_save <- distances_to_points

  # the loop that makes the bags
  for (i in 1:no_rf) {
    # initialize new count for this bag
    count_bag <- setNames(rep(0, length(distances_to_points)),
                          names(distances_to_points))
    # this is here to ensure that the each genomes appears only "max_per_bag"
    # amount of times in the bag. it also ensures that misslabelled strains
    # does not overtake the non-misslabelled
    for (j in 1:bag_size) {
      distances_to_points <- distances_to_points_save
      # draw the sample
      sample_id <- base::sample(1:length(distances_to_points), size = 1L, prob = distances_to_points)
      # add that to the bag
      bags[[i]][length(bags[[i]]) + 1] <- names(distances_to_points)[sample_id]
      # count the sample into the two counters
      index_count <- which(names(count) == bags[[i]][length(bags[[i]])])
      count[index_count] <- count[index_count] + 1
      count_bag_index <- which(names(count_bag) == bags[[i]][length(bags[[i]])])
      count_bag[count_bag_index] <- count_bag[count_bag_index] + 1


      # remove those samples that are in bag more than max_per_bag times
      if (count_bag[count_bag_index] == max_per_bag) {
        name <- names(count_bag)[count_bag_index]
        # update distances_to_points
        index <- which(names(distances_to_points) == name)
        distances_to_points <- distances_to_points[-index]
      }

      # remove misslabelled strains that are in bag more than max_misslablel times
      if (is.numeric(max_misslablel)) {
        name <- names(count_bag)[count_bag_index]
        if (count_bag[count_bag_index] == max_misslablel &
            any(misslabeled_str == name)) {
          index <- which(names(distances_to_points) == name)
          distances_to_points <- distances_to_points[-index]
        }
      }
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
          count_bag[names(count_bag) == rmv_strain] <- count_bag[names(count_bag) == rmv_strain] - 1
          bags[[i]][replace_str] <- misslabeled_str[j]
        }
      }
    }
  }
  return(list(bags = bags, count = count))
}
