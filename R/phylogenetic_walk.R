#' Phylogenetic walk algorithm that obtains training datasets
#'
#' @noRd
phylogenetic_walk <- function(bags,
                              tree_class,
                              misslabeled_str,
                              bag_size,
                              max_per_bag,
                              max_misslablel = FALSE,
                              no_rounds) {
  # this will measure how many times a certain strain is repeated in all bags
  count <- setNames(rep(0, ape::Ntip(tree_class)), tree_class$tip.label)

  # first scale the tree branch lengths so they are between the values 1 and 1000
  # this is important for phylogenetic walk cuz it uses sqrt to reduce large probs
  # Also, if en edge has length 0 now it will have some non zero length

  min_val <- tree_class$edge.length[which.min(tree_class$edge.length)]
  max_val <- tree_class$edge.length[which.max(tree_class$edge.length)]
  tree_class$edge.length <- (tree_class$edge.length - min_val) / (max_val - min_val) * 999 + 1
  # in rare cases all the edge lengts could be the same which then returns a set of NaN values as the edge lengths
  # give the tree a uniform length
  if (any(is.nan(tree_class$edge.length))) {tree_class$edge.length <- rep(100, length(tree_class$edge.length))}

  tree_class_save <- tree_class
  # the loop that makes the bags
  for (i in 1:no_rounds) {
    # refresh the tree
    tree_class <- tree_class_save

    # initialize new count for this bag
    count_bag <- setNames(rep(0, length(tree_class$tip.label)),
                          tree_class$tip.label)
    # this is here to ensure that the each genomes appears only "max_per_bag"
    # amount of times in the bag. it also ensures that misslabelled strains
    # does not overtake the non-misslabelled

    for (j in 1:bag_size) {
      rootie <- tree_class$edge[1,1]
      left_subtree <- tree_class$edge[,2][tree_class$edge[,1] == rootie][1]
      if (left_subtree <= length(tree_class$tip.label)) {
        # if the new subtree constains only a tip run this
        index <- which(tree_class$edge[,2] == left_subtree)
        left_tree_prob <- tree_class$edge.length[index]
        left_tree <- tree_class$tip.label[left_subtree]
      } else {
        # if the subtree is still large run this
        distance_from_root <- tree_class$edge.length[tree_class$edge[,1] == rootie][1]
        left_tree <- ape::extract.clade(tree_class, node = left_subtree)
        left_tree_prob <- sum(left_tree$edge.length) + distance_from_root
      }

      # do the same but with the right side of the tree
      right_subtree <- tree_class$edge[,2][tree_class$edge[,1] == rootie][2]
      if (right_subtree <= length(tree_class$tip.label)) {
        index <- which(tree_class$edge[,2] == right_subtree)
        right_tree_prob <- tree_class$edge.length[index]
        right_tree <- tree_class$tip.label[right_subtree]
      } else {
        distance_from_root <- tree_class$edge.length[tree_class$edge[,1] == rootie][2]
        right_tree <- ape::extract.clade(tree_class, node = right_subtree)
        right_tree_prob <- sum(right_tree$edge.length) + distance_from_root
      }

      repeat {
        # decide if the next step will be on the right of left of the tree
        step <- base::sample(c("left","right"), size = 1, prob = c(left_tree_prob, right_tree_prob))
        if (step == "left") {
          if (is.character(left_tree) == TRUE) {
            # add that to the bag
            bags[[i]][length(bags[[i]]) + 1] <- left_tree

            # count the sample into the two counters
            index_count <- which(names(count) == left_tree)
            count[index_count] <- count[index_count] + 1
            count_bag_index <- which(names(count_bag) == left_tree)
            count_bag[count_bag_index] <- count_bag[count_bag_index] + 1

            # remove those samples that are in bag more than max_per_bag times
            if (count_bag[count_bag_index] == max_per_bag) {
              # remove the strain from the tree
              tree_class <- ape::drop.tip(tree_class, tip = left_tree)
            }

            # remove misslabelled strains that are in bag more than max_misslablel times
            if (is.numeric(max_misslablel)) {
              name <- names(count_bag)[count_bag_index]
              if (count_bag[count_bag_index] == max_misslablel &
                  any(misslabeled_str == name)) {
                tree_class <- ape::drop.tip(tree_class, tip = name)
              }
            }
            break
          }

          # if we did not yet reached the tip, continue
          tree_tmp <- left_tree
          rootie <- tree_tmp$edge[1,1]

          left_subtree <- tree_tmp$edge[,2][tree_tmp$edge[,1] == rootie][1]
          if (left_subtree <= length(tree_tmp$tip.label)) {
            # if the new subtree constains only a tip run this
            index <- which(tree_tmp$edge[,2] == left_subtree)
            left_tree_prob <- tree_tmp$edge.length[index]
            left_tree <- tree_tmp$tip.label[left_subtree]
          } else {
            # if the subtree is still large run this
            distance_from_root <- tree_tmp$edge.length[tree_tmp$edge[,1] == rootie][1]
            left_tree <- ape::extract.clade(tree_tmp, node = left_subtree)
            left_tree_prob <- sum(left_tree$edge.length) + distance_from_root
          }

          # do the same but with the right side of the tree
          right_subtree <- tree_tmp$edge[,2][tree_tmp$edge[,1] == rootie][2]
          if (right_subtree <= length(tree_tmp$tip.label)) {
            index <- which(tree_tmp$edge[,2] == right_subtree)
            right_tree_prob <- tree_tmp$edge.length[index]
            right_tree <- tree_tmp$tip.label[right_subtree]
          } else {
            distance_from_root <- tree_tmp$edge.length[tree_tmp$edge[,1] == rootie][2]
            right_tree <- ape::extract.clade(tree_tmp, node = right_subtree)
            right_tree_prob <- sum(right_tree$edge.length) + distance_from_root
          }

        } else if (step == "right") {
          # if right_tree is character it means we already reached the tip of the tree
          if (is.character(right_tree) == TRUE) {
            # add that to the bag
            bags[[i]][length(bags[[i]]) + 1] <- right_tree

            # count the sample into the two counters
            index_count <- which(names(count) == right_tree)
            count[index_count] <- count[index_count] + 1
            count_bag_index <- which(names(count_bag) == right_tree)
            count_bag[count_bag_index] <- count_bag[count_bag_index] + 1

            # remove those samples that are in bag more than max_per_bag times
            if (count_bag[count_bag_index] == max_per_bag) {
              # remove the strain from the tree
              tree_class <- ape::drop.tip(tree_class, tip = right_tree)
            }

            # remove misslabelled strains that are in bag more than max_misslablel times
            if (is.numeric(max_misslablel)) {
              name <- names(count_bag)[count_bag_index]
              if (count_bag[count_bag_index] == max_misslablel &
                  any(misslabeled_str == name)) {
                tree_class <- ape::drop.tip(tree_class, tip = right_tree)
              }
            }
            break
          }

          # if we did not yet reached the tip, continue
          tree_tmp <- right_tree
          rootie <- tree_tmp$edge[1,1]

          left_subtree <- tree_tmp$edge[,2][tree_tmp$edge[,1] == rootie][1]
          if (left_subtree <= length(tree_tmp$tip.label)) {
            # if the new subtree constains only a tip run this
            index <- which(tree_tmp$edge[,2] == left_subtree)
            left_tree_prob <- tree_tmp$edge.length[index]
            left_tree <- tree_tmp$tip.label[left_subtree]
          } else {
            # if the subtree is still large run this
            distance_from_root <- tree_tmp$edge.length[tree_tmp$edge[,1] == rootie][1]
            left_tree <- ape::extract.clade(tree_tmp, node = left_subtree)
            left_tree_prob <- sum(left_tree$edge.length) + distance_from_root
          }

          # do the same but with the right side of the tree
          right_subtree <- tree_tmp$edge[,2][tree_tmp$edge[,1] == rootie][2]
          if (right_subtree <= length(tree_tmp$tip.label)) {
            index <- which(tree_tmp$edge[,2] == right_subtree)
            right_tree_prob <- tree_tmp$edge.length[index]
            right_tree <- tree_tmp$tip.label[right_subtree]
          } else {
            distance_from_root <- tree_tmp$edge.length[tree_tmp$edge[,1] == rootie][2]
            right_tree <- ape::extract.clade(tree_tmp, node = right_subtree)
            right_tree_prob <- sum(right_tree$edge.length) + distance_from_root
          }
        }
      }
    }
    # its possible that phylo walk did not select the misslabelled strians. If thats
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
