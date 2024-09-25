#' Split features that were grouped
#'
#' @noRd
split_gfs <- function(df, pangenome = FALSE, DRAM = FALSE) {
  index <- c()
  gfs_names_vect <- c()
  for (i in 1:nrow(df)) {
    if (!grepl("___", df[i,1])) {next}
    # split the old gfs name
    gfs_names <- strsplit(df[i, 1], split = "___")[[1]]
    index <- c(index, rep(i, length(gfs_names)))
    gfs_names_vect <- c(gfs_names_vect, gfs_names)
  }

  missing <- setdiff(1:nrow(df), index)
  names(missing) <- df[missing, 1]
  names(index) <- gfs_names_vect
  all_indexes <- c(missing, index)
  all_indexes <- all_indexes[order(all_indexes)]

  df <- df[all_indexes, ]
  df[,1] <- names(all_indexes)
  if (pangenome == TRUE) {
    # replace "..." to "~~~"
    df[,1] <- gsub("\\.", "~", df[,1])
  } else if (DRAM == TRUE) {
    # correct the gfs names if the user is using DRAM input
    df[,1] <- gsub("~", " ", df[,1])
    df[,1] <- gsub(" Woodcroft ", "(Woodcroft)", df[,1])
  }

  return(df)
}
