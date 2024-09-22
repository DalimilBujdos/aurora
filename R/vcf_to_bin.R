#' Retrieves data from vcf file
#'
#' @noRd
vcf_to_bin <- function(bin_mat_snps, pheno_mat, which_snps) {
  custom_skip <- function(line) {
    substr(line, 1, 2) == "##"
  }

  # Read lines from the file, skipping lines starting with ##
  filtered_lines <- readLines(bin_mat_snps, warn = FALSE)
  filtered_lines <- filtered_lines[sapply(filtered_lines, custom_skip)]

  data <- read.table(bin_mat_snps,
                     sep = "\t",
                     quote = "",
                     comment.char = "",
                     colClasses = "character",
                     header = TRUE,
                     skip = length(filtered_lines))
  rownames(data) <- paste(data[,1], data[,2], sep = "_")
  # remove all snips where we are not sure about the alelle
  data <- data[!grepl("[*N<>]", data$ALT), ]
  bin_mat <- data[,-1:-9]
  index_remove_rows <- which(grepl(",", data$ALT))
  remove_rows <- data[index_remove_rows, ]
  bin_mat <- bin_mat[-index_remove_rows, ]
  if (which_snps == "all_alleles") {
    new_rows <- function(index) {
      row <- remove_rows[index, ]
      no_new_rows <- sum(unlist(gregexpr(",", row$ALT)) >= 0) + 1
      new_names <- paste(rep(rownames(row), no_new_rows),
                         strsplit(row$ALT, ",")[[1]], sep = "_")
      row <- as.numeric(row[1, -1:-9])

      # Find unique values in the vector
      unique_values <- unique(row)[unique(row) != 0]
      unique_values <- order(unique_values, decreasing = F)
      # Initialize empty binary vectors for each unique value
      binary_vectors <- lapply(unique_values, function(unique_value) {
        # Create a binary vector indicating the presence of the unique value
        binary_vector <- as.integer(row == unique_value)
        return(binary_vector)
      })
      return(list(vctrs = binary_vectors, names = new_names))
    }
    list_new_rows <- lapply(1:nrow(remove_rows), new_rows)
    df_new_rows <- do.call(rbind, lapply(list_new_rows, function(entry) {
      vectors <- unlist(entry$vctrs)
      row_names <- entry$names
      data.frame(matrix(vectors, nrow = length(row_names), byrow = TRUE), row.names = row_names)
    }))
    colnames(df_new_rows) <- colnames(bin_mat)
    bin_mat <- rbind(bin_mat, df_new_rows)
  }
  fut_rows <- rownames(bin_mat)
  bin_mat <- apply(bin_mat, 2, as.numeric)
  rownames(bin_mat) <- fut_rows
  bin_mat <- as.data.frame(t(bin_mat))
  # convert the colnames to a standardised format
  colnames(bin_mat) <- gsub("[^A-Za-z0-9_]", "_", colnames(bin_mat))
  colnames(bin_mat) <- gsub("^\\d", "x\\0", colnames(bin_mat))
  colnames(bin_mat) <- gsub("___", "_", colnames(bin_mat))
  bin_mat[match(pheno_mat$ids, rownames(bin_mat)), ]
  return(bin_mat)
}
