#' Function sets up conda environment for aurora
#'
#' @return path to the installed conda environment
#'
#' @export
#'
#' @examples
#' setup_conda()
setup_conda <- function() {
  tryCatch(
    {
      output <- system("conda info --envs", intern = TRUE)
    },
    error = function(e) {
      cat("Error message: You likely did not install conda. \n")
    }
  )
  if(any(grepl("aurora3", output))) {
    env_line <- grep("aurora3", output, value = TRUE)
    path <- gsub(".*?(C:\\S+).*", "\\1", env_line)
    path <- gsub("\\\\", "/", path)

    cat(path)
    return(path)
  } else {
    tryCatch(
      {
        # Run the system command
        result <- system(paste0("conda env create -f ", paste(system.file(package="aurora"), ".inst/aurora3.yaml", sep="/")), intern = TRUE)

        # Check if the result is 127
        if (result == 127) {
          cat("Error message: Command 'conda' not found or not executable.\n")
          cat("Error details:", result, "\n")
        } else if (result != 0) {
          cat("Error message: An error occurred while running the command.\n")
          cat("Error details:", result, "\n")
        }
      },
      error = function(e) {
        cat("Error message: An error occurred with the system() function.\n")
        cat("R Error message:", conditionMessage(e), "\n")
      }
    )
    output <- system("conda info --envs", intern = TRUE)
    env_line <- grep("aurora3", output, value = TRUE)
    path <- gsub(".*?(C:\\S+).*", "\\1", env_line)
    path <- gsub("\\\\", "/", path)

    cat(path)
    return(path)
  }
}
