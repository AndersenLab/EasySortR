#' Read all of the plate files from a specific directory
#'
#' This function reads in all of the plates within a single directory.
#'
#' @param directory The directory to be read
#' @param tofmin The minimum time of flight value allowed. Defaults to 60.
#' @param tofmin The minimum time of flight value allowed. Defaults to 2000.
#' @param extmin The minimum extinction value allowed. Defaults to 0.
#' @param extmin The maximum extinction value allowed. Defaults to 10000.
#' @param SVM Boolean specifying whether or not to use the support vector machine to separate worms and bubbles.
#' @return A single data frame with all of the plates from the read directory.
#' @export

readDirectory <- function(directory, tofmin=60, tofmax=2000, extmin=0, extmax=10000, SVM=TRUE) {
    dirFiles <- dir(directory)
    files <- subset(dirFiles, grepl(".txt$", dirFiles))
    filePaths <- file.path(directory, files)
    plates <- lapply(filePaths, function(x){readFile(x, tofmin, tofmax, extmin, extmax, SVM)})
    data <- do.call(rbind, plates)
    return(data)
}