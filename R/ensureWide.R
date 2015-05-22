#' Ensure a data frame is in wide form
#'
#'
#' @param data A data frame
#' @return A data frame in wide form
#' @import tidyr
#' @export

ensureWide <- function(data){
    if("n" %in% colnames(data)){
        return(data)
    } else {
        wideData <- tidyr::spread(data, trait, phenotype)
        return(wideData)
    }
}