#' Ensure a data frame is in long form
#'
#'
#' @param data A data frame
#' @return A data frame in long form
#' @import tidyr
#' @export

ensureLong <- function(data){
    if("trait" %in% colnames(data)){
        return(data)
    } else {
        longData <- tidyr::gather(data, trait, phenotype, -c(date, experiment, round, assay, plate, condition, control, strain, row, col))
        return(longData)
    }
}