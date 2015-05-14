#' summarizeSorted
#'
#'
#' @param plates A list consisting of all of the score plates in the first element and the setup plates in the second elements.
#' @param strains A character vector of length equal to the number of wells in in the plate being summarized. Defaults to NULL.
#' @param quantiles Boolean indicating whether or not quantile values (q10, q25, q75, q90) should be included in the summarized data frame. Defaults to FALSE.
#' @param log Boolean indicating whether or not log-transformed values should be included in the summarized data frame. Defaults to FALSE.
#' @param ends Boolean indicating whether or not min and max values should be included in the summarized data frame. Defaults to FALSE.
#' @export

summarizeSorted <- function(plates, strains=NULL, quantiles=FALSE, log=FALSE, ends=FALSE){
    score <- summarizePlate_worms(plates[[1]], strains, quantiles, log, ends) %>% arrange(plate, row, col)
    setup <- plates[[2]] %>% group_by(date, experiment, round, assay, plate, drug, row, col) %>% summarize(control.n.sorted = sum(sort==6))
    data <- left_join(score, setup) %>% mutate(norm.n = n/control.n.sorted) %>% select(-n.sorted, -control.n.sorted)
    return(data)
}