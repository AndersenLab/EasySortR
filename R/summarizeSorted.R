#' Summarize a plate and a add a column for the normalized population
#'
#' This function takes a list of plates as input and returns a single data frame
#' with data summarized for each well/condition pair. As input, this function
#' takes a two element list, with the first element being a data frame
#' consisting of all of the score plates to be summarized. The second element
#' in the list should be a data frame consisting of all of the setup plates to
#' be summarized, from which the number sorted to each well will be extracted.
#'
#' @param plates A list consisting of all of the score plates in the first
#' element and the setup plates in the second elements.
#' @param strains A character vector of length equal to the number of wells in
#' the plate being summarized. Defaults to NULL.
#' @param quantiles Boolean indicating whether or not quantile values (q10, q25,
#' q75, q90) should be included in the summarized data frame. Defaults to FALSE.
#' @param log Boolean indicating whether or not log-transformed values should be
#' included in the summarized data frame. Defaults to FALSE.
#' @param ends Boolean indicating whether or not min and max values should be
#' included in the summarized data frame. Defaults to FALSE.
#' @export

summarize_sorted <- function(plates, strains=NULL, quantiles=FALSE, log=FALSE,
                            ends=FALSE){
    score <- summarize_plates(plates[[1]], strains, quantiles, log, ends) %>%
        dplyr::arrange(plate, row, col)
    setup <- plates[[2]] %>%
        dplyr::group_by(date, experiment, round, assay,
                        plate, drug, row, col) %>% dplyr::summarize(
                            control.n.sorted = sum(sort == 6))
    data <- dplyr::left_join(score, setup) %>%
        dplyr::mutate(norm.n = n / control.n.sorted) %>%
        dplyr::select(-n.sorted, -control.n.sorted)
    return(data)
}