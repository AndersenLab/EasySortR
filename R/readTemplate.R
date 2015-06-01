#' Reads a template file for labeling of individual wells
#'
#' @param templateFile The template file to be read.
#' @param type The type of template being read. Can be either "strains",
#' "conditions", "controls."
#' @return Returns a data frame with 3 columns corresponing to the row, column,
#' and value for either of the three types.
#' @import tidyr
#' @export


read_template <- function(templatefile, type){
    template <- read.csv(templatefile, check.names=FALSE)
    melttemplate <- tidyr::gather(template, col, variable, -row)
    if(type == "strains"){
        colnames(melttemplate) <- c("row", "col", "strain")
    } else if (type == "conditions"){
        colnames(melttemplate) <- c("row", "col", "condition")
    } else if (type == "controls"){
        colnames(melttemplate) <- c("row", "col", "control")
    } else if (type == "contam") {
        colnames(melttemplate) <- c("row", "col", "contamination")
    }
    return(melttemplate)
}