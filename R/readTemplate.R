#' Reads a template file for labeling of individual wells 
#' 
#' @param templateFile The template file to be read.
#' @param type The type of template being read. Can be either "strains", "conditions", "controls."
#' @return Returns a data frame with 3 columns corresponing to the row, column, and value for either of the three types.
#' @export


readTemplate <- function(templateFile, type){
    template <- read.csv(templateFile, check.names=FALSE)
    meltTemplate <- melt(template, id.vars = "row")
    if(type=="strains"){
        colnames(meltTemplate) <- c("row", "col", "strain")
    } else if(type=="conditions"){
        colnames(meltTemplate) <- c("row", "col", "condition")
    } else if(type=="controls"){
        colnames(meltTemplate) <- c("row", "col", "control")
    }
    return(meltTemplate)
}