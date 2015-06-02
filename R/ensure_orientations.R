# Ensure that the data is in long format using tidyr

ensure_long <- function(data){
    if("trait" %in% colnames(data)){
        return(data)
    } else {
        longdata <- tidyr::gather(data, trait, phenotype, -c(date, experiment,
                                  round, assay, plate, condition, control,
                                  strain, row, col))
        return(longdata)
    }
}

# Ensure that the data is in wide format using tidyr

ensure_wide <- function(data){
    if("n" %in% colnames(data)){
        return(data)
    } else {
        widedata <- tidyr::spread(data, trait, phenotype)
        return(widedata)
    }
}