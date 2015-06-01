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

ensure_wide <- function(data){
    if("n" %in% colnames(data)){
        return(data)
    } else {
        widedata <- tidyr::spread(data, trait, phenotype)
        return(widedata)
    }
}