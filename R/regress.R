#' Regress control values and, optionally, assay values from experimental values
#'
#'
#' @param dataFrame A data frame in long form that contains all of the experimental and control values.
#' @param assay Boolean stating whether to regress out assay as well as controls.
#' @return A data frame in long form 
#' @import broom
#' @import dplyr
#' @import tidyr
#' @export

regress <- function(dataFrame, assay=FALSE){
    dataFrame <- ensureLong(dataFrame)
    data <- dataFrame %>% dplyr::filter(!is.na(control))
    controls <- dataFrame %>% dplyr::filter(is.na(control))
    controls$control <- controls$condition
    moltenData <- data
    moltenControls <- controls %>% group_by(strain, control, trait) %>% summarize(controlPhenotype = mean(phenotype, na.rm=TRUE))
    fusedMoltenData <- left_join(moltenData, moltenControls, by=c("strain", "control", "trait"))
    if(assay){
        modelData <- fusedMoltenData %>% group_by(condition, trait) %>% filter(!is.na(phenotype), !is.na(controlPhenotype)) %>% do(broom::augment(lm(phenotype~controlPhenotype+assay, data=.)))
    } else {
        modelData <- fusedMoltenData %>% group_by(condition, trait) %>% filter(!is.na(phenotype), !is.na(controlPhenotype)) %>% do(broom::augment(lm(phenotype~controlPhenotype, data=.)))
    }
    resids <- modelData %>% select(condition, trait, phenotype, controlPhenotype, .resid) %>% rename(resid = .resid)
    regressedFrame <- left_join(fusedMoltenData, resids) %>% group_by(row, col, trait) %>% distinct()
    return(regressedFrame)
}