#' Regress control values and, optionally, assay values from experimental values
#'
#'
#' @param dataFrame A data frame in long form that contains all of the experimental and control values.
#' @param assay Boolean stating whether to regress out assay as well as controls.
#' @return A data frame in long form 
#' @import dplyr
#' @import reshape2
#' @import broom
#' @export

regress <- function(dataFrame, assay=FALSE){
    data <- dataFrame %>% dplyr::filter(!is.na(control))
    controls <- dataFrame %>% dplyr::filter(is.na(control))
    controls$control <- controls$condition
    moltenData <- reshape2::melt(data, value.name="phenotype", id.vars=c("date", "experiment", "round", "assay", "plate", "condition", "control", "strain", "row", "col"))
    moltenControls <- reshape2::melt(controls, value.name="controlPhenotype", id.vars=c("date", "experiment", "round", "assay", "plate", "condition", "control", "strain", "row", "col")) %>% select(strain, control, variable, controlPhenotype) %>% group_by(strain, control, variable) %>% summarize(controlPhenotype = mean(controlPhenotype, na.rm=TRUE))
    fusedMoltenData <- left_join(moltenData, moltenControls, by=c("strain", "control", "variable")) %>% rename(trait = variable)
    if(assay){
        modelData <- fusedMoltenData %>% group_by(condition, trait) %>% filter(!is.na(phenotype), !is.na(controlPhenotype)) %>% do(broom::augment(lm(phenotype~controlPhenotype+assay, data=.)))
    } else {
        modelData <- fusedMoltenData %>% group_by(condition, trait) %>% filter(!is.na(phenotype), !is.na(controlPhenotype)) %>% do(broom::augment(lm(phenotype~controlPhenotype, data=.)))
    }
    resids <- modelData %>% select(condition, trait, phenotype, controlPhenotype, .resid) %>% rename(resid = .resid)
    regressedFrame <- left_join(fusedMoltenData, resids) %>% group_by(row, col, trait) %>% distinct()
    return(regressedFrame)
}