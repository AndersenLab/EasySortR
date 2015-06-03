#' Regress control values and, optionally, assay values from experimental values
#'
#'
#' @param dataframe A data frame in long form that contains all of the experimental and control values.
#' @param assay Boolean stating whether to regress out assay as well as controls.
#' @return A data frame in long form
#' @importFrom dplyr %>%
#' @export

regress <- function(dataframe, assay=FALSE){
    dataframe <- ensure_long(dataframe)
    dataframe <- dataframe %>%
        filter(trait != "n.sorted")
    data <- dataframe %>%
            dplyr::filter(!is.na(control))
    controls <- dataframe %>%
                dplyr::filter(is.na(control))
    controls$control <- controls$condition
    moltendata <- data
    moltencontrols <- controls %>% dplyr::group_by(strain, control, trait) %>%
        dplyr::summarize(controlphenotype = mean(phenotype, na.rm=TRUE))
    fusedmoltendata <- dplyr::left_join(moltendata, moltencontrols,
                                 by=c("strain", "control", "trait"))
    if(assay){
        modeldata <- fusedmoltendata %>%
            dplyr::group_by(condition, trait) %>%
            dplyr::filter(!is.na(phenotype), !is.na(controlphenotype)) %>%
            dplyr::do(broom::augment(lm(phenotype ~ controlphenotype + assay,
                        data=.)))
    } else {
        modeldata <- fusedmoltendata %>%
            dplyr::group_by(condition, trait) %>%
            dplyr::filter(!is.na(phenotype), !is.na(controlphenotype)) %>%
            dplyr::do(broom::augment(lm(phenotype ~ controlphenotype, data=.)))
    }
    resids <- modeldata %>%
        dplyr::select(condition, trait, phenotype, controlphenotype, .resid) %>%
        dplyr::rename(resid = .resid)
    arrangedmolten <- arrange(fusedmoltendata, condition, trait, phenotype,
                              controlphenotype)
    arrangedresids <- arrange(resids, condition, trait, phenotype,
                              controlphenotype)
    
    
    if (all.equal(arrangedmolten[,c(5,11,12,16)], arrangedresids[,-5])) {
        regressedframe <- cbind(arrangedmolten, arrangedresids$resid)
    } else {
        stop("The values for the raw data (when arranged on 'condition', 'trait', 'phenotype', and 'controlphenotype' do not match up with same values in the residual data frame output by the broom package's augment function. Quitting to prevent binding of incorrect data")
    }
    
    return(regressedframe)
}