#' Regress control values and, optionally, assay values from experimental values
#'
#'
#' @param dataframe A data frame in long form that contains all of the experimental and control values.
#' @param assay Boolean stating whether to regress out assay as well as controls.
#' @return A data frame in long form
#' @import broom
#' @import dplyr
#' @import tidyr
#' @export

regress <- function(dataframe, assay=FALSE){
    dataframe <- ensure_long(dataframe)
    data <- dataframe %>%
            dplyr::filter(!is.na(control))
    controls <- dataframe %>%
                dplyr::filter(is.na(control))
    controls$control <- controls$condition
    moltendata <- data
    moltencontrols <- controls %>% group_by(strain, control, trait) %>%
                      summarize(controlphenotype = mean(phenotype, na.rm=TRUE))
    fusedmoltendata <- left_join(moltendata, moltencontrols,
                                 by=c("strain", "control", "trait"))
    if(assay){
        modeldata <- fusedmoltendata %>%
                     group_by(condition, trait) %>%
                     filter(!is.na(phenotype), !is.na(controlphenotype)) %>%
                     do(broom::augment(lm(phenotype ~ controlphenotype + assay,
                        data=.)))
    } else {
        modeldata <- fusedmoltendata %>%
                     group_by(condition, trait) %>%
                     filter(!is.na(phenotype), !is.na(controlphenotype)) %>%
                     do(broom::augment(
                        lm(phenotype ~ controlphenotype, data=.)))
    }
    resids <- modeldata %>%
              select(condition, trait, phenotype, controlphenotype, .resid) %>%
              rename(resid = .resid)
    regressedframe <- left_join(fusedmoltendata, resids) %>%
                      group_by(row, col, trait) %>%
                      distinct()
    return(regressedframe)
}