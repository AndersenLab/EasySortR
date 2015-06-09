#' Regress control values and, optionally, assay values from experimental values
#'
#'
#' @param dataframe A data frame in long form that contains all of the experimental and control values.
#' @param assay Boolean stating whether to regress out assay as well as controls.
#' @return A data frame in long form
#' @importFrom dplyr %>%
#' @export

regress <- function(dataframe, assay=FALSE){
    
    # Ensure the data frame is in long format, then filter out all NA strains
    # and nonfinite (NA, NaN, Inf) values
    
    dataframe <- ensure_long(dataframe)
    dataframe <- dplyr::filter(dataframe, is.finite(phenotype), !is.na(strain))
    
    # Get rid of n.sorted phenotype values (always 0)
    
    dataframe <- dataframe %>%
        dplyr::filter(trait != "n.sorted")
    data <- dataframe %>%
        
    # Separate the data from the controls into two different data frames
        
            dplyr::filter(!is.na(control))
    controls <- dataframe %>%
                dplyr::filter(is.na(control))
    controls$control <- controls$condition
    moltendata <- data
    
    # Summarize the controls in case there are replicate controls
    
    moltencontrols <- controls %>% dplyr::group_by(strain, control, trait) %>%
        dplyr::summarize(controlphenotype = mean(phenotype, na.rm=TRUE))
    
    # Join the control values back to the controls
    
    fusedmoltendata <- dplyr::left_join(moltendata, moltencontrols,
                                 by=c("strain", "control", "trait")) %>%
        dplyr::filter(!is.na(phenotype), !is.na(controlphenotype))
    
    # Perform the regression step with or without the assay values and pull out
    # the residuals with the broom package
    
    if(assay){
        modeldata <- fusedmoltendata %>%
            dplyr::group_by(condition, trait) %>%
            dplyr::do(broom::augment(lm(phenotype ~ controlphenotype + assay,
                        data=.)))
    } else {
        modeldata <- fusedmoltendata %>%
            dplyr::group_by(condition, trait) %>%
            dplyr::do(broom::augment(lm(phenotype ~ controlphenotype, data=.)))
    }
    
    # Select only the information that we're interested in and rename .resid to
    # just resid
    
    resids <- modeldata %>%
        dplyr::select(condition, trait, phenotype, controlphenotype, .resid) %>%
        dplyr::rename(resid = .resid)
    
    # Arrange the columns so that they align for a cbind
    
    arrangedmolten <- dplyr::arrange(fusedmoltendata, condition, trait, phenotype,
                              controlphenotype)
    arrangedresids <- dplyr::arrange(resids, condition, trait, phenotype,
                              controlphenotype)
    
    # Make sure that the number of rows between the two arranged data frames are
    # the same and that the values for the matching columns are all the same
    
    tryCatch({
        if (all.equal(arrangedmolten[,c("condition", "trait", "phenotype",
                                        "controlphenotype")],
                      arrangedresids[,c("condition", "trait", "phenotype",
                                        "controlphenotype")])) {
            regressedframe <- cbind(arrangedmolten, arrangedresids$resid) %>%
                dplyr::rename(resid = `arrangedresids$resid`)
        } else {
            stop("The values for the raw data (when arranged on 'condition', 'trait', 'phenotype', and 'controlphenotype' do not match up with same values in the residual data frame output by the broom package's augment function. Quitting to prevent binding of incorrect data")
        }
    }, error = function(e) stop("The number of rows in arrangedmolten and arrangedresids do not match."))
    
    return(regressedframe)
}