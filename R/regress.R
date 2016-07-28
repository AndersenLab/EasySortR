#' Regress control values or assay values from experimental values
#'
#'
#' @param dataframe A data frame in long form that contains all of the experimental and control values.
#' @param assay Boolean stating whether to regress out assay or controls. Assay
#' if \code{TRUE}, controls if \code{FALSE}. Defaults to \code{FALSE}.
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

    # Perform the regression step with or without the assay values and pull out
    # the residuals with the broom package

    if(assay){

        # Remove condition only present in one assay, otherwise regression will
        # fail

        dataframe <- dataframe %>%
            dplyr::group_by(condition) %>%
            dplyr::filter(length(unique(assay)) > 1)

        # Get the model data
        # `deparse` and `substitute` are necessary to get the formula to work
        # inside of a call to `mutate`, see:
        # http://stackoverflow.com/questions/28776792/combining-dplyrdo-with-dplyrmutate

        regressed <- fusedmoltendata %>%
          dplyr::group_by(condition, trait) %>%
          do(fit = lm(phenotype ~ assay - 1, .))
        
        withresids <- regressed%>%
          broom::augment(fit)%>%
          ungroup()%>%
          left_join(fusedmoltendata,.,by=c("condition", "trait", "phenotype", "controlphenotype"))%>%
          distinct(condition, trait, phenotype, controlphenotype,strain,row,col,plate,.keep_all = T)%>%
          rename(resid = .resid)
        
        regressedframe <- withresids %>%
          dplyr::mutate(phenotype = resid) %>%
          dplyr::select(-resid, -controlphenotype)

    } else {
        # Separate the data from the controls into two different data frames

        data <- dataframe %>%
            dplyr::filter(!is.na(control))
        controls <- dataframe %>%
            dplyr::filter(is.na(control) | control == "None")
        controls$control <- controls$condition
        moltendata <- data

        # Summarize the controls in case there are replicate controls

        moltencontrols <- controls %>% dplyr::group_by(strain, control, trait, assay) %>%
            dplyr::summarize(controlphenotype = mean(phenotype, na.rm=TRUE))

        # Join the control values back to the controls

        fusedmoltendata <- dplyr::left_join(moltendata, moltencontrols,
                                            by=c("strain", "control", "trait", "assay")) %>%
            dplyr::filter(!is.na(phenotype), !is.na(controlphenotype))

        # Get the model data
        # `deparse` and `substitute` are necessary to get the formula to work
        # inside of a call to `mutate`, see:
        # http://stackoverflow.com/questions/28776792/combining-dplyrdo-with-dplyrmutate

        regressed <- fusedmoltendata %>%
          dplyr::group_by(condition, trait) %>%
          do(fit = lm(phenotype ~ controlphenotype - 1, .))
        
        withresids <- regressed%>%
          broom::augment(fit)%>%
          ungroup()%>%
          left_join(fusedmoltendata,.,by=c("condition", "trait", "phenotype", "controlphenotype"))%>%
          distinct(condition, trait, phenotype, controlphenotype,strain,row,col,plate,.keep_all = T)%>%
          rename(resid = .resid)
        
        regressedframe <- withresids %>%
          dplyr::mutate(phenotype = resid) %>%
          dplyr::select(-resid, -controlphenotype)

    }

    return(regressedframe)
}