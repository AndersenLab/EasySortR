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

        regressed <- dataframe %>%
          dplyr::group_by(condition, trait) %>%
          do(fit = lm(phenotype ~ assay - 1, .))
        
        withresids <- regressed%>%
          broom::augment(fit)%>%
          ungroup()%>%
          left_join(dataframe,.,by=c("condition", "trait", "phenotype", "assay"))%>%
          distinct(condition, trait, phenotype,strain,row,col,plate,.keep_all = T)%>%
          rename(resid = .resid)
        
        regressedframe <- withresids %>%
          dplyr::mutate(phenotype = resid) %>%
          dplyr::select(-resid, -.fitted, -.se.fit, -.hat, -.sigma, -.cooksd, -.std.resid)

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
          dplyr::do(fit = lm(phenotype ~ controlphenotype - 1, .))
        
        withresids <- regressed%>%
          broom::augment(fit)%>%
          dplyr::ungroup()%>%
          dplyr::left_join(fusedmoltendata,.,by=c("condition", "trait", "phenotype", "controlphenotype"))%>%
          dplyr::distinct(condition, trait, phenotype, controlphenotype,strain,row,col,plate,.keep_all = T)%>%
          dplyr::rename(resid = .resid)

        
        regressedframe <- withresids %>%
          dplyr::mutate(phenotype = resid) %>%
          dplyr::select(-resid, -controlphenotype)

    }

    return(regressedframe)
}


#' Select alternative method for control value regression
#'
#'
#' @param dataframe A data frame in long form that contains all of the experimental and control values.
#' @param method Type of regression desired. Either "deltapheno" (default), which takes each condition
#' data point subtracted from the mean control phenotype per strain, or "fracpheno" which takes each
#' condition data point divided by the mean control phenotype per strain.
#' @return A data frame in long form
#' @importFrom dplyr %>%
#' @export
#' 

altregress <- function(dataframe, method= "deltapheno"){
    # Ensure the data frame is in long format, then filter out all NA strains
    # and nonfinite (NA, NaN, Inf) values
    
    if(method %nin% c("deltapheno","fracpheno")) stop("Regression method not recognized")
    
    dataframe <- ensure_long(dataframe)
    dataframe <- dplyr::filter(dataframe, is.finite(phenotype), !is.na(strain))
    
    # Get rid of n.sorted phenotype values (always 0)
    
    dataframe <- dataframe %>%
        dplyr::filter(trait != "n.sorted")
    
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
    
    #Add columns with the delta phenotype (strain's control average - condition data point) and fraction
    #phenotype (phenotype data point / strain's control average)
    
    regression <- fusedmoltendata %>%
        dplyr::group_by(condition, trait, strain) %>%
        mutate(deltapheno = controlphenotype - phenotype, fracpheno = phenotype / controlphenotype)
    
    if(method == "deltapheno") {
        regressedframe <- regression %>%
            dplyr::mutate(phenotype = deltapheno) %>%
            dplyr::select(-deltapheno, -fracpheno, -controlphenotype)
    } else if(method == "fracpheno") {
        regressedframe <- regression %>%
            dplyr::mutate(phenotype = fracpheno) %>%
            dplyr::select(-deltapheno, -fracpheno, -controlphenotype)
    }
    
    
    return(regressedframe)
}