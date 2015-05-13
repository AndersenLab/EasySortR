regress <- function(data, controls){
    moltenData <- melt(data, value.name="phenotype")
    moltenControls <- melt(controls, value.name="controlPhenotype")
    fusedMoltenData <- left_join(moltenData, moltenControls) %>% rename(trait = variable)
    modelData <- fusedMoltenData %>% group_by(trait) %>% filter(!is.na(phenotype), !is.na(controlPhenotype)) %>% do(augment(lm(.$phenotype ~ .$controlPhenotype, na.action=na.omit)))
    resids <- modelData %>% select(trait, ..phenotype, ..controlPhenotype, .resid) %>% rename(phenotype = ..phenotype, controlPhenotype = ..controlPhenotype, resid = .resid)
    regressedFrame <- left_join(fusedMoltenData, resids) %>% group_by(row, col, trait) %>% distinct() %>% %>% melt()
    test <- regressedFrame %>% select(-controlValue) %>% dcast(row+col~variable + resid)
    
    
    
    mdata <- data %>% gather(n, n.sorted, mean.TOF)
    
    
    
    
    
    
    
    
    
    
    regressedValues <- data.frame(do.call(cbind, lapply(which(colnames(data)=="n"):ncol(data),
                                                        function(x){
                                                            tryCatch({residuals(lm(data[,x] ~ data$assay + controlValues[,which(colnames(controlValues)==colnames(data)[x])], na.action=na.exclude))},
                                                                     error = function(err){return(NA)})
                                                        })))
    
    regressedAssayValues <- data.frame(do.call(cbind, lapply(which(colnames(data)=="n"):ncol(data),
                                                             function(x){
                                                                 tryCatch({residuals(lm(data[,x] ~ data$assay, na.action=na.exclude))},
                                                                          error = function(err){return(NA)})
                                                             })))
    
    
    
    #     reactValues <- data.frame(do.call(cbind, lapply(which(colnames(data)=="n"):ncol(data),
    #                                                     function(x){
    #                                                         reactNorms <- data[,x] - controlValues[,which(colnames(controlValues)==colnames(data)[x])]
    #                                                         if(length(reactNorms)==0){
    #                                                             reactNorms <- NA
    #                                                         }
    #                                                         return(reactNorms)
    #                                                     })))
    
    finalDF <- data.frame(data, regressedValues)
    colnames(finalDF)[(which(colnames(finalDF)=="norm.n")+1):ncol(finalDF)] <- paste0("resid.", colnames(finalDF)[which(colnames(finalDF)=="n"):which(colnames(finalDF)=="norm.n")])
    finalDF <- data.frame(finalDF, regressedAssayValues)
    colnames(finalDF)[(which(colnames(finalDF)=="resid.norm.n")+1):ncol(finalDF)] <- paste0("resid.a.", colnames(finalDF)[which(colnames(finalDF)=="n"):which(colnames(finalDF)=="norm.n")])
    return(finalDF)
}