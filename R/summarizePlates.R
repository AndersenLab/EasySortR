#' summarizePlates
#'
#' This function summarizes plates when the associated setup plates are not included in the analysis. The function takes the \code{plate} argument and returns a summarized version of the plate, with one well per line. The original plate data must be in wide format.
#'
#' @param plate The raw plate data frame to be summarized.
#' @param strains A character vector of length equal to the number of wells in in the plate being summarized.
#' @param quantiles Boolean indicating whether or not quantile values (q10, q25, q75, q90) should be included in the summarized data frame. Defaults to FALSE.
#' @param log Boolean indicating whether or not log-transformed values should be included in the summarized data frame. Defaults to FALSE.
#' @param ends Boolean indicating whether or not min and max values should be included in the summarized data frame. Defaults to FALSE.
#' @import dplyr
#' @import COPASutils
#' @export

summarizePlates <- function(plate, strains=NULL, quantiles=FALSE, log=FALSE, ends=FALSE) {
    plate <- plate[plate$call50=="object" | plate$TOF == -1 | is.na(plate$call50),]
    plate <- COPASutils::fillWells(plate)
    processed <- plate %>% group_by(date, experiment, round, assay, plate, condition, control, strain, row, col) %>% summarise(n=ifelse(length(TOF[!is.na(TOF)])==0, NA, length(TOF[!is.na(TOF)])),
                                                            n.sorted=sum(sort==6),
                                                            
                                                            mean.TOF=mean(TOF, na.rm=TRUE),
                                                            min.TOF=as.numeric(quantile(TOF, na.rm=TRUE)[1]),
                                                            q10.TOF=as.numeric(quantile(TOF, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.TOF=as.numeric(quantile(TOF, probs=0.25, na.rm=TRUE)[1]),
                                                            median.TOF=median(TOF, na.rm=TRUE),
                                                            q75.TOF=as.numeric(quantile(TOF, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.TOF=as.numeric(quantile(TOF, probs=0.90, na.rm=TRUE)[1]),
                                                            max.TOF=as.numeric(quantile(TOF, na.rm=TRUE)[5]),
                                                            
                                                            mean.EXT=mean(EXT, na.rm=TRUE),
                                                            min.EXT=as.numeric(quantile(EXT, na.rm=TRUE)[1]),
                                                            q10.EXT=as.numeric(quantile(EXT, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.EXT=as.numeric(quantile(EXT, probs=0.25, na.rm=TRUE)[1]),
                                                            median.EXT=median(EXT, na.rm=TRUE),
                                                            q75.EXT=as.numeric(quantile(EXT, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.EXT=as.numeric(quantile(EXT, probs=0.90, na.rm=TRUE)[1]),
                                                            max.EXT=as.numeric(quantile(EXT, na.rm=TRUE)[5]),
                                                            
                                                            mean.red=mean(red, na.rm=TRUE),
                                                            min.red=as.numeric(quantile(red, na.rm=TRUE)[1]),
                                                            q10.red=as.numeric(quantile(red, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.red=as.numeric(quantile(red, probs=0.25, na.rm=TRUE)[1]),
                                                            median.red=median(red, na.rm=TRUE),
                                                            q75.red=as.numeric(quantile(red, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.red=as.numeric(quantile(red, probs=0.9, na.rm=TRUE)[1]),
                                                            max.red=as.numeric(quantile(red, na.rm=TRUE)[5]),
                                                            
                                                            mean.green=mean(green, na.rm=TRUE),
                                                            min.green=as.numeric(quantile(green, na.rm=TRUE)[1]),
                                                            q10.green=as.numeric(quantile(green, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.green=as.numeric(quantile(green, probs=0.25, na.rm=TRUE)[1]),
                                                            median.green=median(green, na.rm=TRUE),
                                                            q75.green=as.numeric(quantile(green, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.green=as.numeric(quantile(green, probs=0.9, na.rm=TRUE)[1]),
                                                            max.green=as.numeric(quantile(green, na.rm=TRUE)[5]),
                                                            
                                                            mean.yellow=mean(yellow, na.rm=TRUE),
                                                            min.yellow=as.numeric(quantile(yellow, na.rm=TRUE)[1]),
                                                            q10.yellow=as.numeric(quantile(yellow, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.yellow=as.numeric(quantile(yellow, probs=0.25, na.rm=TRUE)[1]),
                                                            median.yellow=median(yellow, na.rm=TRUE),
                                                            q75.yellow=as.numeric(quantile(yellow, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.yellow=as.numeric(quantile(yellow, probs=0.9, na.rm=TRUE)[1]),
                                                            max.yellow=as.numeric(quantile(yellow, na.rm=TRUE)[5]),
                                                            
                                                            mean.norm.EXT=mean(norm.EXT, na.rm=TRUE),
                                                            min.norm.EXT=as.numeric(quantile(norm.EXT, na.rm=TRUE)[1]),
                                                            q10.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.EXT=median(norm.EXT, na.rm=TRUE),
                                                            q75.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.EXT=as.numeric(quantile(norm.EXT, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.EXT=as.numeric(quantile(norm.EXT, na.rm=TRUE)[5]),
                                                            
                                                            mean.norm.red=mean(norm.red, na.rm=TRUE),
                                                            min.norm.red=as.numeric(quantile(norm.red, na.rm=TRUE)[1]),
                                                            q10.norm.red=as.numeric(quantile(norm.red, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.red=as.numeric(quantile(norm.red, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.red=median(norm.red, na.rm=TRUE),
                                                            q75.norm.red=as.numeric(quantile(norm.red, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.red=as.numeric(quantile(norm.red, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.red=as.numeric(quantile(norm.red, na.rm=TRUE)[5]),
                                                            
                                                            mean.norm.green=mean(norm.green, na.rm=TRUE),
                                                            min.norm.green=as.numeric(quantile(norm.green, na.rm=TRUE)[1]),
                                                            q10.norm.green=as.numeric(quantile(norm.green, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.green=as.numeric(quantile(norm.green, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.green=median(norm.green, na.rm=TRUE),
                                                            q75.norm.green=as.numeric(quantile(norm.green, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.green=as.numeric(quantile(norm.green, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.green=as.numeric(quantile(norm.green, na.rm=TRUE)[5]),
                                                            
                                                            mean.norm.yellow=mean(norm.yellow, na.rm=TRUE),
                                                            min.norm.yellow=as.numeric(quantile(norm.yellow, na.rm=TRUE)[1]),
                                                            q10.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.1, na.rm=TRUE)[1]),
                                                            q25.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.25, na.rm=TRUE)[1]),
                                                            median.norm.yellow=median(norm.yellow, na.rm=TRUE),
                                                            q75.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.75, na.rm=TRUE)[1]),
                                                            q90.norm.yellow=as.numeric(quantile(norm.yellow, probs=0.9, na.rm=TRUE)[1]),
                                                            max.norm.yellow=as.numeric(quantile(norm.yellow, na.rm=TRUE)[5]),
                                                            
                                                            mean.log.EXT=mean(log(EXT), na.rm=TRUE),
                                                            min.log.EXT=as.numeric(quantile(log(EXT), na.rm=TRUE)[1]),
                                                            q10.log.EXT=as.numeric(quantile(log(EXT), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.EXT=as.numeric(quantile(log(EXT), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.EXT=median(log(EXT), na.rm=TRUE),
                                                            q75.log.EXT=as.numeric(quantile(log(EXT), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.EXT=as.numeric(quantile(log(EXT), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.EXT=as.numeric(quantile(log(EXT), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.red=mean(log(red), na.rm=TRUE),
                                                            min.log.red=as.numeric(quantile(log(red), na.rm=TRUE)[1]),
                                                            q10.log.red=as.numeric(quantile(log(red), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.red=as.numeric(quantile(log(red), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.red=median(log(red), na.rm=TRUE),
                                                            q75.log.red=as.numeric(quantile(log(red), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.red=as.numeric(quantile(log(red), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.red=as.numeric(quantile(log(red), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.green=mean(log(green), na.rm=TRUE),
                                                            min.log.green=as.numeric(quantile(log(green), na.rm=TRUE)[1]),
                                                            q10.log.green=as.numeric(quantile(log(green), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.green=as.numeric(quantile(log(green), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.green=median(log(red), na.rm=TRUE),
                                                            q75.log.green=as.numeric(quantile(log(green), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.green=as.numeric(quantile(log(green), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.green=as.numeric(quantile(log(green), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.yellow=mean(log(yellow), na.rm=TRUE),
                                                            min.log.yellow=as.numeric(quantile(log(yellow), na.rm=TRUE)[1]),
                                                            q10.log.yellow=as.numeric(quantile(log(yellow), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.yellow=as.numeric(quantile(log(yellow), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.yellow=median(log(yellow), na.rm=TRUE),
                                                            q75.log.yellow=as.numeric(quantile(log(yellow), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.yellow=as.numeric(quantile(log(yellow), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.yellow=as.numeric(quantile(log(yellow), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.EXT=mean(log(norm.EXT), na.rm=TRUE),
                                                            min.log.norm.EXT=as.numeric(quantile(log(norm.EXT), na.rm=TRUE)[1]),
                                                            q10.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.EXT=median(log(norm.EXT), na.rm=TRUE),
                                                            q75.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.EXT=as.numeric(quantile(log(norm.EXT), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.EXT=as.numeric(quantile(log(norm.EXT), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.red=mean(log(norm.red), na.rm=TRUE),
                                                            min.log.norm.red=as.numeric(quantile(log(norm.red), na.rm=TRUE)[1]),
                                                            q10.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.red=median(log(norm.red), na.rm=TRUE),
                                                            q75.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.red=as.numeric(quantile(log(norm.red), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.red=as.numeric(quantile(log(norm.red), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.green=mean(log(norm.green), na.rm=TRUE),
                                                            min.log.norm.green=as.numeric(quantile(log(norm.green), na.rm=TRUE)[1]),
                                                            q10.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.green=median(log(norm.red), na.rm=TRUE),
                                                            q75.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.green=as.numeric(quantile(log(norm.green), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.green=as.numeric(quantile(log(norm.green), na.rm=TRUE)[5]),
                                                            
                                                            mean.log.norm.yellow=mean(log(norm.yellow), na.rm=TRUE),
                                                            min.log.norm.yellow=as.numeric(quantile(log(norm.yellow), na.rm=TRUE)[1]),
                                                            q10.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.1, na.rm=TRUE)[1]),
                                                            q25.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.25, na.rm=TRUE)[1]),
                                                            median.log.norm.yellow=median(log(norm.yellow), na.rm=TRUE),
                                                            q75.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.75, na.rm=TRUE)[1]),
                                                            q90.log.norm.yellow=as.numeric(quantile(log(norm.yellow), probs=0.90, na.rm=TRUE)[1]),
                                                            max.log.norm.yellow=as.numeric(quantile(log(norm.yellow), na.rm=TRUE)[5]),
                                                            
                                                            var.TOF=var(TOF),
                                                            cv.TOF=(sd(TOF,na.rm=TRUE)/mean(TOF,na.rm=TRUE)),
                                                            iqr.TOF=quantile(TOF, na.rm=TRUE, probs=.75)-quantile(TOF, na.rm=TRUE, probs=.25),
                                                            var.EXT=var(EXT),
                                                            cv.EXT=(sd(EXT,na.rm=TRUE)/mean(EXT,na.rm=TRUE)),
                                                            iqr.EXT=quantile(EXT, na.rm=TRUE, probs=.75)-quantile(EXT, na.rm=TRUE, probs=.25),
                                                            var.red=var(red),
                                                            cv.red=(sd(red,na.rm=TRUE)/mean(red,na.rm=TRUE)),
                                                            iqr.red=quantile(red, na.rm=TRUE, probs=.75)-quantile(red, na.rm=TRUE, probs=.25),
                                                            var.green=var(green),
                                                            cv.green=(sd(green,na.rm=TRUE)/mean(green,na.rm=TRUE)),
                                                            iqr.green=quantile(green, na.rm=TRUE, probs=.75)-quantile(green, na.rm=TRUE, probs=.25),
                                                            var.yellow=var(yellow),
                                                            cv.yellow=(sd(yellow,na.rm=TRUE)/mean(yellow,na.rm=TRUE)),
                                                            iqr.yellow=quantile(yellow, na.rm=TRUE, probs=.75)-quantile(yellow, na.rm=TRUE, probs=.25),
                                                            
                                                            f.L1 = length(which(stage == "L1"))/length(stage),
                                                            f.L2L3 = length(which(stage == "L2/L3"))/length(stage),
                                                            f.L4 = length(which(stage == "L4"))/length(stage),
                                                            f.ad = length(which(stage == "adult"))/length(stage))
    
    if(!ends){
        processed <- processed[,-(grep("min", colnames(processed)))]
        processed <- processed[,-(grep("max", colnames(processed)))]
    }
    if(!quantiles){
        processed <- processed[,-(grep("q", colnames(processed)))]
    }
    if(!log){
        processed <- processed[,-(grep("log", colnames(processed)))]
    }
    if(is.null(strains)){
        analysis <- processed
        analysis <- analysis[order(analysis$row, analysis$col),]
    } else {
        analysis <- data.frame(strain = as.character(strains), processed)
        analysis <- analysis[order(analysis$strain),]
        analysis <- analysis[order(analysis$row, analysis$col),]
    }
    analysis[analysis$mean.TOF==-1 | is.na(analysis$mean.TOF),which(colnames(analysis)=="n"):ncol(analysis)] <- NA
    return(analysis)
}