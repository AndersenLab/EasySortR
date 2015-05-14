#' readfile
#'
#' Reads sorter data files into R. Built for use exlusively with worms. Depends on COPASutils for the readSorter function.
#' @param file The file to be read.
#' @param tofmin The minimum time of flight value allowed. Defaults to 60.
#' @param tofmin The minimum time of flight value allowed. Defaults to 2000.
#' @param extmin The minimum extinction value allowed. Defaults to 0.
#' @param extmin The maximum extinction value allowed. Defaults to 10000.
#' @param SVM Boolean specifying whether or not to use the support vector machine to separate worms and bubbles.
#' @return Returns a single data frame for a single plate file.
#' @export

readFile <- function(file, tofmin=60, tofmax=2000, extmin=0, extmax=10000, SVM=TRUE, levels=2){
    plate <- COPASutils::readSorter(file, tofmin, tofmax, extmin, extmax)
    modplate <- with(plate,
                     data.frame(row=Row,
                                col=as.factor(Column),
                                sort=Status.sort,
                                TOF=TOF,
                                EXT=EXT,
                                time=Time.Stamp,
                                green=Green,
                                yellow=Yellow,
                                red=Red))
    modplate <- modplate %>%
        dplyr::group_by(row, col) %>%
        dplyr::do(COPASutils::extractTime(.))
    modplate <- data.frame(modplate)
    modplate[,10:13] <- apply(modplate[,c(5, 7:9)], 2, function(x){x/modplate$TOF})
    colnames(modplate)[10:13] <- c("norm.EXT", "norm.green", "norm.yellow", "norm.red")
    if(SVM){
        plateprediction <- kernlab::predict(COPASutils::bubbleSVMmodel_noProfiler, modplate[,3:length(modplate)], type="probabilities")
        modplate$object <- plateprediction[,"1"]
        modplate$call50 <- factor(as.numeric(modplate$object>0.5), levels=c(1,0), labels=c("object", "bubble"))
    }
    modplate$stage <- ifelse(modplate$TOF>=60 & modplate$TOF<90, "L1", 
                             ifelse(modplate$TOF>=90 & modplate$TOF<200, "L2/L3",
                                    ifelse(modplate$TOF>=200 & modplate$TOF<300, "L4",
                                           ifelse(modplate$TOF>=300, "adult", NA))))
    modplate[,as.vector(which(lapply(modplate, class) == "integer"))] <- lapply(modplate[,as.vector(which(lapply(modplate, class) == "integer"))], as.numeric)
    modplate <- cbind(info(file, levels), modplate)
    return(modplate)
}