#' Remove outliers by multivariate pruning
#'
#' 
#' @import mvoutlier
#' @export

mvPrune <- function(data, mvCols = c("f.L1", "f.ad", "mean.TOF", "mean.EXT")){
    data <- ensureWide(data)
    noNAs <- data %>% ungroup() %>% na.omit()
    outlierFrame <- noNAs %>% select_(.dots=mvCols)
    outlierFrameWithOutliers <- outlierFrame %>% mutate(mvOutlier=pcout(outlierFrame, explvar = 0.99, crit.M1 = 0.35, crit.c1=2.5, critM2=0.5, crit.c2= 0.995, outbound = 0.1, cs = 0.45, makeplot = FALSE)$wfinal01)
    noNAsWithOutliers <- left_join(outlierFrameWithOutliers, noNAs)
    output <- left_join(data, noNAsWithOutliers)
    return(output)
}