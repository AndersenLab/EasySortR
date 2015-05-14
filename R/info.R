#' info
#'
#'
#' @param filePath The path of the file
#' @param levels The number of levels above the file name that identifying information is found. For "~/Dropbox/HTA/Results/20140317_GWAS1a/score/p01_control.txt" this would be 2, since 'score' is one level up and '20140317_GWAS1a' is 2 levels up.
#' @return A one row data frame with a column for the date, experiment, round, assay, plate, and drug. 
#' @export

info <- function(filePath, levels = 1){
    splitfp <- strsplit(filePath,"/")
    dirName <- splitfp[[1]][(length(splitfp[[1]])-levels)]
    
    date <- strsplit(dirName,"_")[[1]][1]
    
    details <- strsplit(dirName,"_")[[1]][2]
    
    experiment <- strsplit(details,"[0-9]+")[[1]][1]
    round <- as.numeric(strsplit(details,"(?i)[a-z]+")[[1]][2])
    assay <- strsplit(details,"[0-9]+")[[1]][2]
    
    split <- strsplit(splitfp[[1]][(length(splitfp[[1]]))],"_")[[1]]
    condition <- strsplit(split[2],"\\.")[[1]][1]
    plate <- as.numeric(strsplit(split[1],"p")[[1]][2])
    
    frame <- data.frame(date,experiment,round,assay,plate,condition)
    
    return(frame)
}