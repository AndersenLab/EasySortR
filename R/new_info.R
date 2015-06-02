new_info <- function(filepath, levels = 1){
    splitfp <- strsplit(filepath,"/")
    dirname <- splitfp[[1]][(length(splitfp[[1]]) - levels)]

    date <- strsplit(dirname,"_")[[1]][1]

    details <- strsplit(dirname,"_")[[1]][2]

    experiment <- strsplit(details,"[0-9]+")[[1]][1]
    round <- as.numeric(strsplit(details,"(?i)[a-z]+")[[1]][2])
    assay <- strsplit(details,"[0-9]+")[[1]][2]

    split <- strsplit(splitfp[[1]][(length(splitfp[[1]]))],"_")[[1]]
    plate <- as.numeric(strsplit(split[1],"p")[[1]][2])

    straintemplate <- split[2]
    conditiontemplate <- split[3]
    controltemplate <- strsplit(split[4],"\\.")[[1]][1]

    frame <- data.frame(date, experiment, round, assay, plate, straintemplate,
                        conditiontemplate, controltemplate)

    return(frame)
}