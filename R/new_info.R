# Pulls in all information given to a file based on directory structure and
# naming conventions

new_info <- function(filepath, levels = 1){
    #Split up the file path and go up the specified number of levels
    
    splitfp <- strsplit(filepath,"/")
    dirname <- splitfp[[1]][(length(splitfp[[1]]) - levels)]

    # Pull out the date
    
    date <- strsplit(dirname,"_")[[1]][1]
    
    # Get the string with all of the experiment details
    
    details <- strsplit(dirname,"_")[[1]][2]
    
    # Get the experiment name, round number, and assay value
    
    experiment <- strsplit(details,"[0-9]+")[[1]][1]
    round <- as.numeric(strsplit(details,"(?i)[a-z]+")[[1]][2])
    assay <- strsplit(details,"[0-9]+")[[1]][2]
    
    # Split the file name and get the plate number
    
    split <- strsplit(splitfp[[1]][(length(splitfp[[1]]))],"_")[[1]]
    plate <- as.numeric(strsplit(split[1],"p")[[1]][2])

    # Get the template info from the file name
    
    straintemplate <- split[2]
    conditiontemplate <- split[3]
    controltemplate <- strsplit(split[4],"\\.")[[1]][1]
    
    # Return all of the information as a data frame
    
    frame <- data.frame(date, experiment, round, assay, plate, straintemplate,
                        conditiontemplate, controltemplate)

    return(frame)
}