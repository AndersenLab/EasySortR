# Pulls in all information given to a file based on directory structure and
# naming conventions

new_info <- function(filepath, levels = 1){
    
    # Split up the file path and go up the specified number of levels
    
    splitfp <- strsplit(filepath,"/")
    dirname <- splitfp[[1]][(length(splitfp[[1]]) - levels)]
    
    # Check if the direcory name is provides the correct information
    
    if (!grepl("[0-9]{8}_[A-Za-z]+[0-9]*[a-z]{0,1}", dirname)) {
        warning("Directory name does not match predefined structure. Information about the date, name, round , and assay of the experiment has been set to NA. This may be an issue with the number of levels up from the individual file.")
        date <- NA
        experiment <- NA
        round <- NA
        assay <- NA
    } else {
        # Pull out the date
        
        date <- strsplit(dirname,"_")[[1]][1]
        
        # Get the string with all of the experiment details
        
        details <- strsplit(dirname,"_")[[1]][2]
        
        # Get the experiment name, round number, and assay value
        
        experiment <- strsplit(details,"[0-9]+")[[1]][1]
        round <- as.numeric(strsplit(details,"(?i)[a-z]+")[[1]][2])
        assay <- strsplit(details,"[0-9]+")[[1]][2]
    }
    
    # Check to see if the file name is the correct structure
    
    filename <- splitfp[[1]][(length(splitfp[[1]]))]
    
    if (!grepl("p[0-9]{2}_([0-9a-zA-Z\-]*_){2,}[0-9a-zA-Z\-]*\.txt",
               filename)) {
        stop("File name does not conform to the standard 'pXX_straintemplate_conditiontemplate_controltemplate")
    }
    
    # Split the file name and get the plate number
    
    split <- strsplit(filename, "_")[[1]]
    plate <- as.numeric(strsplit(split[1], "p")[[1]][2])

    # Get the template info from the file name
    
    straintemplate <- split[2]
    conditiontemplate <- split[3]
    controltemplate <- strsplit(split[4], "\\.")[[1]][1]
    
    # Return all of the information as a data frame
    
    frame <- data.frame(date, experiment, round, assay, plate, straintemplate,
                        conditiontemplate, controltemplate)

    return(frame)
}