#' Read plate data
#'
#' This is the primary fuction for reading data into R with this package. It is
#' built exlusively for use with worm data. This function can read in
#' individual plates, plates from a single directory (e.g. setup or score
#' directory), or plates from an experimental directory (with both score and
#' setup sub-directories). The files/directories to be read in (individual
#' file/directory/experiment directory) are detected and handled automatically
#' by the function. This function also utilizes an SVM to determine whether
#' objects are more likely to be bubbles or worms and labels them accordingly.
#' It also estimates the lifestage of each worm.
#' 
#' @param filedir The file or directory to be read. This can be an individual
#' file, a directory containing .txt files to be read, or the directory
#' containing the score and setup directories.
#' @param tofmin The minimum time of flight value allowed. Defaults to 60.
#' @param tofmax The maximum time of flight value allowed. Defaults to 2000.
#' @param extmin The minimum extinction value allowed. Defaults to 0.
#' @param extmax The maximum extinction value allowed. Defaults to 10000.
#' @param SVM Boolean specifying whether or not to use the support vector
#' machine to separate worms and bubbles.
#' @param levels The number of levels above the individual files for the
#' directory containing information about experiment date, name, round, and
#' assay. Defaults to 2, which is the standard for Andersen Lab file structure
#' as of June 2015.
#' @return If a single file is given, a single data frame for only the provided
#' plate is returned. If an experiment directory is given, a list of two data
#' frames will be returned. The first element in the list will be a single data
#' frame for all of the score plates. The second element will be all of the
#' setup plates. If there are no setup plates, NA will fill the second list
#' element. If a setup or score directory is given, a single data frame of all
#' of the plates in the directory will be returned.
#' @export

read_data <- function(filedir, tofmin = 60, tofmax = 2000, extmin = 0,
                      extmax = 10000, SVM = TRUE, levels = 2) {

    # NOTE: 'The next two lines were added to get known issue with
    # dplyr::left_join treating character NAs incorrectly. This is currently
    # (6/9/2015) issue number 965 for the dplyr package. It should be corrected
    # in dplyr version 0.5. At that time, remove the two lines below as well as
    # the reseting of the option at the end of the function and test the
    # function to ensure that the strain names are still present after the point
    # of the final in bamf_prune.
    
    saf <- getOption("stringsAsFactors")
    options(stringsAsFactors = TRUE)
    
    
    # Remove trailing '/' if present in the file path
    
    if (grepl("/$", filedir)){
        filedir <- substr(filedir, 1, nchar(filedir) - 1)
    }
    
    # If an individual file is given, read it in
    
    if (length(dir(filedir)) == 0){
        data <- read_file(filedir, tofmin, tofmax, extmin, extmax, SVM, levels)
    } else if ("score" %in% dir(filedir) & "setup" %in% dir(filedir)) {
        
        # If both setup and score directories are subdirectories of the given
        # directory, read in both directories and return them as a list
        
        scorepath <- file.path(filedir, "score")
        setuppath <- file.path(filedir, "setup")
        
        # If the 'no files in directory' error is thrown by read_directory, 
        # make it a bit more specific when it is thrown from read_data
        
        score <- read_directory(scorepath, tofmin, tofmax, extmin, extmax, SVM,
                                levels)
        setup <- read_directory(setuppath, tofmin, tofmax, extmin, extmax, SVM,
                                levels)
        data <- list(score, setup)
    } else if ("score" %in% dir(filedir) & !("setup" %in% dir(filedir))) {
        
        # If in an experiment directory but there is no setup directory, read in
        # everything from the score directory
        
        scorepath <- file.path(filedir, "score")
        score <- read_directory(scorepath, tofmin, tofmax, extmin, extmax, SVM,
                                levels)
        data <- score
    } else {
        
        # Otherwise, just read in the given directory
        
        data <- read_directory(filedir, tofmin, tofmax, extmin, extmax, SVM,
                               levels)
    }
    
    # NOTE: Remove the following line after dplyr updates to > v0.5. 
    
    options(stringsAsFactors = saf)
    
    return(data)
}

read_directory <- function(directory, tofmin = 60, tofmax = 2000, extmin = 0,
                           extmax = 10000, SVM = TRUE, levels = 2) {
    
    # Get all of the txt files from adirectory and read them in individually.
    # Then rbind them and return a data frame.
    
    dirfiles <- dir(directory)
    files <- subset(dirfiles, grepl(".txt$", dirfiles))
    
    # Throw an error if there are no files in the given directory
    
    if (length(files) == 0) stop("There are no text files in the given directory.")
    
    filepaths <- file.path(directory, files)
    plates <- lapply(filepaths, function(x){
        read_file(x, tofmin, tofmax, extmin, extmax,
                  SVM, levels)
    })
    data <- do.call(rbind, plates)
    return(data)
}

read_file <- function(file, tofmin = 60, tofmax = 2000, extmin = 0,
                      extmax = 10000, SVM = TRUE, levels = 2){
    
    # Read the raw sorter files and make the row names
    
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
    
    # Extract the time so that it is realtive to the first worm sorted
    
    modplate <- modplate %>%
        dplyr::group_by(row, col) %>%
        dplyr::do(COPASutils::extractTime(.))
    modplate <- data.frame(modplate)
    
    # Normalize the optical values by time of flight
    
    modplate[, 10:13] <- apply(modplate[, c(5, 7:9)], 2,
                               function(x) x / modplate$TOF)
    colnames(modplate)[10:13] <- c("norm.EXT", "norm.green", "norm.yellow",
                                   "norm.red")
    
    # Handle the SVM predictions if requested
    
    if(SVM){
        plateprediction <- kernlab::predict(
            COPASutils::bubbleSVMmodel_noProfiler,
            modplate[,3:length(modplate)],
            type="probabilities")
        modplate$object <- plateprediction[, "1"]
        modplate$call50 <- factor(as.numeric(modplate$object > 0.5),
                                  levels=c(1, 0), labels=c("object", "bubble"))
    }
    
    # Calculate the life stage values based on the size of the worms
    
    modplate$stage <- ifelse(modplate$TOF >= 60 & modplate$TOF < 90, "L1",
                             ifelse(modplate$TOF >= 90 & modplate$TOF < 200,
                                    "L2/L3",
                                    ifelse(modplate$TOF >= 200
                                           & modplate$TOF < 300, "L4",
                                           ifelse(modplate$TOF >= 300,
                                                  "adult", NA))))
    
    # Convert integer values to numerics
    
    modplate[, as.vector(which(lapply(modplate, class) == "integer"))] <- lapply(
        modplate[, as.vector(which(lapply(modplate, class) == "integer"))],
        as.numeric)
    
    # Get info about the plate using the new_info function 
    
    plateinfo <- new_info(file, levels)
    
    # Get the template base directory
    
    templatedir <- strsplit(file, "/")[[1]]
    templatedir <- templatedir[-c(length(templatedir), length(templatedir) - 1)]
    templatedir <- paste0(templatedir, collapse = "/")
    templatedir <- paste0(templatedir, "/")
    
    # Get the template file paths
    
    strainsfile <- paste0(templatedir, "strains/",
                          plateinfo$straintemplate[1], ".csv")
    conditionsfile <- paste0(templatedir, "conditions/",
                             plateinfo$conditiontemplate[1],
                             ".csv")
    controlsfile <- paste0(templatedir, "controls/",
                           plateinfo$controltemplate[1], ".csv")
    contamfile <- paste0(templatedir,
                         "contamination/",
                         sprintf("p%02d", plateinfo$plate[1]),
                         "_contamination.csv")
    
    # Read all of the templates
    
    strains <- read_template(strainsfile, type="strains")
    conditions  <- read_template(conditionsfile, type="conditions")
    controls <- read_template(controlsfile, type="controls")
    contam <- read_template(contamfile, type="contam")
    
    # Join all of the metadata and template info to the plate data
    
    modplate <- cbind(plateinfo[,1:5], modplate)
    modplate <- dplyr::left_join(modplate, strains, by = c("row", "col"))
    modplate <- dplyr::left_join(modplate, conditions, by = c("row", "col"))
    modplate <- dplyr::left_join(modplate, controls, by = c("row", "col"))
    modplate <- dplyr::left_join(modplate, contam, by = c("row", "col"))
    
    return(modplate)
}

read_template <- function(templatefile, type){
    
    if (!file.exists(templatefile)) {
        stop(paste("The", type, "template file at", templatefile, "could not be
                   found. Please check to ensure that file exists and is in the
                   right directory structure, then try again."))
    }
    
    # Read in the teamplate file and melt it with tidyr
    
    template <- read.csv(templatefile, check.names=FALSE)
    melttemplate <- tidyr::gather(template, col, variable, -row)
    
    # Change the column names based on the type of template and return the data
    
    if(type == "strains"){
        colnames(melttemplate) <- c("row", "col", "strain")
    } else if (type == "conditions"){
        colnames(melttemplate) <- c("row", "col", "condition")
    } else if (type == "controls"){
        colnames(melttemplate) <- c("row", "col", "control")
    } else if (type == "contam") {
        colnames(melttemplate) <- c("row", "col", "contamination")
        melttemplate$contamination[is.na(melttemplate$contamination)] <- FALSE
    }
    
    return(melttemplate)
}