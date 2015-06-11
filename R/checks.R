#' Check directory naming
#'
#' @param dirname Path to the directory to be checked
#' @export

checkdir <- function(dirname) {

    # Remove trailing '/' if present in the file path
    
    if (grepl("/$", dirname)){
        dirname <- substr(dirname, 1, nchar(dirname) - 1)
    }
    
    cat("\n")
    cat("\n")
    cat(paste("TESTING EXPERIMENT DIRECTORY:", dirname, "\n"))
    cat("\n")
    topview <- dir(dirname)
    
    # Make sure that all necessary subdirectories are present
    
    cat("Checking that all necessary subdirectories are present...\n")
    cat("\n")
    
    if (!("score" %in% topview)) {
        cat(">> ERROR: No 'score' data subdirectory\n")
    }
    if (!("setup" %in% topview)) {
        cat(">> WARNING: No 'setup' data subdirectory\n")
    }
    if (!("conditions" %in% topview)) {
        cat(">> ERROR: No 'conditions' templates subdirectory\n")
    }
    if (!("contamination" %in% topview)) {
        cat(">> ERROR: No 'contamination' templates subdirectory\n")
    }
    if (!("controls" %in% topview)) {
        cat(">> ERROR: No 'controls' templates subdirectory\n")
    }
    if (!("strains" %in% topview)) {
        cat(">> ERROR: No 'strains' templates subdirectory\n")
    }
    
    # Check that all template names are valid
    
    cat("\n")
    cat("Checking that all template names are valid...\n")
    cat("\n")
    
    conditionsdir <- file.path(dirname, 'conditions')
    contaminationdir <- file.path(dirname, 'contamination')
    controlsdir <- file.path(dirname, 'controls')
    strainsdir <- file.path(dirname, 'strains')
    
    checktemplatenames(conditionsdir)
    checktemplatenames(contaminationdir)
    checktemplatenames(controlsdir)
    checktemplatenames(strainsdir)
    
    # Check that all data file names are correct
    
    cat("\n")
    cat("Checking that all data file names are valid...\n")
    cat("\n")
    
    setupdir <- file.path(dirname, 'setup')
    scoredir <- file.path(dirname, 'score')
    
    cat("Checking setup directory...\n")
    cat("\n")
    
    if (dir.exists(setupdir)) {
        checkdatanames(setupdir, conditionsdir, contaminationdir, controlsdir, strainsdir)
    }
    
    cat("\n")
    cat("Checking score directory...\n")
    cat("\n")
    
    checkdatanames(scoredir, conditionsdir, contaminationdir, controlsdir, strainsdir)
    
    cat("\n")
    cat("######################### DONE #########################\n")
}

checktemplatenames <- function(templatedir) {
    cat(paste("TESTING TEMPLATE SUBDIRECTORY:", templatedir, "\n"))
    cat("\n")
    templates <- dir(templatedir)
    for (template in templates) {
        split <- strsplit(template, "\\.")
        if (length(split[[1]]) != 2) {
            cat(paste0(">> ERROR: Template '",
                         template,
                         "' does not conform to the template naming scheme\n"))
        }
        if (split[[1]][length(split[[1]])] != "csv") {
            cat(paste0(">> ERROR: Template '", template, "' is not a .csv file\n"))
        }
#         if (grepl("_", split[[1]])){
#             cat(paste0("ERROR: Template '", template, "' cannot contain underscores\n"))
#             cat("\n")
#         }
    }
}

checkdatanames <- function(datadir, conditions, contamination, controls, strains) {
    
    dirfiles <- dir(datadir)
    files <- subset(dirfiles, grepl(".txt$", dirfiles))
    
    for (file in files) {
        split1 <- strsplit(file, "\\.")
        split2 <- strsplit(split1[[1]][1], "_")
        if (length(split2[[1]]) < 4) {
            cat(paste0(">> ERROR: File named '", file, "' does not contain a field for all of the folowing: contamination, strains, conditions, controls\n"))
        }
        if (!(paste0(split2[[1]][1], "_contamination.csv") %in% dir(contamination))) {
            cat(paste0(">> ERROR: File named '", file, "' does not match a contamination template file\n"))
        }
        if (!(paste0(split2[[1]][2], ".csv") %in% dir(strains))) {
            cat(paste0(">> ERROR: File named '", file, "' does not match a strains template file\n"))
        }
        if (!(paste0(split2[[1]][3], ".csv") %in% dir(conditions))) {
            cat(paste0(">> ERROR: File named '", file, "' does not match a conditions template file\n"))
        }
        if (!(paste0(split2[[1]][4], ".csv") %in% dir(controls))) {
            cat(paste0(">> ERROR: File named '", file, "' does not match a controls template file\n"))
        }
    }
    
}
