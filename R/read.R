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
#' @return If a single file is given, a single data frame for only the provided
#' plate is returned. If an experiment directory is given, a list of two data
#' frames will be returned. The first element in the list will be a single data
#' frame for all of the score plates. The second element will be all of the
#' setup plates. If there are no setup plates, NA will fill the second list
#' element. If a setup or score directory is given, a single data frame of all
#' of the plates in the directory will be returned.
#' @export

read_data <- function(filedir, tofmin=60, tofmax=2000, extmin=0, extmax=10000,
    SVM=TRUE) {
    if (grep("/$", filedir) == 1){
        filedir <- substr(filedir, 1, nchar(filedir) - 1)
    }
    if (length(dir(filedir)) == 0){
        data <- read_file(filedir, tofmin, tofmax, extmin, extmax, SVM)
    } else if ("score" %in% dir(filedir) & "setup" %in% dir(filedir)) {
        scorepath <- file.path(filedir, "score")
        setuppath <- file.path(filedir, "setup")
        score <- read_directory(scorepath, tofmin, tofmax, extmin, extmax, SVM)
        setup <- read_directory(setuppath, tofmin, tofmax, extmin, extmax, SVM)
        data <- list(score, setup)
    } else if ("score" %in% dir(filedir) & !("setup" %in% dir(filedir))) {
        scorepath <- file.path(filedir, "score")
        score <- read_directory(scorepath, tofmin, tofmax, extmin, extmax, SVM)
        data <- score
    } else {
        data <- read_directory(filedir, tofmin, tofmax, extmin, extmax, SVM)
    }
    return(data)
}

read_directory <- function(directory, tofmin=60, tofmax=2000, extmin=0,
                           extmax=10000, SVM=TRUE) {
    dirfiles <- dir(directory)
    files <- subset(dirfiles, grepl(".txt$", dirfiles))
    filepaths <- file.path(directory, files)
    plates <- lapply(filepaths, function(x){
        read_file(x, tofmin, tofmax, extmin, extmax,
                  SVM)
    })
    data <- do.call(rbind, plates)
    return(data)
}

read_file <- function(file, tofmin=60, tofmax=2000, extmin=0, extmax=10000,
                      SVM=TRUE, levels=2){
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
    modplate[, 10:13] <- apply(modplate[, c(5, 7:9)], 2,
                               function(x) x / modplate$TOF)
    colnames(modplate)[10:13] <- c("norm.EXT", "norm.green", "norm.yellow",
                                   "norm.red")
    if(SVM){
        plateprediction <- kernlab::predict(
            COPASutils::bubbleSVMmodel_noProfiler,
            modplate[,3:length(modplate)],
            type="probabilities")
        modplate$object <- plateprediction[, "1"]
        modplate$call50 <- factor(as.numeric(modplate$object > 0.5),
                                  levels=c(1, 0), labels=c("object", "bubble"))
    }
    modplate$stage <- ifelse(modplate$TOF >= 60 & modplate$TOF < 90, "L1",
                             ifelse(modplate$TOF >= 90 & modplate$TOF < 200,
                                    "L2/L3",
                                    ifelse(modplate$TOF >= 200
                                           & modplate$TOF < 300, "L4",
                                           ifelse(modplate$TOF >= 300,
                                                  "adult", NA))))
    modplate[, as.vector(which(lapply(modplate, class) == "integer"))] <- lapply(
        modplate[, as.vector(which(lapply(modplate, class) == "integer"))],
        as.numeric)
    
    plateinfo <- new_info(file, levels)
    
    templatedir <- strsplit(file, "/")[[1]]
    templatedir <- templatedir[-c(length(templatedir), length(templatedir) - 1)]
    templatedir <- paste0(templatedir, collapse = "/")
    templatedir <- paste0(templatedir, "/")
    
    #### Put in new file path to templates once known
    strainsfile <- paste0(templatedir, "strains",
                          plateinfo$straintemplate[1], ".csv")
    conditionsfile <- paste0(templatedir, "conditions",
                             plateinfo$conditiontemplate[1],
                             ".csv")
    controlsfile <- paste0(templatedir, "controls",
                           plateinfo$controltemplate[1], ".csv")
    contamfile <- paste0(templatedir,
                         "contamination",
                         sprintf("p%02d", plateinfo$plate[1]),
                         "_contamination.csv")
    
    strains <- read_template(strainsfile, type="strains")
    conditions  <- read_template(conditionsfile, type="conditions")
    controls <- read_template(controlsfile, type="controls")
    contam <- read_template(contamfile, type="contam")
    
    modplate <- cbind(plateinfo[,1:5], modplate)
    modplate <- dplyr::left_join(modplate, strains, by = c("row", "col"))
    modplate <- dplyr::left_join(modplate, conditions, by = c("row", "col"))
    modplate <- dplyr::left_join(modplate, controls, by = c("row", "col"))
    modplate <- dplyr::left_join(modplate, contam, by = c("row", "col"))
    
    return(modplate)
}

read_template <- function(templatefile, type){
    template <- read.csv(templatefile, check.names=FALSE)
    melttemplate <- tidyr::gather(template, col, variable, -row)
    if(type == "strains"){
        colnames(melttemplate) <- c("row", "col", "strain")
    } else if (type == "conditions"){
        colnames(melttemplate) <- c("row", "col", "condition")
    } else if (type == "controls"){
        colnames(melttemplate) <- c("row", "col", "control")
    } else if (type == "contam") {
        colnames(melttemplate) <- c("row", "col", "contamination")
    }
    return(melttemplate)
}