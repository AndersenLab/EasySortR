#' Reads single plate data
#'
#' Reads individual sorter data files into R. Built exlusively for use with worms.
#'
#' @param file The file to be read.
#' @param tofmin The minimum time of flight value allowed. Defaults to 60.
#' @param tofmin The minimum time of flight value allowed. Defaults to 2000.
#' @param extmin The minimum extinction value allowed. Defaults to 0.
#' @param extmin The maximum extinction value allowed. Defaults to 10000.
#' @param SVM Boolean specifying whether or not to use the support vector machine to separate worms and bubbles.
#' @return Returns a single data frame for a single plate file.
#' @import kernlab
#' @import COPASutils
#' @import dplyr
#' @export

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
    strainsfile <- paste0(templatedir, plateinfo$straintemplate[1], ".csv")
    conditionsfile <- paste0(templatedir, plateinfo$conditiontemplate[1],
                             ".csv")
    controlsfile <- paste0(templatedir, plateinfo$controltemplate[1], ".csv")
    contamfile <- paste0(templatedir, sprintf("p%02d", plateinfo$plate[1]),
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