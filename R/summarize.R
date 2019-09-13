summarize_plates <- function(plate, quantiles = FALSE, log = FALSE,
                             ends = FALSE, long = FALSE) {
    
    # Select only objects, missing rows, or calls not made
    
    plate <- plate[plate$call50 == "object"
                   | plate$TOF == -1
                   | is.na(plate$call50),]
    
    # Fill missing wells with NA so that all rows are represented
    
    plate <- COPASutils::fillWells(plate)
    
    # Calculate all of the summary statistics for each well
    
    processed <- plate %>%
        dplyr::group_by(date, experiment, round, assay, plate, condition,
                 control, strain, row, col) %>%
        dplyr::summarise(
            
            # Number and number sorted
            
            n = ifelse(length(TOF[!is.na(TOF)]) == 0, NA,
                       length(TOF[!is.na(TOF)])),
            n.sorted=sum(sort == 6),

            # Time of flight values
            
            mean.TOF = mean(TOF, na.rm = TRUE),
            min.TOF = as.numeric(quantile(TOF, na.rm = TRUE)[1]),
            q10.TOF = as.numeric(quantile(TOF, probs = 0.1, na.rm = TRUE)[1]),
            q25.TOF = as.numeric(quantile(TOF, probs = 0.25, na.rm = TRUE)[1]),
            median.TOF = median(TOF, na.rm = TRUE),
            q75.TOF = as.numeric(quantile(TOF, probs = 0.75, na.rm = TRUE)[1]),
            q90.TOF = as.numeric(quantile(TOF, probs = 0.90, na.rm = TRUE)[1]),
            max.TOF = as.numeric(quantile(TOF, na.rm = TRUE)[5]),
            
            # Extinction values

            mean.EXT = mean(EXT, na.rm = TRUE),
            min.EXT = as.numeric(quantile(EXT, na.rm = TRUE)[1]),
            q10.EXT = as.numeric(quantile(EXT, probs = 0.1, na.rm = TRUE)[1]),
            q25.EXT = as.numeric(quantile(EXT, probs = 0.25, na.rm = TRUE)[1]),
            median.EXT = median(EXT, na.rm = TRUE),
            q75.EXT = as.numeric(quantile(EXT, probs = 0.75, na.rm = TRUE)[1]),
            q90.EXT = as.numeric(quantile(EXT, probs = 0.90, na.rm = TRUE)[1]),
            max.EXT = as.numeric(quantile(EXT, na.rm = TRUE)[5]),
            
            # Red values

            mean.red = mean(red, na.rm = TRUE),
            min.red = as.numeric(quantile(red, na.rm = TRUE)[1]),
            q10.red = as.numeric(quantile(red, probs = 0.1, na.rm = TRUE)[1]),
            q25.red = as.numeric(quantile(red, probs = 0.25, na.rm = TRUE)[1]),
            median.red = median(red, na.rm = TRUE),
            q75.red = as.numeric(quantile(red, probs = 0.75, na.rm = TRUE)[1]),
            q90.red = as.numeric(quantile(red, probs = 0.9, na.rm = TRUE)[1]),
            max.red = as.numeric(quantile(red, na.rm = TRUE)[5]),
            
            # Green values

            mean.green = mean(green, na.rm = TRUE),
            min.green = as.numeric(quantile(green, na.rm = TRUE)[1]),
            q10.green = as.numeric(quantile(green, probs = 0.1,
                                            na.rm = TRUE)[1]),
            q25.green = as.numeric(quantile(green, probs = 0.25,
                                            na.rm = TRUE)[1]),
            median.green = median(green, na.rm = TRUE),
            q75.green = as.numeric(quantile(green, probs = 0.75,
                                            na.rm = TRUE)[1]),
            q90.green = as.numeric(quantile(green, probs = 0.9,
                                            na.rm = TRUE)[1]),
            max.green = as.numeric(quantile(green, na.rm = TRUE)[5]),
            
            # Yellow values

            mean.yellow = mean(yellow, na.rm = TRUE),
            min.yellow = as.numeric(quantile(yellow, na.rm = TRUE)[1]),
            q10.yellow = as.numeric(quantile(yellow, probs = 0.1,
                                             na.rm = TRUE)[1]),
            q25.yellow = as.numeric(quantile(yellow, probs = 0.25,
                                             na.rm = TRUE)[1]),
            median.yellow = median(yellow, na.rm = TRUE),
            q75.yellow = as.numeric(quantile(yellow, probs = 0.75,
                                             na.rm = TRUE)[1]),
            q90.yellow = as.numeric(quantile(yellow, probs = 0.9,
                                             na.rm = TRUE)[1]),
            max.yellow = as.numeric(quantile(yellow, na.rm = TRUE)[5]),
            
            # Normalized extinction values

            mean.norm.EXT = mean(norm.EXT, na.rm = TRUE),
            min.norm.EXT = as.numeric(quantile(norm.EXT, na.rm = TRUE)[1]),
            q10.norm.EXT = as.numeric(quantile(norm.EXT, probs = 0.1,
                                             na.rm = TRUE)[1]),
            q25.norm.EXT = as.numeric(quantile(norm.EXT, probs = 0.25,
                                             na.rm = TRUE)[1]),
            median.norm.EXT = median(norm.EXT, na.rm = TRUE),
            q75.norm.EXT = as.numeric(quantile(norm.EXT, probs = 0.75,
                                             na.rm = TRUE)[1]),
            q90.norm.EXT = as.numeric(quantile(norm.EXT, probs = 0.9,
                                             na.rm = TRUE)[1]),
            max.norm.EXT = as.numeric(quantile(norm.EXT, na.rm = TRUE)[5]),
            
            # Normalized red values

            mean.norm.red = mean(norm.red, na.rm = TRUE),
            min.norm.red = as.numeric(quantile(norm.red, na.rm = TRUE)[1]),
            q10.norm.red = as.numeric(quantile(norm.red, probs = 0.1,
                                             na.rm = TRUE)[1]),
            q25.norm.red = as.numeric(quantile(norm.red, probs = 0.25,
                                             na.rm = TRUE)[1]),
            median.norm.red = median(norm.red, na.rm = TRUE),
            q75.norm.red = as.numeric(quantile(norm.red, probs = 0.75,
                                             na.rm = TRUE)[1]),
            q90.norm.red = as.numeric(quantile(norm.red, probs = 0.9,
                                             na.rm = TRUE)[1]),
            max.norm.red = as.numeric(quantile(norm.red, na.rm = TRUE)[5]),
            
            # Normalized green values

            mean.norm.green = mean(norm.green, na.rm = TRUE),
            min.norm.green = as.numeric(quantile(norm.green, na.rm = TRUE)[1]),
            q10.norm.green = as.numeric(quantile(norm.green, probs = 0.1,
                                               na.rm = TRUE)[1]),
            q25.norm.green = as.numeric(quantile(norm.green, probs = 0.25,
                                               na.rm = TRUE)[1]),
            median.norm.green = median(norm.green, na.rm = TRUE),
            q75.norm.green = as.numeric(quantile(norm.green, probs = 0.75,
                                               na.rm = TRUE)[1]),
            q90.norm.green = as.numeric(quantile(norm.green, probs = 0.9,
                                               na.rm = TRUE)[1]),
            max.norm.green = as.numeric(quantile(norm.green, na.rm = TRUE)[5]),
            
            # Normalized yellow values

            mean.norm.yellow = mean(norm.yellow, na.rm = TRUE),
            min.norm.yellow = as.numeric(quantile(norm.yellow,
                                                  na.rm = TRUE)[1]),
            q10.norm.yellow = as.numeric(quantile(norm.yellow, probs = 0.1,
                                                na.rm = TRUE)[1]),
            q25.norm.yellow = as.numeric(quantile(norm.yellow, probs = 0.25,
                                                na.rm = TRUE)[1]),
            median.norm.yellow = median(norm.yellow, na.rm = TRUE),
            q75.norm.yellow = as.numeric(quantile(norm.yellow, probs = 0.75,
                                                na.rm = TRUE)[1]),
            q90.norm.yellow = as.numeric(quantile(norm.yellow, probs = 0.9,
                                                na.rm = TRUE)[1]),
            max.norm.yellow = as.numeric(quantile(norm.yellow,
                                                  na.rm = TRUE)[5]),
            
            # Log-transformed extinction values

            mean.log.EXT = mean(log(EXT), na.rm = TRUE),
            min.log.EXT = as.numeric(quantile(log(EXT), na.rm = TRUE)[1]),
            q10.log.EXT = as.numeric(quantile(log(EXT), probs = 0.1,
                                            na.rm = TRUE)[1]),
            q25.log.EXT = as.numeric(quantile(log(EXT), probs = 0.25,
                                            na.rm = TRUE)[1]),
            median.log.EXT = median(log(EXT), na.rm = TRUE),
            q75.log.EXT = as.numeric(quantile(log(EXT), probs = 0.75,
                                            na.rm = TRUE)[1]),
            q90.log.EXT = as.numeric(quantile(log(EXT), probs = 0.90,
                                            na.rm = TRUE)[1]),
            max.log.EXT = as.numeric(quantile(log(EXT), na.rm = TRUE)[5]),
            
            # Log-transformed red values

            mean.log.red = mean(log(red), na.rm = TRUE),
            min.log.red = as.numeric(quantile(log(red), na.rm = TRUE)[1]),
            q10.log.red = as.numeric(quantile(log(red), probs = 0.1,
                                            na.rm = TRUE)[1]),
            q25.log.red = as.numeric(quantile(log(red), probs = 0.25,
                                            na.rm = TRUE)[1]),
            median.log.red = median(log(red), na.rm = TRUE),
            q75.log.red = as.numeric(quantile(log(red), probs = 0.75,
                                            na.rm = TRUE)[1]),
            q90.log.red = as.numeric(quantile(log(red), probs = 0.90,
                                            na.rm = TRUE)[1]),
            max.log.red = as.numeric(quantile(log(red), na.rm = TRUE)[5]),
            
            # Log-transformed green values

            mean.log.green = mean(log(green), na.rm = TRUE),
            min.log.green = as.numeric(quantile(log(green), na.rm = TRUE)[1]),
            q10.log.green = as.numeric(quantile(log(green), probs = 0.1,
                                              na.rm = TRUE)[1]),
            q25.log.green = as.numeric(quantile(log(green), probs = 0.25,
                                              na.rm = TRUE)[1]),
            median.log.green = median(log(red), na.rm = TRUE),
            q75.log.green = as.numeric(quantile(log(green), probs = 0.75,
                                              na.rm = TRUE)[1]),
            q90.log.green = as.numeric(quantile(log(green), probs = 0.90,
                                              na.rm = TRUE)[1]),
            max.log.green = as.numeric(quantile(log(green), na.rm = TRUE)[5]),
            
            # Log-transformed yellow values

            mean.log.yellow = mean(log(yellow), na.rm = TRUE),
            min.log.yellow = as.numeric(quantile(log(yellow), na.rm = TRUE)[1]),
            q10.log.yellow = as.numeric(quantile(log(yellow), probs = 0.1,
                                               na.rm = TRUE)[1]),
            q25.log.yellow = as.numeric(quantile(log(yellow), probs = 0.25,
                                               na.rm = TRUE)[1]),
            median.log.yellow = median(log(yellow), na.rm = TRUE),
            q75.log.yellow = as.numeric(quantile(log(yellow), probs = 0.75,
                                               na.rm = TRUE)[1]),
            q90.log.yellow = as.numeric(quantile(log(yellow), probs = 0.90,
                                               na.rm = TRUE)[1]),
            max.log.yellow = as.numeric(quantile(log(yellow), na.rm = TRUE)[5]),
            
            # Log-transformed, normalized extinction values

            mean.log.norm.EXT = mean(log(norm.EXT), na.rm = TRUE),
            min.log.norm.EXT = as.numeric(quantile(log(norm.EXT),
                                                   na.rm = TRUE)[1]),
            q10.log.norm.EXT = as.numeric(quantile(log(norm.EXT), probs = 0.1,
                                                 na.rm = TRUE)[1]),
            q25.log.norm.EXT = as.numeric(quantile(log(norm.EXT), probs = 0.25,
                                                 na.rm = TRUE)[1]),
            median.log.norm.EXT = median(log(norm.EXT), na.rm = TRUE),
            q75.log.norm.EXT = as.numeric(quantile(log(norm.EXT), probs = 0.75,
                                                 na.rm = TRUE)[1]),
            q90.log.norm.EXT = as.numeric(quantile(log(norm.EXT), probs = 0.90,
                                                 na.rm = TRUE)[1]),
            max.log.norm.EXT = as.numeric(quantile(log(norm.EXT),
                                                   na.rm = TRUE)[5]),
            
            # Log-transformed, normalized red values

            mean.log.norm.red = mean(log(norm.red), na.rm = TRUE),
            min.log.norm.red = as.numeric(quantile(log(norm.red),
                                                   na.rm = TRUE)[1]),
            q10.log.norm.red = as.numeric(quantile(log(norm.red), probs = 0.1,
                                                 na.rm = TRUE)[1]),
            q25.log.norm.red = as.numeric(quantile(log(norm.red), probs = 0.25,
                                                 na.rm = TRUE)[1]),
            median.log.norm.red = median(log(norm.red), na.rm = TRUE),
            q75.log.norm.red = as.numeric(quantile(log(norm.red), probs = 0.75,
                                                 na.rm = TRUE)[1]),
            q90.log.norm.red = as.numeric(quantile(log(norm.red), probs = 0.90,
                                                 na.rm = TRUE)[1]),
            max.log.norm.red = as.numeric(quantile(log(norm.red),
                                                   na.rm = TRUE)[5]),
            
            # Log-transformed, normalized green values

            mean.log.norm.green = mean(log(norm.green), na.rm = TRUE),
            min.log.norm.green = as.numeric(quantile(log(norm.green),
                                                   na.rm = TRUE)[1]),
            q10.log.norm.green = as.numeric(quantile(log(norm.green),
                                                     probs = 0.1,
                                                   na.rm = TRUE)[1]),
            q25.log.norm.green = as.numeric(quantile(log(norm.green),
                                                     probs = 0.25,
                                                   na.rm = TRUE)[1]),
            median.log.norm.green = median(log(norm.red), na.rm = TRUE),
            q75.log.norm.green = as.numeric(quantile(log(norm.green),
                                                     probs = 0.75,
                                                   na.rm = TRUE)[1]),
            q90.log.norm.green = as.numeric(quantile(log(norm.green),
                                                     probs = 0.90,
                                                   na.rm = TRUE)[1]),
            max.log.norm.green = as.numeric(quantile(log(norm.green),
                                                   na.rm = TRUE)[5]),
            
            # Log-transformed, normalized yellow values

            mean.log.norm.yellow = mean(log(norm.yellow), na.rm = TRUE),
            min.log.norm.yellow = as.numeric(quantile(log(norm.yellow),
                                                    na.rm = TRUE)[1]),
            q10.log.norm.yellow = as.numeric(quantile(log(norm.yellow),
                                                      probs = 0.1,
                                                      na.rm = TRUE)[1]),
            q25.log.norm.yellow = as.numeric(quantile(log(norm.yellow),
                                                    probs = 0.25,
                                                    na.rm = TRUE)[1]),
            median.log.norm.yellow = median(log(norm.yellow), na.rm = TRUE),
            q75.log.norm.yellow = as.numeric(quantile(log(norm.yellow),
                                                    probs = 0.75,
                                                    na.rm = TRUE)[1]),
            q90.log.norm.yellow = as.numeric(quantile(log(norm.yellow),
                                                    probs = 0.90,
                                                    na.rm = TRUE)[1]),
            max.log.norm.yellow = as.numeric(quantile(log(norm.yellow),
                                                    na.rm = TRUE)[5]),
            
            # Variance, CV, and IQR traits

            var.TOF = var(TOF),
            cv.TOF = (sd(TOF,na.rm = TRUE) / mean(TOF,na.rm = TRUE)),
            iqr.TOF = quantile(TOF, na.rm = TRUE, probs = .75) -
                quantile(TOF, na.rm = TRUE, probs = .25),
            var.EXT = var(EXT),
            cv.EXT = (sd(EXT, na.rm = TRUE) / mean(EXT, na.rm = TRUE)),
            iqr.EXT = quantile(EXT, na.rm = TRUE, probs = .75) -
                quantile(EXT, na.rm = TRUE, probs = .25),
            var.red = var(red),
            cv.red = (sd(red, na.rm = TRUE) / mean(red, na.rm = TRUE)),
            iqr.red = quantile(red, na.rm = TRUE, probs = .75) -
                quantile(red, na.rm = TRUE, probs = .25),
            var.green = var(green),
            cv.green = (sd(green, na.rm = TRUE) / mean(green, na.rm = TRUE)),
            iqr.green = quantile(green, na.rm = TRUE, probs = .75) -
                quantile(green, na.rm = TRUE, probs = .25),
            var.yellow = var(yellow),
            cv.yellow = (sd(yellow, na.rm = TRUE) / mean(yellow, na.rm = TRUE)),
            iqr.yellow = quantile(yellow, na.rm = TRUE, probs = .75) -
                quantile(yellow, na.rm = TRUE, probs = .25),
            
            # Life stage fractions

            f.L1 = length(which(stage == "L1")) / length(stage),
            f.L2L3 = length(which(stage == "L2/L3")) / length(stage),
            f.L4 = length(which(stage == "L4")) / length(stage),
            f.ad = length(which(stage == "adult")) / length(stage))

    # Remove the min and max values if not requested
    
    if(!ends){
        processed <- processed[, -(grep("min", colnames(processed)))]
        processed <- processed[, -(grep("max", colnames(processed)))]
    }
    
    # Remove the quantile values if not requested
    
    if(!quantiles){
        processed <- processed[, -(grep("q", colnames(processed)))]
    }
    
    # Remove the log transformed values if not requested
    
    if(!log){
        processed <- processed[, -(grep("log", colnames(processed)))]
    }
    
    # Order the data and NA out missing data from platestitcher.py
    analysis <- processed
    analysis <- analysis[order(analysis$row, analysis$col), ]
    analysis[analysis$mean.TOF == -1 | is.na(analysis$mean.TOF),
             which(colnames(analysis) == "n"):ncol(analysis)] <- NA
    
    # Melt the data with tidyr, if requested by user
    
    if(long){
        analysis <- tidyr::gather(analysis, trait, phenotype,
                                  -c(date, experiment, round, assay, plate,
                                     condition, control, strain, row, col))
    }

    return(analysis)
}


#' Summarize a plate or multiple plates
#'
#' This function takes a list of plates as input and returns a single data frame
#' with data summarized for each well/condition pair. As input, this function
#' takes a two element list, with the first element being a data frame
#' consisting of all of the score plates to be summarized. The second element
#' in the list should be a data frame consisting of all of the setup plates to
#' be summarized, from which the number sorted to each well will be extracted.
#' 
#' Edit: 20190913 KSE - add v3_assay flag. Default = false (no changes to original)
#' If v3_assay = TRUE; expect only a score plate and no sort data.
#'
#' @param plates A single data frame of all data or a list consisting of all of
#' the score plates in the first element and the setup plates in the second
#' elements. Alternatively, a list of data frames can be supplied (the output
#' from \code{read_data} if multiple directories are supplied) if the
#' \code{directories} flag is set to \code{TRUE}.
#' @param directories Set to \code{TRUE} if the supplied data comes from
#' mutliple directories (e.g. the RIAILs1 or RIAILs2 data). Defaults to
#' \code{FALSE}.
#' @param quantiles Boolean indicating whether or not quantile values (q10, q25,
#' q75, q90) should be included in the summarized data frame. Defaults to FALSE.
#' @param log Boolean indicating whether or not log-transformed values should be
#' included in the summarized data frame. Defaults to FALSE.
#' @param ends Boolean indicating whether or not min and max values should be
#' included in the summarized data frame. Defaults to FALSE.
#' @param long Boolean stating whether to output the data in long format.
#' @param v3_assay Boolean stating whether assay was V2 (FALSE) or V3 (TRUE)
#' @importFrom dplyr %>%
#' @export

sumplate <- function(plates, picked_plates = FALSE, directories = FALSE, quantiles = FALSE,
                     log = FALSE, ends = FALSE, long = FALSE, v3_assay = FALSE) {
    
  
  #   if no sort day, generate dummy data (also for v3 assay)
  
  if(picked_plates == TRUE | v3_assay == TRUE){
    sort_plate <- plates 
    sort_plate$sort <- 3
    
    plates <- list(plates, sort_plate)
  }
  
    # If directories iss not set to true and the list looks like it came from
    # multiple directories, alert the user.
    
    if (length(plates) > 2 & !directories) {
        stop("It appears that you may be processing a data set made up of
             multiple directories, but you have not set the `directories` flag
             to `TRUE`. Please try again with `directories = TRUE`.")
    }
    
    if (directories) {
        data <- lapply(plates, function(x) {
            sumplate(x, quantiles = quantiles, log = log, ends = ends,
                     long = long)
        })
        data <- dplyr::rbind_all(data)
        return(data)
    } else {
        # If plates is a list of data frames, summarize n.sorted of setup plate,
        # join it to the score plate, and calculate normalized n
        
        if (length(plates) == 2) {
            score <- summarize_plates(plates[[1]],
                                      quantiles, log, ends, long) %>%
                dplyr::arrange(plate, row, col)
            setup <- plates[[2]] %>%
                dplyr::group_by(date, experiment, round, assay,
                                plate, condition, row, col) %>%
                dplyr::summarize(control.n.sorted = sum(sort == 6))
            
            # if v3_assay = T, norm.n is not calculated
            if(v3_assay == TRUE) {
                data <- dplyr::left_join(score, setup,
                                         by = c("date", "experiment", "round",
                                                "assay", "plate", "condition",
                                                "row", "col")) %>%
                    dplyr::select(-n.sorted, -control.n.sorted)
            } else {
                data <- dplyr::left_join(score, setup,
                                         by = c("date", "experiment", "round",
                                                "assay", "plate", "condition",
                                                "row", "col")) %>%
                    dplyr::mutate(norm.n = n / control.n.sorted) %>%
                    dplyr::select(-n.sorted, -control.n.sorted)
            }
        } else {
            
            # Otherwise, just summarize the plate
            
            data <- summarize_plates(plates, quantiles, log, ends, long) %>%
                dplyr::arrange(plate, row, col)
        }
        return(data)
    }
}


#' Remove contamination from plate data frame or list of plate data frames. Also removes bubble calls from data.
#' Edit KSE 20190913 - for V3 (L1-L4 48h assay), if data is simply a data.frame, remove contaminations without errors
#'
#' @param data Data object from read_data
#' @return The original data object with all contaminated wells removed
#' @export

remove_contamination <- function(data){
  if (class(data[[1]]) == "list") {
    for(i in 1:length(data)){
      data[[i]][[1]] <- dplyr::filter(data[[i]][[1]], !contamination)%>%
        dplyr::filter(call50 != "bubble")
      data[[i]][[2]] <- dplyr::filter(data[[i]][[2]], !contamination)%>%
        dplyr::filter(call50 != "bubble")
    }
  }
  else if(class(data) == "data.frame") {
      data <- dplyr::filter(data, !contamination)%>%
          dplyr::filter(call50 != "bubble")
    }
    else {
        data[[1]] <- dplyr::filter(data[[1]], !contamination)%>%
          dplyr::filter(call50 != "bubble")
        data[[2]] <- dplyr::filter(data[[2]], !contamination)%>%
          dplyr::filter(call50 != "bubble")
      }
  return(data)
}
