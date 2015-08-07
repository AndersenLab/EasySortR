#' Binned Anomaly, Mitigate and Fit Pruning
#'
#' Prune based on binned IQR multiples emanating from the center of the
#' distribution. Anomalies are not removed if they account for more than 5\% of
#' the population.
#'
#' @param data A melted data frame to be analyzed for outliers.
#' @param drop A boolean stating whether observations from outlier data points
#' should be dropped from the returned data set. Defaults to \code{FALSE}.
#' @return A data frame either with (\code{drop = TRUE}) all of the outlier data
#' points removed from the data frame or (\code{drop = FALSE}) three additional
#' columns stating whether the data point was identified as an outlier in any of
#' the three conditions.
#' @importFrom dplyr %>%
#' @export

bamf_prune <- function(data, drop = FALSE) {
    
    # Make sure that the data being fed into the pruning function is in long
    # format
    
    data <- ensure_long(data)
    
    napheno <- data[is.na(data[, "phenotype"]), ] %>%
        dplyr::mutate(bamfoutlier1 = NA, bamfoutlier2 = NA, bamfoutlier3 = NA)
    
    datacuts <- data %>%
        
        # Filter out all of the wash and/or empty wells
        
        dplyr::filter(!is.na(strain)) %>%
        
        # Group on condition and trait, the, calculate the first and third
        # quartiles for each of the traits
        
        dplyr::group_by(condition, trait) %>%
        dplyr::summarise(iqr = IQR(phenotype, na.rm = TRUE),
                         q1 = quantile(phenotype, probs = .25, na.rm = TRUE),
                         q3 = quantile(phenotype, probs = .75,
                                       na.rm = TRUE)) %>%
        
        # Add a column for the boundaries of each of the bins
        
        dplyr::mutate(cut1h = q3 + (iqr * 2),
                      cut1l =q1 - (iqr * 2),
                      cut2h = q3 + (iqr * 3),
                      cut2l =q1 - (iqr * 3),
                      cut3h = q3 + (iqr * 4),
                      cut3l =q1 - (iqr * 4),
                      cut4h = q3 + (iqr * 5),
                      cut4l =q1 - (iqr * 5),
                      cut5l = q1 - (iqr * 7),
                      cut5h = q3 + (iqr * 7),
                      cut6l = q1 - (iqr * 10),
                      cut6h = q3 + (iqr * 10)) %>%
        
        # Join the bin boundaries back to the original data frame
        
        dplyr::left_join(data, ., by=c("condition", "trait"))
    
    # Add columns tallying the total number of points in each of the bins
    
    datawithoutliers <- datacuts %>%
        dplyr::mutate(onehs = as.numeric(cut2h > phenotype &
                                             phenotype >= cut1h)) %>%
        dplyr::mutate(onels = as.numeric(cut2l < phenotype &
                                             phenotype <= cut1l)) %>%
        dplyr::mutate(twohs = as.numeric(cut3h > phenotype &
                                             phenotype >= cut2h)) %>%
        dplyr::mutate(twols = as.numeric(cut3l < phenotype &
                                             phenotype <= cut2l)) %>%
        dplyr::mutate(threehs = as.numeric(cut4h > phenotype &
                                               phenotype >= cut3h)) %>%
        dplyr::mutate(threels = as.numeric(cut4l < phenotype &
                                               phenotype <= cut3l)) %>%
        dplyr::mutate(fourhs = as.numeric(cut5h > phenotype &
                                              phenotype >= cut4h)) %>%
        dplyr::mutate(fourls = as.numeric(cut5l < phenotype &
                                              phenotype <= cut4l)) %>%
        dplyr::mutate(fivehs = as.numeric(cut6h > phenotype &
                                              phenotype >= cut5h)) %>%
        dplyr::mutate(fivels = as.numeric(cut6l < phenotype &
                                              phenotype <= cut5l)) %>%
        dplyr::mutate(sixhs = as.numeric(phenotype >= cut6h)) %>%
        dplyr::mutate(sixls = as.numeric(phenotype <= cut6l))
    
    # Group on condition and trait, then sum the total number of data points
    # in each of the IQR multiple bins
    
    datawithoutliers <- datawithoutliers %>%
        dplyr::group_by(condition, trait) %>%
        dplyr::mutate(s1h = sum(onehs, na.rm = TRUE),
                      s2h = sum(twohs, na.rm = TRUE),
                      s3h = sum(threehs, na.rm = TRUE),
                      s4h = sum(fourhs, na.rm = TRUE),
                      s5h = sum(fivehs, na.rm = TRUE),
                      s1l = sum(onels, na.rm = TRUE),
                      s2l = sum(twols, na.rm = TRUE),
                      s3l = sum(threels, na.rm = TRUE),
                      s4l = sum(fourls, na.rm = TRUE),
                      s5l = sum(fivels, na.rm = TRUE),
                      s6h = sum(sixhs, na.rm = TRUE),
                      s6l = sum(sixls, na.rm = TRUE))
        
    # Group on condition and trait, then check to see if the number of
    # points in each bin is more than 5% of the total number of data points
    
    datawithoutliers <- datawithoutliers %>%
        dplyr::group_by(condition, trait) %>%
        dplyr::mutate(p1h = ifelse(sum(onehs, na.rm = TRUE) / n() >= .05,1,0),
                      p2h = ifelse(sum(twohs, na.rm = TRUE) / n() >= .05,1,0),
                      p3h = ifelse(sum(threehs, na.rm = TRUE) / n() >= .05,1,0),
                      p4h = ifelse(sum(fourhs, na.rm = TRUE) / n() >= .05,1,0),
                      p5h = ifelse(sum(fivehs, na.rm = TRUE) / n() >= .05,1,0),
                      p6h = ifelse(sum(sixhs, na.rm = TRUE) / n() >= .05,1,0),
                      p1l = ifelse(sum(onels, na.rm = TRUE) / n() >= .05,1,0),
                      p2l = ifelse(sum(twols, na.rm = TRUE) / n() >= .05,1,0),
                      p3l = ifelse(sum(threels, na.rm = TRUE) / n() >= .05,1,0),
                      p4l = ifelse(sum(fourls, na.rm = TRUE) / n() >= .05,1,0),
                      p5l = ifelse(sum(fivels, na.rm = TRUE) / n() >= .05,1,0),
                      p6l = ifelse(sum(sixls,
                                       na.rm = TRUE) / n() >= .05,1,0))
        
    # Count the number of observations in each condition/trait combination
    
    datawithoutliers <- datawithoutliers %>%
        dplyr::mutate(numst = n())
        
    # Group on condition and trait, then filter out NAs in any of the added
    # columns
        
    datawithoutliers <- datawithoutliers %>%
        dplyr::group_by(condition, trait) %>%
        dplyr::filter(!is.na(trait), !is.na(phenotype), !is.na(iqr), !is.na(q1),
                      !is.na(q3), !is.na(cut1h), !is.na(cut1l), !is.na(cut2h),
                      !is.na(cut2l), !is.na(cut3h), !is.na(cut3l),
                      !is.na(cut4h), !is.na(cut4l), !is.na(cut5l),
                      !is.na(cut5h), !is.na(cut6l), !is.na(cut6h),
                      !is.na(onehs), !is.na(onels), !is.na(twohs),
                      !is.na(twols), !is.na(threehs), !is.na(threels),
                      !is.na(fourhs), !is.na(fourls), !is.na(fivehs),
                      !is.na(fivels), !is.na(sixhs), !is.na(sixls),
                      !is.na(s1h), !is.na(s2h), !is.na(s3h), !is.na(s4h),
                      !is.na(s5h), !is.na(s1l), !is.na(s2l), !is.na(s3l),
                      !is.na(s4l), !is.na(s5l), !is.na(s6h), !is.na(s6l),
                      !is.na(p1h), !is.na(p2h), !is.na(p3h), !is.na(p4h),
                      !is.na(p5h), !is.na(p6h), !is.na(p1l), !is.na(p2l),
                      !is.na(p3l), !is.na(p4l), !is.na(p5l), !is.na(p6l),
                      !is.na(numst))
        
    # Add three columns stating whether the observation is an outlier
    # based the three outlier detection functions below
    
    datawithoutliers <- datawithoutliers %>%
        dplyr::ungroup()
    
    datawithoutliers <- datawithoutliers %>%
        dplyr::mutate(cuts = categorize1(.),
                      cuts1 = categorize2(.),
                      cuts2 = categorize3(.))
    
    # Necessary to prevent segfaults
    
    datawithoutliers <- data.frame(datawithoutliers)
    
    # Select the necessary columns, rename the outlier calling columns, and
    # arrange by condition, row, and columns. Generally clean things up and bind
    # it back to the original data frame to regain lost data.
    
    primaryselect <- datawithoutliers %>%
        dplyr::select(date, experiment, round, assay, condition, plate, row,
                      col, trait, phenotype, cuts, cuts1, cuts2) 
    renamedcols <- primaryselect %>%
        dplyr::rename(bamfoutlier1 = cuts, bamfoutlier2 = cuts1,
                      bamfoutlier3 = cuts2) 
    joineddata <- dplyr::left_join(renamedcols, data,
                                   by = c("date", "experiment", "round",
                                          "assay", "condition", "plate", "row",
                                          "col", "trait", "phenotype")) 
    arrangedcols <- joineddata %>%
        dplyr::arrange(plate, row, as.numeric(col), trait)
    output <- arrangedcols %>%
        dplyr::select(date, experiment, round, assay, condition, control, plate,
                      row, col, strain, trait, phenotype, bamfoutlier1,
                      bamfoutlier2, bamfoutlier3)
    
    # Necessary to prevent segfaults
    
    output <- data.frame(output)
    
    # If the user wants to drop the outliers, filter against all three columns
    
    if (drop) {
        output <- output %>%
            dplyr::filter(!bamfoutlier1 & !bamfoutlier2 & !bamfoutlier3) %>%
            dplyr::select(-bamfoutlier1, -bamfoutlier2, -bamfoutlier3)
    }
    
    # Return the output data frame
    
    return(output)
}

# If the observation is in the sixth bin (>10x IQR outside the distribution)
# and the three outermost bins are discontinuous and make up less than 5% of the
# distribution, mark the observation as an outlier.

categorize1 <- function(data) {
    with(data,
         (sixhs >= 1 & ( (s6h + s5h + s4h ) / numst) <= .05
          & (s5h == 0 | s4h == 0))
         | (sixls >= 1 & ( (s6l + s5l + s4l) / numst) <= .05
            & (s5l == 0 | s4l == 0))
    )
}

# If the 5 innermost bins are discontinuous by more than a 1 bin gap, the
# observation is in the fifth bin (between 7 and 10x IQR outside the
# distribution), and the four outermost bins make up less than 5% of the
# population, mark the observation an outlier

categorize2 <- function(data) {
    with(data,
         ( ( ! ( (s5h >= 1 & s4h >= 1 & s3h >= 1 & s2h >= 1 & s1h >= 1)
                 | (s5h >= 1 & s3h >= 1 & s2h >= 1 & s1h >= 1)
                 | (s5h >= 1 & s4h >= 1 & s2h >= 1 & s1h >= 1)
                 | (s5h >= 1 & s4h >= 1 & s3h >= 1 & s1h >= 1)
                 | (s5h >= 1 & s4h >= 1 & s3h >= 1 & s2h >= 1)))
           & (fivehs == 1 & ( (s6h + s5h + s4h + s3h) / numst) <= .05))
         | ( ( ! ( (s5h >= 1 & s4l >= 1 & s3l >= 1 & s2l >= 1 & s1l >= 1)
                   | (s5h >= 1 & s3l >= 1 & s2l >= 1 & s1l >= 1)
                   | (s5h >= 1 & s4l >= 1 & s2l >= 1 & s1l >= 1)
                   | (s5h >= 1 & s4l >= 1 & s3l >= 1 & s1l >= 1)
                   | (s5h >= 1 & s4l >= 1 & s3l >= 1 & s2l >= 1)))
             & (fivels == 1 & ( (s6l + s5l + s4l + s3l) / numst) <= .05))
    )
}

# If the 4 innermost bins are discontinuous by more than a 1 bin gap, the
# observation is in the fourth bin (between 5 and 7x IQR outside the
# distribution), and the four outermost bins make up less than 5% of the
# population, mark the observation an outlier

categorize3 <- function(data) {
    with(data,
         ( ( ! ( (s4h >= 1 & s3h >= 1 & s2h >= 1 & s1h >= 1)
                 | (s4h >= 1 & s2h >= 1 & s1h >= 1)
                 | (s4h >= 1 & s3h >= 1 & s1h >= 1)
                 | (s4h >= 1 & s3h >= 1 & s2h >= 1)))
           & (fourhs == 1 & (s5h + s4h + s3h + s2h) / numst <= .05))
         | ( ( ! ( (s4l >= 1 & s3l >= 1 & s2l >= 1 & s1l >= 1)
                   | (s4l >= 1 & s2l >= 1 & s1l >= 1)
                   | (s4l >= 1 & s3l >= 1 & s1l >= 1)
                   | (s4l >= 1 & s3l >= 1 & s2l >= 1)))
             & (fourls == 1 & (s5l + s4l + s3l + s2l) / numst <= .05))
    )
}

#' Prune data based on biological cutoffs
#'
#' Prune the data based on brood size and/or normalized brood size
#'
#' @param data A melted data frame to be analyzed for outliers.
#' @param drop A boolean stating whether observations from outlier data points
#' should be dropped from the returned data set. Defaults to \code{FALSE}.
#' @return A data frame either with (\code{drop = TRUE}) all of the outlier data
#' points removed from the data frame or (\code{drop = FALSE}) three additional
#' columns stating whether the data point was identified as an outlier in any of
#' the three conditions.
#' @importFrom dplyr %>%
#' @export

bioprune <- function(data){
    if ("norm.n" %in% colnames(data)){
        biopruneddata <- data %>%
            dplyr::filter(n > 5, n < 1000, norm.n < 350)
    } else {
        biopruneddata <- data %>%
            dplyr::filter(n > 5, n < 1000)
    }
    return(biopruneddata)
}