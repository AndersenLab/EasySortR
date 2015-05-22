#' Binned Anomaly, Mitigate and Fit Pruning
#'
#' Prune based on binned IQR multiples eminating from the center of the distribution. Anomalies not removd if they account for more than 5% of the population.
#' 
#' @param data A melted data frame to be analyzed for outliers
#' @return A single data frame for a single plate file with the outlier data NA'd out.
#' @import dplyr
#' @export

bamfPrune <- function(data) {
    data <- ensureLong(data)
    dataWithOutliers <- data %>%
        dplyr::filter(!is.na(strain)) %>%
        dplyr::group_by(condition, trait)%>%
        dplyr::summarise(iqr = IQR(phenotype, na.rm = TRUE),
                         q1 = quantile(phenotype ,probs=.25, na.rm = TRUE),
                         q3 = quantile(phenotype ,probs=.75, na.rm = TRUE))%>%
        dplyr::mutate(cut1h = q3+(iqr*2),
               cut1l =q1-(iqr*2),
               cut2h = q3+(iqr*3),
               cut2l =q1-(iqr*3),
               cut3h = q3+(iqr*4),
               cut3l =q1-(iqr*4),
               cut4h = q3+(iqr*5),
               cut4l =q1-(iqr*5),
               cut5l = q1-(iqr*7),
               cut5h = q3+(iqr*7),
               cut6l = q1-(iqr*10),
               cut6h = q3+(iqr*10))%>%
        dplyr::left_join(data, ., by=c("condition", "trait"))%>%
        dplyr::mutate(onehs = ifelse( cut2h > phenotype & phenotype >= cut1h,1,0),
               onels = ifelse( cut2l < phenotype & phenotype <= cut1l,1,0),
               twohs = ifelse( cut3h > phenotype & phenotype >= cut2h,1,0),
               twols = ifelse( cut3l < phenotype & phenotype <= cut2l,1,0),
               threehs = ifelse(cut4h > phenotype & phenotype >= cut3h,1,0),
               threels = ifelse(cut4l < phenotype & phenotype <= cut3l,1,0),
               fourhs = ifelse(cut5h > phenotype &  phenotype >= cut4h,1,0),
               fourls = ifelse(cut5l < phenotype &  phenotype <= cut4l,1,0),
               fivehs = ifelse(cut6h > phenotype & phenotype >= cut5h,1,0),
               fivels = ifelse(cut6l < phenotype & phenotype <= cut5l,1,0),
               sixhs = ifelse(phenotype >= cut6h,1,0),
               sixls = ifelse(phenotype <= cut6l,1,0))%>%
        dplyr::group_by(condition, trait)%>%
        dplyr::mutate(s1h = sum(onehs, na.rm=TRUE),
               s2h = sum(twohs, na.rm=TRUE),
               s3h = sum(threehs, na.rm=TRUE),
               s4h = sum(fourhs, na.rm=TRUE),
               s5h = sum(fivehs, na.rm=TRUE),
               s1l = sum(onels, na.rm=TRUE),
               s2l = sum(twols, na.rm=TRUE),
               s3l = sum(threels, na.rm=TRUE),
               s4l = sum(fourls, na.rm=TRUE),
               s5l = sum(fivels, na.rm=TRUE),
               s6h = sum(sixhs, na.rm=TRUE),
               s6l = sum(sixls, na.rm=TRUE))%>%
        dplyr::group_by(condition, trait)%>%
        dplyr::mutate(p1h = ifelse(sum(onehs, na.rm=TRUE)/n()>=.05,1,0),
               p2h = ifelse(sum(twohs, na.rm=TRUE)/n()>=.05,1,0),
               p3h = ifelse(sum(threehs, na.rm=TRUE)/n()>=.05,1,0),
               p4h = ifelse(sum(fourhs, na.rm=TRUE)/n()>=.05,1,0),
               p5h = ifelse(sum(fivehs, na.rm=TRUE)/n()>=.05,1,0),
               p6h = ifelse(sum(sixhs, na.rm=TRUE)/n()>=.05,1,0),
               p1l = ifelse(sum(onels, na.rm=TRUE)/n()>=.05,1,0),
               p2l = ifelse(sum(twols, na.rm=TRUE)/n()>=.05,1,0),
               p3l = ifelse(sum(threels, na.rm=TRUE)/n()>=.05,1,0),
               p4l = ifelse(sum(fourls, na.rm=TRUE)/n()>=.05,1,0),
               p5l = ifelse(sum(fivels, na.rm=TRUE)/n()>=.05,1,0),
               p6l = ifelse(sum(sixls, na.rm=TRUE)/n()>=.05,1,0))%>%
        dplyr::mutate(numst = n())%>%
        dplyr::group_by(condition, trait)%>%
        dplyr::filter(!is.na(trait), !is.na(phenotype), !is.na(iqr), !is.na(q1), !is.na(q3), !is.na(cut1h), !is.na(cut1l), !is.na(cut2h), !is.na(cut2l), !is.na(cut3h), !is.na(cut3l), !is.na(cut4h), !is.na(cut4l), !is.na(cut5l), !is.na(cut5h), !is.na(cut6l), !is.na(cut6h), !is.na(onehs), !is.na(onels), !is.na(twohs), !is.na(twols), !is.na(threehs), !is.na(threels), !is.na(fourhs), !is.na(fourls), !is.na(fivehs), !is.na(fivels), !is.na(sixhs), !is.na(sixls), !is.na(s1h), !is.na(s2h), !is.na(s3h), !is.na(s4h), !is.na(s5h), !is.na(s1l), !is.na(s2l), !is.na(s3l), !is.na(s4l), !is.na(s5l), !is.na(s6h), !is.na(s6l), !is.na(p1h), !is.na(p2h), !is.na(p3h), !is.na(p4h), !is.na(p5h), !is.na(p6h), !is.na(p1l), !is.na(p2l), !is.na(p3l), !is.na(p4l), !is.na(p5l), !is.na(p6l), !is.na(numst)) %>%
        dplyr::mutate(cuts = ifelse(sixhs == 1 & ((s6h + s5h + s4h )/numst) <= .05 & (s5h==0|s4h==0), TRUE,
                             ifelse(sixls == 1 & ((s6l + s5l +s4l)/numst) <= .05 & (s5l==0|s4l==0), TRUE, FALSE)),
               cuts1 = ifelse( ((s6h + s5h + s4h )/numst) >= .05 & s6h>=1&s5h>=1&s4h>=1, FALSE,
                               ifelse((s5h>=1&s4h>=1&s3h>=1&s2h>=1&s1h>=1) | (s5h>=1&s3h>=1&s2h>=1&s1h>=1)|(s5h>=1&s4h>=1&s2h>=1&s1h>=1)| (s5h>=1&s4h>=1&s3h>=1&s1h>=1) |(s5h>=1&s4h>=1&s3h>=1&s2h>=1),FALSE,
                                      ifelse((s5h>=1&s4l>=1&s3l>=1&s2l>=1&s1l>=1) |(s5h>=1&s3l>=1&s2l>=1&s1l>=1) | (s5h>=1&s4l>=1&s2l>=1&s1l>=1)|(s5h>=1&s4l>=1&s3l>=1&s1l>=1) | (s5h>=1&s4l>=1&s3l>=1&s2l>=1) ,FALSE,
                                             ifelse(fivehs == 1 & ((s6h+ s4h + s5h + s3h)/numst) <= .05, TRUE,
                                                    ifelse(fivels == 1 & (s6l +s4l + s5l + s3l)/numst <= .05, TRUE, FALSE))))),
               cuts2 = ifelse( ((s6h + s5h + s4h )/numst) >= .05 & s6h>=1&s5h>=1&s4h>=1, FALSE,
                               ifelse((s4h>=1&s3h>=1&s2h>=1&s1h>=1) | (s4h>=1&s2h>=1&s1h>=1)| (s4h>=1&s3h>=1&s1h>=1) |(s4h>=1&s3h>=1&s2h>=1), FALSE,
                                      ifelse((s4l>=1&s3l>=1&s2l>=1&s1l>=1) | (s4l>=1&s2l>=1&s1l>=1)|(s4l>=1&s3l>=1&s1l>=1) | (s4l>=1&s3l>=1&s2l>=1), FALSE,
                                             ifelse(fourhs == 1 & fivehs == 0 & (s5h + s4h + s3h + s2h)/numst <= .05, TRUE,
                                                    ifelse(fourls == 1  & fivels == 0 & (s5l + s4l + s3l + s2l)/numst <= .05, TRUE, FALSE))))))
    
    output <- dataWithOutliers %>% dplyr::select(date, experiment, round, assay, condition, row, col, cuts, cuts1, cuts2) %>% dplyr::rename(bamfOutlier1 = cuts, bamfOutlier2 = cuts1, bamfOutlier3 = cuts2) %>% dplyr::left_join(data, .) %>% dplyr::arrange(data, row, col)
    return(output)
}