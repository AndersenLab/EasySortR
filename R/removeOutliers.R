# Takes long format data frame

removeOutliers <- function(df){
  cutPheno <- df %>%
    group_by(ids)%>%
    summarise(iqr = IQR(value, na.rm = T),
              q1 = quantile(value,probs=.25),
              q3 = quantile(value,probs=.75))%>%
    group_by(ids) %>%
    mutate(cut1h = q3+(iqr*2), 
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
    left_join(df,.,by=c("ids"))%>%
    mutate(onehs = ifelse( cut2h > value & value >= cut1h,1,0),
           onels = ifelse( cut2l < value & value <= cut1l,1,0),
           twohs = ifelse( cut3h > value & value >= cut2h,1,0),
           twols = ifelse( cut3l < value & value <= cut2l,1,0),
           threehs = ifelse(cut4h > value & value >= cut3h,1,0),
           threels = ifelse(cut4l < value & value <= cut3l,1,0),
           fourhs = ifelse(cut5h > value &  value >= cut4h,1,0),
           fourls = ifelse(cut5l < value &  value <= cut4l,1,0),
           fivehs = ifelse(cut6h > value & value >= cut5h,1,0),
           fivels = ifelse(cut6l < value & value <= cut5l,1,0),
           sixhs = ifelse(value >= cut6h,1,0),
           sixls = ifelse(value <= cut6l,1,0))%>%
    group_by(ids)%>%
    mutate(s1h = sum(onehs),
           s2h = sum(twohs),
           s3h = sum(threehs),
           s4h = sum(fourhs),
           s5h = sum(fivehs),
           s1l = sum(onels),
           s2l = sum(twols),
           s3l = sum(threels),
           s4l = sum(fourls),
           s5l = sum(fivels),
           s6h = sum(sixhs),
           s6l = sum(sixls))%>%
    group_by(ids)%>%
    mutate(p1h = ifelse(sum(onehs)/n()>=.05,1,0),
           p2h = ifelse(sum(twohs)/n()>=.05,1,0),
           p3h = ifelse(sum(threehs)/n()>=.05,1,0),
           p4h = ifelse(sum(fourhs)/n()>=.05,1,0),
           p5h = ifelse(sum(fivehs)/n()>=.05,1,0),
           p6h = ifelse(sum(sixhs)/n()>=.05,1,0),
           p1l = ifelse(sum(onels)/n()>=.05,1,0),
           p2l = ifelse(sum(twols)/n()>=.05,1,0),
           p3l = ifelse(sum(threels)/n()>=.05,1,0),
           p4l = ifelse(sum(fourls)/n()>=.05,1,0),
           p5l = ifelse(sum(fivels)/n()>=.05,1,0),
           p6l = ifelse(sum(sixls)/n()>=.05,1,0))%>%
    group_by(ids)%>%
    mutate(numst = n())%>%
    group_by(ids,strain)%>%
    mutate(cuts = ifelse(sixhs == 1 & ((s6h + s5h + s4h )/numst) <= .05 & (s5h==0|s4h==0), "remove",
                         ifelse(sixls == 1 & ((s6l + s5l +s4l)/numst) <= .05 & (s5l==0|s4l==0), "remove", "keep")),
           cuts1 = ifelse( ((s6h + s5h + s4h )/numst) >= .05 & s6h>=1&s5h>=1&s4h>=1, "keep",
                           ifelse((s5h>=1&s4h>=1&s3h>=1&s2h>=1&s1h>=1) | (s5h>=1&s3h>=1&s2h>=1&s1h>=1)|(s5h>=1&s4h>=1&s2h>=1&s1h>=1)| (s5h>=1&s4h>=1&s3h>=1&s1h>=1) |(s5h>=1&s4h>=1&s3h>=1&s2h>=1),"keep",
                                  ifelse((s5h>=1&s4l>=1&s3l>=1&s2l>=1&s1l>=1) |(s5h>=1&s3l>=1&s2l>=1&s1l>=1) | (s5h>=1&s4l>=1&s2l>=1&s1l>=1)|(s5h>=1&s4l>=1&s3l>=1&s1l>=1) | (s5h>=1&s4l>=1&s3l>=1&s2l>=1) ,"keep",
                                         ifelse(fivehs == 1 & ((s6h+ s4h + s5h + s3h)/numst) <= .05, "remove",
                                                ifelse(fivels == 1 & (s6l +s4l + s5l + s3l)/numst <= .05, "remove","keep"))))),
           cuts2 = ifelse( ((s6h + s5h + s4h )/numst) >= .05 & s6h>=1&s5h>=1&s4h>=1, "keep",
                           ifelse((s4h>=1&s3h>=1&s2h>=1&s1h>=1) | (s4h>=1&s2h>=1&s1h>=1)| (s4h>=1&s3h>=1&s1h>=1) |(s4h>=1&s3h>=1&s2h>=1),"keep",
                                  ifelse((s4l>=1&s3l>=1&s2l>=1&s1l>=1) | (s4l>=1&s2l>=1&s1l>=1)|(s4l>=1&s3l>=1&s1l>=1) | (s4l>=1&s3l>=1&s2l>=1) ,"keep",
                                         ifelse(fourhs == 1 & fivehs == 0 & (s5h + s4h + s3h + s2h)/numst <= .05, "remove",
                                                ifelse(fourls == 1  & fivels == 0 & (s5l + s4l + s3l + s2l)/numst <= .05, "remove","keep"))))))%>%
    filter(cuts != "remove")%>%
    filter(cuts1 != "remove")%>%
    filter(cuts2 != "remove")%>%
    select(drug,strain,contCond,trait, value, ids)
  
  return(cutPheno)
}

