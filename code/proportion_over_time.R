####################################################
## James Hay <jhay@hsph.harvard.edu>
## Created: 2022-01-10
## Summary: Script to read in the tidy Ct data and find the proportion with Ct values < 30 over time, stratified
## by how soon after a last negative test individuals were detected

## Preamble
library(tidyverse)
library(zoo)
library(lubridate)
library(patchwork)

setwd("~/Documents/GitHub/Ct_Omicron")
min_date <- as.Date("2021-07-05")
low_ct_threshold <- 30
inconclusive_sensitivity <- FALSE

filename_base <- paste0("figures/preprint/Ct_threshold_",low_ct_threshold,"_inconclusives",inconclusive_sensitivity)
if(!file.exists(filename_base)) dir.create(filename_base)

if(inconclusive_sensitivity){
    load("data/ct_dat_subset_figure1_using_inconclusives.RData")
} else {
    load("data/ct_dat_subset_figure1.RData")
}

dat <- dat %>% filter(TestResult != "No sample")

## Find when individuals were detected relative to their last negative PCR result
dat_detections <- dat %>% 
    filter(DaysSinceDetection == 0) %>%
    group_by(PersonID) %>%
    mutate(DetectionSpeed = ifelse(DaysSinceNegative >= 2 | is.na(DaysSinceNegative),
                                   ">=2 days",ifelse(DaysSinceNegative < 2, "<=1 days",NA))) %>% 
    mutate(DetectionSpeed=ifelse(is.na(DetectionSpeed),
                                 ">=2 days",DetectionSpeed)) %>%
    select(PersonID, TestDate, DetectionSpeed, DaysSinceNegative)

#dat_detections <- classify_detection_speed(dat)
dat <- dat %>% left_join(dat_detections)%>% 
    fill(DetectionSpeed, .direction="down")

## Get infection histories as unique lineages
infection_histories <- dat %>% dplyr::select(PersonID, Lineage)  %>% 
    distinct() %>% 
    nest(all_lineages=c(Lineage))  %>%
    group_by(PersonID)
dat <- left_join(dat, infection_histories)


# Omicron -----------------------------------------------------------------
## Get only those infected with omicron at some point
omicron_ids <- infection_histories %>% 
    group_by(PersonID) %>% 
    filter("BA.1" %in% unlist(all_lineages) | "Suspected Omicron" %in% unlist(all_lineages)) %>% 
    dplyr::pull(PersonID)

## Subset data by those who had omicron at some point
## If suspected lineage is listed as Omicron but the sample collection date is before November 2021, then set lineage to None
#dat <- dat %>% mutate(Lineage=ifelse(Lineage %in% c("BA.1","Suspected Omicron") & TestDate < "2021-11-01", "None",Lineage))
dat_omi <- dat %>% filter(PersonID %in% omicron_ids)
dat_omi <- dat_omi %>% group_by(PersonID, CumulativeInfectionNumber) %>% 
    mutate(has_omicron = any(Lineage %in% c("BA.1","Suspected Omicron"))) %>% 
    filter(has_omicron == TRUE)

## Merge this indexing data subset and then get Ct values over the whole infection course
dat_omi_subset <- dat_omi %>% 
    ungroup() %>%
    left_join(dat_detections %>% dplyr::select(PersonID, DetectionSpeed,DaysSinceNegative)) %>% 
    dplyr::select(PersonID, Lineage, TestDate, TestResult, CtT1, CtT2, MostRecentDetection, DetectionSpeed,
                  DaysSinceNegative, Detected)%>%
    distinct()

## Only look at positives between -5 days and 15 days between the first positive
dat_omi_subset <- dat_omi_subset %>% group_by(PersonID) %>% filter(TestResult == "Positive" | TestDate >= (MostRecentDetection - 5))
## Days since first positive
dat_omi_subset <- dat_omi_subset %>% ungroup() %>% mutate(TimeSinceFirstPos=TestDate-MostRecentDetection)
## Only plot time points between 0 and 15 days post first positive
dat_omi_subset_tmp <- dat_omi_subset %>% filter(TimeSinceFirstPos >= 0,TimeSinceFirstPos <= 15) %>% 
filter(TestDate > as.numeric(as.Date("2021-11-01") - min_date))

## Sample sizes
dat_omi_subset_tmp %>% dplyr::select(PersonID, DetectionSpeed) %>% distinct() %>% group_by(DetectionSpeed) %>% tally() 
dat_omi_subset_tmp %>% distinct() %>% group_by(DetectionSpeed) %>% tally()

## Find proportion of Cts < 30
dat_omi_subset_tmp <- dat_omi_subset_tmp %>% mutate(CtT1=ifelse(is.na(CtT1),40,CtT1),CtT2=ifelse(is.na(CtT2),40,CtT2))
dat_omi_subset_tmp <- dat_omi_subset_tmp %>% mutate(low_ct1 = CtT1 < low_ct_threshold, low_ct2 = CtT2 < low_ct_threshold)

tmp <- dat_omi_subset_tmp %>% group_by(TimeSinceFirstPos, DetectionSpeed) %>% 
    summarize(n_low=sum(low_ct1,na.rm=TRUE),N=n()) %>%
    mutate(prop_low=n_low/N)

## Proportion with Ct < 30 on day 5 or later
dat_omi_subset_tmp %>% mutate(AfterDay5 = TimeSinceFirstPos >= 5) %>% 
    filter(AfterDay5 == TRUE) %>%
    group_by(PersonID, DetectionSpeed) %>% 
    summarize(AnyLowCtAfterDay5 = any(low_ct1)) %>%
    select(PersonID, AnyLowCtAfterDay5, DetectionSpeed) %>%
    group_by(DetectionSpeed) %>% 
    summarize(NLowDay5Onward=sum(AnyLowCtAfterDay5,na.rm=TRUE),N=n()) %>%
    mutate(PLowDay5Onward=NLowDay5Onward/N)

## Proportion with Ct < 30 at any point
dat_omi_subset_tmp %>% 
    group_by(PersonID, DetectionSpeed) %>% 
    summarize(AnyLowCt = any(low_ct1)) %>%
    select(PersonID, AnyLowCt, DetectionSpeed) %>%
    group_by(DetectionSpeed) %>% 
    summarize(NLowCT=sum(AnyLowCt,na.rm=TRUE),N=n()) %>%
    mutate(PLowAny=NLowCT/N)

omicron_p1_key <- c("<=1 days"="Frequent testing (≤1 days since last non-positive PCR); n=27",
                    ">=2 days"="Testing due to symptoms or contact \n(≥2 days since last non-positive PCR); n=70")
tmp$DetectionSpeed <- omicron_p1_key[tmp$DetectionSpeed]
dat_omi_subset_tmp$DetectionSpeed <- omicron_p1_key[dat_omi_subset_tmp$DetectionSpeed]
dat_omi_summary <- dat_omi_subset_tmp %>% group_by(DetectionSpeed, TimeSinceFirstPos) %>% summarize(mean_ct=mean(CtT1,na.rm=TRUE)) 
dat_omi_summary2 <- dat_omi_subset_tmp %>% filter(CtT1 < 40) %>% group_by(DetectionSpeed, TimeSinceFirstPos) %>% summarize(mean_ct=mean(CtT1,na.rm=TRUE)) 
## Plot trajectories with mean Ct trajectory
p_traj_omicron <-  ggplot()  + 
    geom_line(data=dat_omi_subset_tmp, aes(x=TimeSinceFirstPos,y=CtT1,group=PersonID),alpha=0.25,size=0.25,col="red") +
    geom_line(data=dat_omi_summary, aes(x=TimeSinceFirstPos,y=mean_ct,linetype="All tests"),size=0.75,col="red") + 
    geom_line(data=dat_omi_summary2, aes(x=TimeSinceFirstPos,y=mean_ct,linetype="Only positive"),size=0.75,col="red",linetype="dashed") + 
    scale_y_continuous(trans="reverse") + facet_wrap(~DetectionSpeed) + theme_minimal() + xlab("Days since first positive") + ylab("Ct value")+ theme(plot.background = element_rect(fill="white",color=NA)) +
    geom_hline(yintercept=low_ct_threshold,linetype="solid",col="black")+
    scale_linetype_manual(name="Mean of",values=c("All tests"="solid","Only positive"="dashed")) +
    theme(legend.position=c(0.8,0.8)) +
    ggtitle("Omicron")

## Plot proportion detectable over time
p_detect_omicron <- ggplot(tmp) + geom_line(aes(x=TimeSinceFirstPos, y=prop_low,col=DetectionSpeed)) + 
    scale_y_continuous(limits=c(0,1)) + 
    ylab(paste0("Proportion Ct<",low_ct_threshold)) + 
    xlab("Days since first positive") + 
    scale_color_manual(name="Time since last non-positive PCR",values=c("Frequent testing (≤1 days since last non-positive PCR); n=27"="dark green","Testing due to symptoms or contact \n(≥2 days since last non-positive PCR); n=70"="orange")) + theme_minimal() + theme(plot.background = element_rect(fill="white",color=NA),legend.position=c(0.7,0.8)) + 
    geom_vline(xintercept=5) + scale_x_continuous(expand=c(0,0), limits=c(0,15), breaks=seq(0,15,by=1))

p_detect_omicronN <- ggplot(tmp) + 
    geom_line(aes(x=TimeSinceFirstPos, y=N,col=DetectionSpeed,linetype="Total")) + 
    geom_line(aes(x=TimeSinceFirstPos, y=n_low,col=DetectionSpeed,linetype="N Ct<30")) + 
    ylab("N") + 
    xlab("Days since first positive") + 
    scale_y_continuous(breaks=seq(0,80,by=5)) +
    scale_color_manual(name=NULL,values=c("Frequent testing (≤1 days since last non-positive PCR); n=27"="dark green","Testing due to symptoms or contact \n(≥2 days since last non-positive PCR); n=70"="orange")) + 
    scale_linetype_manual(name=NULL,values=c("Total"="solid","N Ct<30"="dashed"))+
    theme_minimal() + 
    theme(plot.background = element_rect(fill="white",color=NA),legend.position="bottom", legend.text=element_text(size=7)) + 
    geom_vline(xintercept=5) + scale_x_continuous(expand=c(0,0), limits=c(0,15), breaks=seq(0,15,by=1))


# Delta -------------------------------------------------------------------
## Get only those infected with omicron at some point
delta_lineages <- unique(dat$Lineage)
library(data.table)
delta_lineages <- delta_lineages[delta_lineages %like% "AY" | delta_lineages == "B.1.617.2"]
delta_ids <- infection_histories %>% 
    group_by(PersonID) %>% 
    filter(any(delta_lineages %in% unlist(all_lineages) )) %>% 
    dplyr::pull(PersonID)

## Subset data by those who had Delta at some point
dat_delta <- dat %>% filter(PersonID %in% delta_ids)
dat_delta <- dat_delta %>% group_by(PersonID, CumulativeInfectionNumber) %>% 
    mutate(has_delta = any(Lineage %in% delta_lineages)) %>% 
    filter(has_delta == TRUE)

## Merge this indexing data subset and then get Ct values over the whole infection course
dat_delta_subset <- dat_delta %>% 
    left_join(dat_detections %>% dplyr::select(PersonID, DetectionSpeed,DaysSinceNegative)) %>% 
    dplyr::select(PersonID, Lineage, TestDate, TestResult, CtT1, CtT2, MostRecentDetection, DetectionSpeed,
                  DaysSinceNegative, Detected) %>%
    distinct()

dat_delta_subset <- dat_delta_subset %>% mutate(DetectionSpeed = ifelse(is.na(DetectionSpeed),">=2 days",DetectionSpeed ))
## Only look at positives between -5 days and 15 days between the first positive
dat_delta_subset <- dat_delta_subset %>% group_by(PersonID) %>% filter(TestResult == "Positive" | TestDate >= (MostRecentDetection - 5))
## Days since first positive
dat_delta_subset <- dat_delta_subset %>% ungroup() %>% mutate(TimeSinceFirstPos=TestDate-MostRecentDetection)
## Only plot time points between 0 and 15 days post first positive
dat_delta_subset_tmp <- dat_delta_subset %>% filter(TimeSinceFirstPos >= 0,TimeSinceFirstPos <= 15) 

## Sample sizes
## Note PersonID 429 only has one positive, so is filtered out around here
dat_delta_subset_tmp %>% dplyr::select(PersonID, DetectionSpeed) %>% distinct() %>% group_by(DetectionSpeed) %>% tally()
dat_delta_subset_tmp %>% distinct() %>% group_by(DetectionSpeed) %>% tally()

## Find proportion of Cts < 30
dat_delta_subset_tmp <- dat_delta_subset_tmp %>% mutate(CtT1=ifelse(is.na(CtT1),40,CtT1),CtT2=ifelse(is.na(CtT2),40,CtT2))
dat_delta_subset_tmp <- dat_delta_subset_tmp %>% mutate(low_ct1 = CtT1 < 30, low_ct2 = CtT2 < 30)

tmp_delta <- dat_delta_subset_tmp %>% group_by(TimeSinceFirstPos, DetectionSpeed) %>% 
    summarize(n_low=sum(low_ct1,na.rm=TRUE),N=n()) %>%
    mutate(prop_low=n_low/N)

delta_p1_key <- c("<=1 days"="Frequent testing (≤1 days since last non-positive PCR); n=6",
                    ">=2 days"="Testing due to symptoms or contact \n(≥2 days since last non-positive PCR); n=101")
tmp_delta$DetectionSpeed <- delta_p1_key[tmp_delta$DetectionSpeed]
dat_delta_subset_tmp$DetectionSpeed <- delta_p1_key[dat_delta_subset_tmp$DetectionSpeed]


dat_delta_summary <- dat_delta_subset_tmp %>% group_by(DetectionSpeed, TimeSinceFirstPos) %>% summarize(mean_ct=mean(CtT1,na.rm=TRUE)) 
dat_delta_summary2 <- dat_delta_subset_tmp %>% filter(CtT1 < 40) %>% group_by(DetectionSpeed, TimeSinceFirstPos) %>% summarize(mean_ct=mean(CtT1,na.rm=TRUE)) 

## Plot trajectories with mean Ct trajectory
p_traj_delta <- ggplot()  + 
    geom_line(data=dat_delta_subset_tmp, aes(x=TimeSinceFirstPos,y=CtT1,group=PersonID),alpha=0.25,size=0.25,col="blue") +
    geom_line(data=dat_delta_summary,aes(x=TimeSinceFirstPos,y=mean_ct,linetype="All tests"),size=0.75,col="blue") +
    geom_line(data=dat_delta_summary2,aes(x=TimeSinceFirstPos,y=mean_ct,linetype="Only positive"),size=0.75,col="blue") + 
    scale_y_continuous(trans="reverse") + 
    facet_wrap(~DetectionSpeed) + 
    theme_minimal() + 
    xlab("Days since first positive") + 
    ylab("Ct value")+ 
    theme(plot.background = element_rect(fill="white",color=NA)) +
    geom_hline(yintercept=low_ct_threshold,linetype="solid",col="black") +
    scale_linetype_manual(name="Mean of",values=c("All tests"="solid","Only positive"="dashed")) +
    theme(legend.position=c(0.8,0.8)) +
    ggtitle("Delta")

## Plot proportion detectable over time
p_detect_delta <- ggplot(tmp_delta) + geom_line(aes(x=TimeSinceFirstPos, y=prop_low,col=DetectionSpeed)) + 
    scale_y_continuous(limits=c(0,1)) +
    ylab(paste0("Proportion Ct<",low_ct_threshold)) + 
    xlab("Days since first positive") + 
    scale_color_manual(name="Time since last non-positive PCR",values=c("Frequent testing (≤1 days since last non-positive PCR); n=6"="dark green","Testing due to symptoms or contact \n(≥2 days since last non-positive PCR); n=101"="orange")) + theme_minimal() + theme(plot.background = element_rect(fill="white",color=NA),legend.position=c(0.7,0.8)) + 
    geom_vline(xintercept=5) + scale_x_continuous(expand=c(0,0), limits=c(0,15), breaks=seq(0,15,by=1))

p_detect_deltaN <- ggplot(tmp_delta) + 
    geom_line(aes(x=TimeSinceFirstPos, y=N,col=DetectionSpeed,linetype="Total")) + 
    geom_line(aes(x=TimeSinceFirstPos, y=n_low,col=DetectionSpeed,linetype="N Ct<30")) + 
    ylab("N") + 
    xlab("Days since first positive") + 
    scale_y_continuous(breaks=seq(0,100,by=5)) +
    scale_color_manual(name=NULL,values=c("Frequent testing (≤1 days since last non-positive PCR); n=6"="dark green","Testing due to symptoms or contact \n(≥2 days since last non-positive PCR); n=101"="orange")) + 
    scale_linetype_manual(name=NULL,values=c("Total"="solid","N Ct<30"="dashed"))+
    theme_minimal() + 
    theme(plot.background = element_rect(fill="white",color=NA),legend.position="bottom", legend.text=element_text(size=7)) + 
    geom_vline(xintercept=5) + scale_x_continuous(expand=c(0,0), limits=c(0,15), breaks=seq(0,15,by=1))


p_traj <- p_traj_omicron/p_traj_delta

ggsave(paste0(filename_base,"/prop_detectable_omicron.png"),p_detect_omicron,width=8,height=4,units="in",dpi=300)
ggsave(paste0(filename_base,"/N_detectable_omicron.png"),p_detect_omicronN,width=8,height=4,units="in",dpi=300)
ggsave(paste0(filename_base,"/trajectories_omicron.png"),p_traj_omicron,width=8,height=4,units="in",dpi=300)

ggsave(paste0(filename_base,"/prop_detectable_delta.png"),p_detect_delta,width=8,height=4,units="in",dpi=300)
ggsave(paste0(filename_base,"/N_detectable_delta.png"),p_detect_deltaN,width=8,height=4,units="in",dpi=300)
ggsave(paste0(filename_base,"/trajectories_delta.png"),p_traj_delta,width=8,height=4,units="in",dpi=300)



## Sample sizes
dat_omi_subset_tmp %>% dplyr::select(PersonID, DetectionSpeed) %>% distinct() %>% group_by(DetectionSpeed) %>% tally()
dat_delta_subset_tmp %>% dplyr::select(PersonID, DetectionSpeed) %>% distinct() %>% group_by(DetectionSpeed) %>% tally()

tmp %>% filter(TimeSinceFirstPos == 5)
tmp %>% filter(TimeSinceFirstPos == 6)
tmp %>% filter(TimeSinceFirstPos == 7)
tmp %>% filter(TimeSinceFirstPos == 8)
tmp %>% filter(TimeSinceFirstPos == 10)

tmp_delta %>% filter(TimeSinceFirstPos == 5)
tmp_delta %>% filter(TimeSinceFirstPos == 6)
tmp_delta %>% filter(TimeSinceFirstPos == 7)
tmp_delta %>% filter(TimeSinceFirstPos == 8)



## Plot trajectories with mean Ct trajectory
dat_traj_comb <- bind_rows(dat_omi_subset_tmp %>% mutate(Lineage="Omicron"),dat_delta_subset_tmp %>% mutate(Lineage="Delta"))
dat_traj_comb$Lineage <- factor(as.character(dat_traj_comb$Lineage), levels=unique(dat_traj_comb$Lineage))
dat_traj_summary <- dat_traj_comb %>% group_by(DetectionSpeed, TimeSinceFirstPos, Lineage) %>% summarize(mean_ct=mean(CtT1,na.rm=TRUE))%>% mutate(`Mean from`="All tests")
dat_traj_summary2 <- dat_traj_comb %>% filter(CtT1 < 40) %>% group_by(DetectionSpeed, TimeSinceFirstPos, Lineage) %>% summarize(mean_ct=mean(CtT1,na.rm=TRUE)) %>% mutate(`Mean from`="Only positives")
#dat_traj_summary <- bind_rows(dat_traj_summary, dat_traj_summary2)

dat_traj_summary$Lineage <- factor(dat_traj_summary$Lineage, levels=unique(dat_traj_summary$Lineage))

p_traj_comb <-  ggplot()  + 
    geom_hline(yintercept=low_ct_threshold,linetype="dotted",col="black") +
    geom_line(data=dat_traj_comb, aes(x=TimeSinceFirstPos,y=CtT1,group=PersonID,col=Lineage),alpha=0.25,size=0.25) +
    geom_line(data=dat_traj_summary,aes(x=TimeSinceFirstPos,y=mean_ct,col=Lineage),size=0.75) + 
    scale_y_continuous(trans="reverse") + 
    facet_wrap(~Lineage + DetectionSpeed,
               labeller = function (labels) {
                   labels <- lapply(labels, as.character)
                   list(do.call(paste, c(labels, list(sep = " infection; "))))
               }) + 
    theme_minimal() + xlab("Days since first positive") + ylab("Ct value")+ 
    theme(plot.background = element_rect(fill="white",color=NA),legend.position="none",
          panel.spacing=unit(0,"lines"),strip.text=element_text(size=8),
          axis.title.x=element_blank()) +
    scale_color_manual(values=c("Delta"="blue","Omicron"="red")) +
    labs(tag="A") +
    theme(plot.tag = element_text(face="bold"))
p_traj_comb

p_fig1 <- (p_traj_comb / p_detect_omicron + labs(tag="B")+
        theme(plot.tag = element_text(face="bold"))) + plot_layout(heights=c(1.7,1),nrow=2)

ggsave(paste0(filename_base,"/trajectories.png"),p_traj_comb,width=9,height=7,units="in",dpi=300)
ggsave(paste0(filename_base,"/fig1.png"),p_fig1,width=9,height=10,units="in",dpi=300)


write_csv(tmp %>% arrange(DetectionSpeed, TimeSinceFirstPos), paste0(filename_base,"/omicron_tableS1.csv"))
write_csv(tmp_delta %>% arrange(DetectionSpeed, TimeSinceFirstPos), paste0(filename_base,"/delta_tableS1.csv"))