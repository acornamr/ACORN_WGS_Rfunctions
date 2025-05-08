# this script is used for metadata checkup an clean up
# written and maintained by Tung Trinh 
# ver 1.0
# Feb 27th 2025
# For more information, please contact to tungts@oucru.org
####################################################################################################
##### required packages and built-in functions ######
library(readxl)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
"%ni%" <- Negate("%in%")
##### read the data #####
wgs <- read_excel("Data/ACORN2_VN001_WGS_Cumulative_V2_AKP.xlsx")
whole_dat  <- read_excel("Data/acorn_data_2025-01-22_04H35.xlsx", sheet = "acorn_dta")
staramr <- read_excel("Data/assemblies/Hanoi/staramr_kleb_plus/results.xlsx")
##### process data #####
# first get the acorn id from WGS 
target_id <- wgs$ACORN
# filter the whole data by selected id 
sub_dat <- whole_dat %>% filter(acorn_id %in% target_id)
### create the date range 
# set day 0 as the date of enrolment 
# calculate day
sub_dat <- sub_dat %>% mutate(start_date = as.numeric(as.Date(date_admission) - as.Date(date_enrolment) ),
                   hospitalised_date = as.numeric(as.Date(date_hospitalisation) - as.Date(date_enrolment)),
                   end_date = as.numeric(as.Date(ho_discharge_date) - as.Date(date_enrolment)),
                   death_date = as.numeric(as.Date(d28_death_date) - as.Date(date_enrolment)),
                   specimen_date = as.numeric(as.Date(specdate) - as.Date(date_enrolment)),
                   onset_date = as.numeric(as.Date(hai_date_symptom_onset) - as.Date(date_enrolment)),
                   ab1_start = as.numeric(as.Date(bsi_antibiotic1_startdate) - as.Date(date_enrolment)),
                   ab2_start = as.numeric(as.Date(bsi_antibiotic2_startdate) - as.Date(date_enrolment)),
                   ab3_start = as.numeric(as.Date(bsi_antibiotic3_startdate) - as.Date(date_enrolment)),
                   ab4_start = as.numeric(as.Date(bsi_antibiotic4_startdate) - as.Date(date_enrolment)),
                   ab5_start = as.numeric(as.Date(bsi_antibiotic5_startdate) - as.Date(date_enrolment)),
                   ab1_end = as.numeric(as.Date(bsi_antibiotic1_enddate) - as.Date(date_enrolment)),
                   ab2_end = as.numeric(as.Date(bsi_antibiotic2_enddate) - as.Date(date_enrolment)),
                   ab3_end = as.numeric(as.Date(bsi_antibiotic3_enddate) - as.Date(date_enrolment)),
                   ab4_end = as.numeric(as.Date(bsi_antibiotic4_enddate) - as.Date(date_enrolment)),
                   ab5_end = as.numeric(as.Date(bsi_antibiotic5_enddate) - as.Date(date_enrolment))
                   )

# subset the data for timeline chart 
timeline <- sub_dat %>% select(redcap_id, start_date,hospitalised_date, end_date, death_date, specimen_date,onset_date,ends_with("_start"),ends_with("_end"))%>% distinct(., .keep_all = T)
# transform the specimen date base the number of specimen times
timeline_ <- timeline %>% group_by(redcap_id) %>% mutate(times_specimen = n()) %>% ungroup %>% arrange((redcap_id)) %>% 
  group_by(redcap_id) %>% mutate(specimen_time_ = 1:max(times_specimen)) %>% ungroup %>% spread(.,key = "specimen_time_",value = "specimen_date") %>% 
  filter(times_specimen > 1) 


colnames(timeline_)[18:24] <- paste0("date_specimen",colnames(timeline_)[18:24])

acorn_id <- unique(timeline_$redcap_id)

jpeg(filename = "Figures/patient_timeline_1.jpg", res = 300, unit = "in",width = 10,height = 7)
ggplot(data = timeline_ %>% filter(redcap_id %in% acorn_id[1:16]), aes(y = 0, x = start_date))+
  geom_point(color = "#a65628", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = end_date), color = "#377eb8", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = death_date), color = "black", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = hospitalised_date), color = "#4daf4a", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = onset_date), color = "#e41a1c", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen1), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen2), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen3), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen4), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen5), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen6), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen7), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_segment(aes(x = ab1_start, y = 0.2, yend = 0.2, xend = ab1_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab2_start, y = -0.2, yend = -0.2, xend = ab2_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab3_start, y = 0.35, yend = 0.35, xend = ab3_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab4_start, y = -0.35, yend = -0.35, xend = ab4_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab5_start, y = 0.6, yend = 0.6, xend = ab5_end),size = 2, color = "darkgreen", alpha = 0.8)+
  facet_wrap(~redcap_id)+xlab("time since enrollment (days)")+ylab("")+ylim(-1,1)+xlim(-20,50)+
  theme_classic()+
  theme_classic()+
  theme(axis.text.y=element_blank())
dev.off()

jpeg(filename = "Figures/patient_timeline_2.jpg", res = 300, unit = "in",width = 10,height = 7)
ggplot(data = timeline_ %>% filter(redcap_id %in% acorn_id[17:32]), aes(y = 0, x = start_date))+
  geom_point(color = "#a65628", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = end_date), color = "#377eb8", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = death_date), color = "black", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = hospitalised_date), color = "#4daf4a", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = onset_date), color = "#e41a1c", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen1), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen2), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen3), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen4), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen5), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen6), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen7), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_segment(aes(x = ab1_start, y = 0.2, yend = 0.2, xend = ab1_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab2_start, y = -0.2, yend = -0.2, xend = ab2_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab3_start, y = 0.35, yend = 0.35, xend = ab3_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab4_start, y = -0.35, yend = -0.35, xend = ab4_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab5_start, y = 0.6, yend = 0.6, xend = ab5_end),size = 2, color = "darkgreen", alpha = 0.8)+
  facet_wrap(~redcap_id)+xlab("time since enrollment (days)")+ylab("")+ylim(-1,1)+xlim(-20,50)+
  theme_classic()+
  theme_classic()+
  theme(axis.text.y=element_blank())
dev.off()

jpeg(filename = "Figures/patient_timeline_3.jpg", res = 300, unit = "in",width = 10,height = 7)
ggplot(data = timeline_ %>% filter(redcap_id %in% acorn_id[33:48]), aes(y = 0, x = start_date))+
  geom_point(color = "#a65628", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = end_date), color = "#377eb8", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = death_date), color = "black", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = hospitalised_date), color = "#4daf4a", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = onset_date), color = "#e41a1c", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen1), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen2), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen3), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen4), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen5), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen6), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen7), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_segment(aes(x = ab1_start, y = 0.2, yend = 0.2, xend = ab1_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab2_start, y = -0.2, yend = -0.2, xend = ab2_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab3_start, y = 0.35, yend = 0.35, xend = ab3_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab4_start, y = -0.35, yend = -0.35, xend = ab4_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab5_start, y = 0.6, yend = 0.6, xend = ab5_end),size = 2, color = "darkgreen", alpha = 0.8)+
  facet_wrap(~redcap_id)+xlab("time since enrollment (days)")+ylab("")+ylim(-1,1)+xlim(-20,50)+
  theme_classic()+
  theme(axis.text.y=element_blank())
dev.off()

jpeg(filename = "Figures/patient_timeline_4.jpg", res = 300, unit = "in",width = 10,height = 7)
ggplot(data = timeline_ %>% filter(redcap_id %in% acorn_id[49:60]), aes(y = 0, x = start_date))+
  geom_point(color = "#a65628", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = end_date), color = "#377eb8", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = death_date), color = "black", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = hospitalised_date), color = "#4daf4a", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = onset_date), color = "#e41a1c", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen1), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen2), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen3), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen4), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen5), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen6), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen7), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_segment(aes(x = ab1_start, y = 0.2, yend = 0.2, xend = ab1_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab2_start, y = -0.2, yend = -0.2, xend = ab2_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab3_start, y = 0.35, yend = 0.35, xend = ab3_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab4_start, y = -0.35, yend = -0.35, xend = ab4_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab5_start, y = 0.6, yend = 0.6, xend = ab5_end),size = 2, color = "darkgreen", alpha = 0.8)+
  facet_wrap(~redcap_id)+xlab("time since enrollment (days)")+ylab("")+ylim(-1,1)+xlim(-20,50)+
  theme_classic()+
  theme(axis.text.y=element_blank())
dev.off()

jpeg(filename = "Figures/patient_timeline_5.jpg", res = 300, unit = "in",width = 10,height = 7)
ggplot(data = timeline_ %>% filter(redcap_id %in% acorn_id[61:72]), aes(y = 0, x = start_date))+
  geom_point(color = "#a65628", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = end_date), color = "#377eb8", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = death_date), color = "black", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = hospitalised_date), color = "#4daf4a", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = onset_date), color = "#e41a1c", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen1), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen2), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen3), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen4), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen5), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen6), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen7), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_segment(aes(x = ab1_start, y = 0.2, yend = 0.2, xend = ab1_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab2_start, y = -0.2, yend = -0.2, xend = ab2_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab3_start, y = 0.35, yend = 0.35, xend = ab3_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab4_start, y = -0.35, yend = -0.35, xend = ab4_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab5_start, y = 0.6, yend = 0.6, xend = ab5_end),size = 2, color = "darkgreen", alpha = 0.8)+
  facet_wrap(~redcap_id)+xlab("time since enrollment (days)")+ylab("")+ylim(-1,1)+xlim(-20,50)+
  theme_classic()+
  theme(axis.text.y=element_blank())
dev.off()

jpeg(filename = "Figures/patient_timeline_6.jpg", res = 300, unit = "in",width = 10,height = 7)
ggplot(data = timeline_ %>% filter(redcap_id %in% acorn_id[73:81]), aes(y = 0, x = start_date))+
  geom_point(color = "#a65628", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = end_date), color = "#377eb8", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = death_date), color = "black", alpha = 0.8,size = 2)+
  geom_point(aes(y = 0, x = hospitalised_date), color = "#4daf4a", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = onset_date), color = "#e41a1c", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen1), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen2), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen3), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen4), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen5), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen6), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_jitter(aes(y = 0, x = date_specimen7), color = "#ff7f00", alpha = 0.8,size = 2)+
  geom_segment(aes(x = ab1_start, y = 0.2, yend = 0.2, xend = ab1_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab2_start, y = -0.2, yend = -0.2, xend = ab2_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab3_start, y = 0.35, yend = 0.35, xend = ab3_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab4_start, y = -0.35, yend = -0.35, xend = ab4_end),size = 2, color = "darkgreen", alpha = 0.8)+
  geom_segment(aes(x = ab5_start, y = 0.6, yend = 0.6, xend = ab5_end),size = 2, color = "darkgreen", alpha = 0.8)+
  facet_wrap(~redcap_id)+xlab("time since enrollment (days)")+ylab("")+ylim(-1,1)+xlim(-20,50)+
  theme_classic()+
  theme(axis.text.y=element_blank())
dev.off()




