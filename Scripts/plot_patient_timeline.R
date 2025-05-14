# this script is used to create the patient timeline from the metadata.
# thus, the input is the metadata file from the the acorn database
# Written and maintained by Tung Trinh
# May 7th 2025
# Version 1.1
# May 8th 2025
# For more information, please contact to tungts@oucru.org 
####################################################################################################
plot_patient_timeline <- function(input = "input data", output_dir = "", optional_text = NULL, start_day = -20, end_day = 50){
  ##### required packages and built-in functions #####
  require(dplyr)
  require(ggplot2)
  require(stringr)
  require(reshape2)
  require(tidyr)
  "%ni%" <- Negate("%in%")
  ##### checking the data #####
  ##### transform the data #####
  # calculate the interval with date of enrolment as day 0 
  input <- input %>% mutate(start_date = as.numeric(as.Date(date_admission) - as.Date(date_enrolment) ),
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
  # subset the data 
  # subset the data for timeline chart 
  timeline <- input %>% select(redcap_id, start_date,hospitalised_date, end_date, death_date, specimen_date,onset_date,ends_with("_start"),ends_with("_end"))%>% distinct(., .keep_all = T)
  # calculate the number of specimen times then transform the date of sampling from long to wide 
  timeline_ <- timeline %>% group_by(redcap_id) %>% mutate(times_specimen = n()) %>% ungroup %>% arrange((redcap_id)) %>% 
    group_by(redcap_id) %>% mutate(specimen_time_ = 1:max(times_specimen)) %>% ungroup %>% spread(.,key = "specimen_time_",value = "specimen_date") %>% 
    filter(times_specimen > 1) 
  # assign column name of variables contains the date of sampling
  colnames(timeline_)[18:24] <- paste0("date_specimen",colnames(timeline_)[18:24])
  # get the id 
  # acorn_id <- unique(timeline_$redcap_id)
  plot <- ggplot(data = timeline_ %>% filter(redcap_id %in% acorn_id[1:16]), aes(y = 0, x = start_date))+
    geom_point(color = "#a65628", alpha = 0.8,size = 2)+
    geom_point(aes(y = 0, x = end_date), color = "#377eb8", alpha = 0.8,size = 2)+
    geom_point(aes(y = 0, x = death_date), color = "black", alpha = 0.8,size = 2)+
    geom_point(aes(y = 0, x = hospitalised_date), color = "#4daf4a", alpha = 0.8,size = 2)+
    geom_point(aes(y = 0.1, x = onset_date), color = "#e41a1c", alpha = 0.8,size = 2)+
    geom_point(aes(y = 0.2, x = date_specimen1), color = "#ff7f00", alpha = 0.8,size = 2)+
    geom_point(aes(y = 0.2, x = date_specimen2), color = "#ff7f00", alpha = 0.8,size = 2)+
    geom_point(aes(y = 0.2, x = date_specimen3), color = "#ff7f00", alpha = 0.8,size = 2)+
    geom_point(aes(y = 0.2, x = date_specimen4), color = "#ff7f00", alpha = 0.8,size = 2)+
    geom_point(aes(y = 0.2, x = date_specimen5), color = "#ff7f00", alpha = 0.8,size = 2)+
    geom_point(aes(y = 0.2, x = date_specimen6), color = "#ff7f00", alpha = 0.8,size = 2)+
    geom_point(aes(y = 0.2, x = date_specimen7), color = "#ff7f00", alpha = 0.8,size = 2)+
    geom_segment(aes(x = ab1_start, y = -0.1, yend = -0.1, xend = ab1_end),size = 2, color = "darkgreen", alpha = 0.8)+
    geom_segment(aes(x = ab2_start, y = -0.15, yend = -0.15, xend = ab2_end),size = 2, color = "darkgreen", alpha = 0.8)+
    geom_segment(aes(x = ab3_start, y = -0.2, yend = -0.2, xend = ab3_end),size = 2, color = "darkgreen", alpha = 0.8)+
    geom_segment(aes(x = ab4_start, y = -0.25, yend = -0.25, xend = ab4_end),size = 2, color = "darkgreen", alpha = 0.8)+
    geom_segment(aes(x = ab5_start, y = -0.3, yend = -0.3, xend = ab5_end),size = 2, color = "darkgreen", alpha = 0.8)+
    facet_wrap(~redcap_id)+xlab("time since enrollment (days)")+ylab("")+ylim(-0.4,0.4)+xlim(start_day,end_day)+
    theme_classic()+
    theme(axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid.major.x = element_line(linetype = 2,colour = "gray90",size = 0.5)
            )
  ##### output ######
  jpeg()
}
