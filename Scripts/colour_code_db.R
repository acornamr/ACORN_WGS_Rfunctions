# this script is used to create the colour code for each catergorical variable.
# written and maintained by Tung Trinh
# Aug 7th 2025
# ver 1.0
# For more information, please contact to tungts@oucru.org 
####################################################################################################
##### required packages and built-in functions #####
library(dplyr)
library(data.table)
##### create database 
# source of sample 
source <- data.frame(variable = "surveillance_category",
                     value = c("CAI","HAI","HCAI",NA),
                     colour = c("#1F8F99FF","#D60C00FF","#6B58EFFF","gray100"))
# sample type 
sample_type <- data.frame(variable = "specgroup",
                          value = c("Blood","CSF","Lower respiratory tract specimen","Other specimens",
                                    "Pleural fluid","Sterile fluids","Stool","Urine",NA),
                          colour = c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#ffff33","#e6ab02","#a6761d","gray100"))
# outcome 
outcome <- data.frame(variable = "d28_status",
                      value = c("Alive - completely recovered","Alive - not back to normal activities","Dead",                                 
                                "Unable to contact"),
                      colour = c("#FCFDBFFF","#F1605DFF","#000004FF","gray50"))
# severity 
# create function
clr_func <- colorRampPalette(c("#bca89f","#a63603"))
# generate color 
severity_pal <- clr_func(6)
severity <- data.frame(variable = "clinical_severity_score",
                       value = c("0","1","2","3","4","6"),
                       colour = severity_pal)
# merge the dataset 
output <- rbindlist(list(outcome,sample_type,source,severity))
##### output the dataset 
saveRDS(output,file = paste0("Data/colour_code_db_",Sys.Date(),".rds"))

