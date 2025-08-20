# this script is used to create the colour code for each catergorical variable.
# written and maintained by Tung Trinh
# Aug 7th 2025
# ver 1.0
# For more information, please contact to tungts@oucru.org 
####################################################################################################
##### required packages and built-in functions #####
library(dplyr)
##### create database 
# source of sample 
source <- data.frame(variable = "surveillance_category",
                     value = c("CAI","HAI","HCAI",NA),
                     colour = c("#1F8F99FF","#D60C00FF","#6B58EFFF","gray80"))
# sample type 
sample_type <- data.frame(variable = "specgroup",
                          value = c("Blood","CSF","Lower respiratory tract specimen","Other specimens",
                                    "Pleural fluid","Sterile fluids","Stool","Urine",NA),
                          colour = c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d","#666666","gray80"))
# outcome 
outcome <- data.frame(variable = "d28_status",
                      value = c("Alive - completely recovered","Alive - not back to normal activities","Dead",                                 
                                "Unable to contact"),
                      colour = c("#FCFDBFFF","#F1605DFF","#000004FF","gray70"))
# 