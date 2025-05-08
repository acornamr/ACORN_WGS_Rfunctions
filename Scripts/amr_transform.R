# this script is used to prepare the data for the amr annotation in ACORN WGS analysis
# written and maintained by Tung Trinh
# ver 1.0
# April 10th 2025
# For more information, please contact to tungts@oucru.org
####################################################################################################
##### required packages and built-in functions #####
library(dplyr)
library(stringr)
library(readxl)
library(tidyr)
library(data.table)
library(rJava)
library(openxlsx)
##### load data #####
# staramr
# staramr <- read_excel("Data/assemblies/Hanoi/staramr_ecol/staramr_ecol_results.xlsx")
# staramr <- read_excel("Data/assemblies/Hanoi/staramr_staph/results.xlsx")
# staramr <- read_excel("Data/assemblies/Hanoi/staramr_kleb_plus/results.xlsx")
staramr <- read_excel("Data/assemblies/Hanoi/staramr_acine/results.xlsx")
# amr database 
amrdb <- readRDS("Data/AMRdb/amr_genotype_db_updated_2025-04-10.rds")
##### transform the data #####
# create the id 
staramr <- staramr %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+")) 
# %>% mutate(seqid = if_else(is.na(seqid),"A3-12",seqid))
# select variables to include into the analysis 
# genotype 
amr <- staramr %>% select(seqid, Genotype,Scheme) %>% separate_longer_delim(Genotype,delim = ",") %>% 
  mutate(Genotype = str_replace(Genotype," ",""))
# merge with the amr database 
amr <- amrdb %>% select(gene_name,amr_gene_group,hydrolysis_group) %>% distinct(.,.keep_all = T) %>% left_join(amr,., by = c("Genotype" = "gene_name")) %>% 
  distinct(.,.keep_all = T)
# transform the list from long to wide 
amr <- amr %>% mutate(presence = T)

hydrolysis <- data.table::dcast(as.data.table(amr),seqid ~ hydrolysis_group, value.var = c("presence"))

hydrolysis <- hydrolysis[,-2]
# colnames(hydrolysis)[-1] <- paste0("hydrolysis_group.",colnames(hydrolysis))[-1]
id <- hydrolysis[,1]
hydrolysis_tw <- hydrolysis[,-1]

hydrolysis_tw[hydrolysis_tw >  0] <- 1
hydrolysis <- cbind(id,hydrolysis_tw) %>% as_tibble()

amr_group <- data.table::dcast(as.data.table(amr),seqid ~ amr_gene_group, value.var = c("presence"))

amr_group <- amr_group[,-2]
id <- amr_group[,1]
amr_group_tw <- amr_group[,-1]

amr_group_tw[amr_group_tw >  0] <- 1
amr_group <- cbind(id,amr_group_tw) %>% as_tibble()

### write out the data
# create a workbook
wb <- createWorkbook()
# add worksheet
# addWorksheet(wb,"ecoli - hydrolysis group")
# addWorksheet(wb,"staph - hydrolysis group")
# addWorksheet(wb,"kleb - hydrolysis group")
addWorksheet(wb,"acine - hydrolysis group")
# add the data into the worksheet 
# writeData(wb,"ecoli - hydrolysis group",hydrolysis)
# writeData(wb,"staph - hydrolysis group",hydrolysis)
# writeData(wb,"kleb - hydrolysis group",hydrolysis)
writeData(wb,"acine - hydrolysis group",hydrolysis)
# add worksheet
# addWorksheet(wb,"staph - amr group")
# addWorksheet(wb,"kleb - amr group")
addWorksheet(wb,"acine - amr group")
# add the data into the worksheet 
# writeData(wb,"kleb - amr group",amr_group)
writeData(wb,"acine - amr group",amr_group)


optional_text = NULL
# save the workbook
# openxlsx::saveWorkbook(wb, file = paste0("Data/acorn_ecoli_transformed_amr_",optional_text,Sys.Date(),".xlsx"), overwrite = T)
# openxlsx::saveWorkbook(wb, file = paste0("Data/acorn_staph_transformed_amr_",optional_text,Sys.Date(),".xlsx"), overwrite = T)
# openxlsx::saveWorkbook(wb, file = paste0("Data/acorn_kleb_transformed_amr_",optional_text,Sys.Date(),".xlsx"), overwrite = T)
openxlsx::saveWorkbook(wb, file = paste0("Data/acorn_acine_transformed_amr_",optional_text,Sys.Date(),".xlsx"), overwrite = T)


