# this script is used to perform the amr transformation. 
# it will add the larger and boarder group of amr ( based on hydrolysis mechanism and antibiotic group)
# the input is the output of staramr 
# the output will be an excel file with two tabs
# written and maintained by Tung Trinh 
# May 12th 2025 
# version 1.0 
# For more information, please contact to tungts@oucru.org 
####################################################################################################
amr_transform <- function(input,outdir = "path to the output folder",optional_text = NULL){
  ##### required packages and built-in function #####
  require(dplyr)
  require(stringr)
  require(readxl)
  require(tidyr)
  require(data.table)
  require(rJava)
  require(openxlsx)
  ##### get the AMR database #####
  amrdb <- readRDS("Data/AMRdb/amr_genotype_db_updated_2025-04-10.rds")
  ##### transform the data #####
  # create the id 
  input <- input %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+")) 
  # %>% mutate(seqid = if_else(is.na(seqid),"A3-12",seqid))
  # select variables to include into the analysis 
  # genotype 
  amr <- input %>% select(seqid, Genotype,Scheme) %>% separate_longer_delim(Genotype,delim = ",") %>% 
    mutate(Genotype = str_replace(Genotype," ",""))
  # merge with the amr database 
  amr <- amrdb %>% select(gene_name,amr_gene_group,hydrolysis_group) %>% distinct(.,.keep_all = T) %>% left_join(amr,., by = c("Genotype" = "gene_name")) %>% 
    distinct(.,.keep_all = T)
  # transform the list from long to wide 
  amr <- amr %>% mutate(presence = T)
  
  # work on hydrolysis group
  hydrolysis <- data.table::dcast(as.data.table(amr),seqid ~ hydrolysis_group, value.var = c("presence"))
  hydrolysis <- hydrolysis[,-2]
  id <- hydrolysis[,1]
  hydrolysis_tw <- hydrolysis[,-1]
  hydrolysis_tw[hydrolysis_tw >  0] <- 1
  hydrolysis <- cbind(id,hydrolysis_tw) %>% as_tibble()
  # work on amr group
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
  addWorksheet(wb,"hydrolysis group")
  # add the data into the worksheet 
  writeData(wb,"hydrolysis group",hydrolysis)
  # add worksheet
  addWorksheet(wb,"amr group")
  writeData(wb,"amr group",amr_group)
  optional_text = NULL
  openxlsx::saveWorkbook(wb, file = paste0(outdir,"/transformed_amr_",optional_text,Sys.Date(),".xlsx"), overwrite = T)
  
  
}