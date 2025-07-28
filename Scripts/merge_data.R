# This script is used to merge the sequencing data with the metadata from ACORN's REDCAP database
# Written and maintained by Tung Trinh
# July 23rd 2025
# Version 1.1
# For more information, please contact to tungts@oucru.org 
####################################################################################################
merge_data <- function(input_redcap = "path to the input redcap data",
                       input_seq = "path to the input sequence metadata",
                       seq_id = "seq id from the input_seq",
                       sample_id = "sample id from the input_seq",
                       optional_text = NULL,
                       output_dir = "path to output directory"){
  ##### required packages and built-in library #####
  require(readxl)
  require(dplyr)
  require(stringr)
  require(rJava)
  require(openxlsx)
  require(purrr)
  ##### read the data #####
  # redcap 
  redcap <- read_excel(input_redcap, sheet = "acorn_dta")
  # seq
  seq <- read_excel(input_seq)
  ##### transform the data #####
  ### transform seq data into usable form 
  ## find the row that contain the acorn id 
  # find the acorn id 
  cor_acornid <- seq %>% map(~grep("acornid",.)) %>% unlist 
  # get the column name 
  colnames(seq) <- seq[cor_acornid,] 
  # exclude unneccsary column 
  if(any(grepl("-",colnames(seq)))){
     seq <- seq[,-grep("-",colnames(seq))] 
  }
  # exclude the column row 
  seq <- seq %>% slice(-cor_acornid)
  # change the variable name of certain variables 
  seq <- seq %>% rename(seq_id = matches(seq_id),
                 specid = matches(sample_id))
  
  ##### merge with the acorn data #####
  #merge
  out <- left_join(seq,redcap, by = c("specid" = "specid")) %>% distinct(specid,.keep_all = T)
  ##### output the file #####
  ### rds file 
  saveRDS(out,file = paste0(output_dir,"/merged_metadata_",optional_text,Sys.Date(),".rds"))
  ### excel file 
  # create a workbook
  wb <- createWorkbook()
  # add worksheet
  addWorksheet(wb,"data")
  # add the data into the worksheet 
  writeData(wb,"data",out)
  openxlsx::saveWorkbook(wb, file = paste0(output_dir,"/merged_metadata_",optional_text,Sys.Date(),".xlsx"), overwrite = T)
}