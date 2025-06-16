# this script is used to create the summary report of GHRU assembly pipeline 
# the input is the path to the folder contains report, tools_summary output from the pipeline
# Written and maintained by Tung Trinh
# June 14th 2025
# Ver 1.0 
# For more information, please contact to tungts@oucru.org
####################################################################################################
create_summary_report <- function(input_dir = "path to the input folder",
                                 output_dir = "path to the output folder",
                                 optional_text = NULL){
  ##### required packages and built-in function #####
  require(dplyr)
  require(data.table)
  require(openxlsx)
  require(rJava)
  ##### process the file #####
  # get the list of file 
  list_file <- list.files(path = input_dir, full.names = T)
  # create an empty list to write out the data
  dat_list <- list()
  # create an indicator 
  i = 1 
  ### for loop to read the file 
  for(file in list_file){
    # read file 
    dat_list[[i]] <- read.csv(file, header = T, sep = "\t")
    # accumulate the indicator 
    i = i + 1 
  }
  # bind the list 
  out <- rbindlist(dat_list, fill = T)
  ### output the file 
  # create wb 
  wb <- createWorkbook()
  # add worksheet
  addWorksheet(wb,"report")
  # add data to the report 
  writeData(wb,"report",out)
  openxlsx::saveWorkbook(wb, file = paste0(output_dir,"/summary_report_",optional_text,Sys.Date(),".xlsx"), overwrite = T)
  
}
