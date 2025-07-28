# this script is used to merge excel files into one uniform dataset
# these excel files should contains a same format and structure
# writen and maintained by Tung Trinh
# version 1.0 
# July 27th 2025
# For more information, please contact to tungts@oucru.org
####################################################################################################
merge_excel <- function(list_file, sheet_name = NULL){
  ##### required packages and built-in functions #####
  require(dplyr)
  require(data.table)
  require(readxl)
  ##### read the file #####
  # create an empty list 
  list_excel <- list()
  # create an indicator 
  i = 1 
  # for loop to read the file 
  for(file in list_file){
    # read the excel file 
    # check the sheet_name 
    if(!is.null(sheet_name)){
      # write out into the list 
      list_excel[[i]] <- read_excel(file, sheet = sheet_name)
    }else{
      list_excel[[i]] <- read_excel(file)
    }
    # accumulate the indicator 
    i = i + 1 
  }
  # merge the excel files 
  out <- rbindlist(list_excel,fill = T)
  ##### output the file #####
  return(out)
}