# this function is used to define the colour code for each categorical group 
# written and maintained by Tung Trinh
# version 1.0 
# Aug 7th 2025
# For more information, please contact to tungts@oucru.org
####################################################################################################
get_colour_code <- function(input_data = "",
                            target_variable = ""){
  ##### required packages and built-in function #####
  require(dplyr)
  ##### get the database #####
  # find the latest colour code database 
  list_file <- list.files(path = "Data",pattern = "colour_code_db", full.names = T)
  file <- file.info(list_file) %>% arrange(mtime) %>% slice(1) %>% rownames()
  # read the file 
  db <- readRDS(file = file)
  ##### get the color code 
  if(!is.null(target_variable)){
    if((target_variable %in% colnames(input_data)) == F){
      stop(paste("Could not find",target_variable,"in the dataset"))
    }else{
      print(paste("The target variable",target_variable,"found in the dataset"))
    }
  }
  # rename the target variable 
  input_data <- input_data %>% rename(target_variable = matches(target_variable))
  # sort the database 
  target_db <- db %>% filter(variable == target_variable)
  # get the list of unique value 
  list_value <- input_data %>% pull(target_variable) %>% unique() %>% sort()
  # get the colour code 
  output <- target_db %>% filter(value %in% list_value) %>% pull(colour)
  ##### output the colour code 
  return(output)
}
input_data <- readRDS("Data/ACORN_Indonesia/merged_metadata_acorn_indonesia2025-08-13.rds")
# target_variable = "d28_status"
