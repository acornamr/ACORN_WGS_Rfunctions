# this script is used to plot the heatmap of the virulence factors of a particular pathogen. 
# the input should be the input matrix with row name is the sample ID/ sequence ID. 
# you can sort the order of the sample by another variable
# written and maintained by Tung Trinh
# May 14th 2025 
# version 1.0 
# For more information, please contact to tungts@oucru.org
####################################################################################################
plot_heatmap_vf <- function(input = "input data", var_id = "variable contains the sample ID", metadata = NULL,
                            sort_by = NULL , outdir = "path to the output directory", optional_text = NULL){
  ###### required packages and built-in functions #####
  require(dplyr)
  require(readxl)
  require(magrittr)
  require(stringr)
  require(tidyr)
  require(ggplot2)
  require(gridExtra)
  require(patchwork)
  require(data.table)
  require(ComplexHeatmap)
  ##### read the virulence factor dictionary and database  #####
  # read the virulence factor dictionary
  vf <- read_excel("Data/VFs.xls")
  # create the column name 
  colnames(vf) <- vf[1,] %>% as.vector()
  # discard the first row 
  vf <- vf[-1,]
  ##### process the data #####
  ### check the input
  if(is.data.frame(input) == F){
    stop("Please use the correct format of the dataset")
  }else{
    print(paste("Start processing the dataset"))
  }
  # check whether the var_id variable available in the dataset or not
  print("Check the target variable available in the dataset or not")
  if((var_id %in% colnames(input)) == F){
    stop(paste("Could not find",var_id,"in the dataset"))
  }else{
    print(paste("The target variable",var_id,"found in the dataset"))
  }
  # if sort_by is supplied, check whether sort_by variable is supplied or not 
  if(!is.null(sort_by)){
    if((sort_by %in% colnames(metadata)) == F){
      stop(paste("Could not find",sort_by,"in the dataset"))
    }else{
      print(paste("The target variable",sort_by,"found in the dataset"))
    }
  }
  # print the summary of the input dataset 
  print(paste("input matrix dimension:",dim(dat)) )
  # rename the var_id 
  input %<>% rename(var_id = matches(var_id))
  # rename the sort_by and var_id 
  if(!is.null(metadata) && !is.null(sort_by)){
    metadata %<>% rename(sort_by = matches(sort_by))
    metadata %<>% rename(var_id = matches(var_id))
    # select the selected variable only from the metadata 
    to_sort_dat <- metadata %>% select(var_id,sort_by)
  }
  ### join the virulence database with the input 
  
  ##### plot the figure #####
  
  ##### output the figure #####
  
}