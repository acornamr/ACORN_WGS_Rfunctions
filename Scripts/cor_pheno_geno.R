# this script is used to create a function to create correlation map between phenotype and genotype of 
# AMR isolated from ACORN project 
# written and maintained by Tung Trinh 
# Jan 22nd 2025 
# version 1.0
# For more information, please contact to tungts@oucru.org 
####################################################################################################
cormap_amr<- function(genotype = "Path to the genotype data ",
                      phenotype = "Path to the phenotype data",
                      pathogen = "target pathogen, either e.col, kleb, etc.",
                      optional_text = NULL){
  ##################### required packages and built-in function ##############################
  require(stringr)
  require(dplyr)
  require(ggplot2)
  ##################### load and check the data #############################################
  ##### Genotype 
  
  ##### Phenotype w
  #################### Process the data and merge to each other ############################
  
  #################### plot the correlation map and output the file ########################
  
  
  
}
