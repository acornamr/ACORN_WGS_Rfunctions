# this script is used to create an additional sample_list to run the leftover samples, in case you are
# not able to use -resume flag 
# Written and maintained by Tung Trinh 
# Aug 20th 2025
# version 1.0
# For more information, please contact to tungts@oucru.org 
####################################################################################################
create_rerun_sample_list <- function(original_samplelist = "path to original file",
                                     completed_samplelist = "path to completed list file",
                                     path_out = "path to output dir",
                                     optional_text = NULL){
  ###### required packages and built-in functions ######
  require(dplyr)
  require(stringr)
  "%ni%" <- Negate("%in%")
  ##### get the files ######
  # read the original list, output from create_sample_list function
  original <- read.csv(original_samplelist, header = T)
  # read the completed list 
  completed <- read.table(completed_samplelist)
  ##### procressing ######
  # fix the id 
  completed <- completed %>% mutate(fixed_id = str_remove_all(V1,".long|.fasta|.short"))
  # find the remained data 
  output <- original %>% filter(sample_id %ni% completed$fixed_id)
  ##### output the file #####
  write.csv(output, file = paste0(path_out,"rerun_sample_list_",optional_text,Sys.Date(),".csv"),
            quote = F, row.names = F)
}
# original_samplelist = "Data/ACORN_Indonesia/sample_list_acorn_indonesia2025-08-13.csv"
# completed_samplelist = "DAta/ACORN_Indonesia/sequenced.tab"
# path_out = "Data/ACORN_Indonesia/"
# optional_text = "acorn_indonesia"
