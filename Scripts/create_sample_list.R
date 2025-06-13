# this script is used to create sample sheet from the sample list provided 
# written and maintained by Tung Trinh 
# May 15th 2025
# Version 2.0
# June 12th 2025
# June 11th 2025
# For more information, please contact to tungts@oucru.org
####################################################################################################
create_sample_list <- function(path = "path to the sample list", optional_text = NULL, 
                               path_fastq = "path to the fastq files's directory",
                               read_short = TRUE, read_long = FALSE,genome_size = "predicted genome size",
                               path_out = "path to the output directory"){
  ##### required packages and built-in function ##### 
  require(dplyr)
  require(stringr) 
  require(data.table)
  ##### load the sample list #####
  input <- read.table(path, sep = "\t")
  ##### process the input #####
  # search whether we have different ID 
  # get the list file 
  list_files <- input$V1
  if(read_short == T){
  # create unique sequencing ID if short_read is true
    short_id <- input$V1 %>% as.vector() %>% str_replace(.,"_001","") %>% str_replace(.,"_R[1-2].fastq.gz","") %>% unique()
  }
  if(any(grepl("[A-Z]|[a-z]",short_id)) == F){
    short_id_fix <- paste0("sample_",short_id)
  }else{
    short_id_fix <- short_id
  }
  if(read_long == T){
    long_id <- input$V1
  }
  ### create a NULL data frame 
  output <- data.frame(sample_id = character(),
                       short_reads1 = character(),
                       short_reads2 = character(),
                       long_reads = character(),
                       genome_size = character())
  # assign the value 
  if(read_short == T){
    # output <- output %>%  add_row(sample_id = short_id_fix)
    for(id in short_id_fix){
      read1 <- paste0(path_fastq,"/",list_files[grep(paste0(id,".*R1"),list_files)])
      read2 <- paste0(path_fastq,"/",list_files[grep(paste0(id,".*R2"),list_files)])
      output <- output %>% add_row(sample_id = id,short_reads1 = read1,short_reads2 = read2)
    }
    output$long_reads <- NULL
    output$genome_size <- genome_size 
  }
  if(read_long == T){
    output <- output %>%  add_row( sample_id = long_id)
    output$long_reads <- paste0(path_fastq,"/",long_id)
    output$genome_size <- genome_size 
  }
  ##### output the file #####
  write.csv(output, file = paste0(path_out,"sample_list_",optional_text,Sys.Date(),".csv"),
            quote = F, row.names = F)
  
}