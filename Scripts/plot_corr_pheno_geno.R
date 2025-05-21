# this script is to create a function to plot the correlation between the phenotype and genotype 
# of antibiotic resistance among whole genome sequencing isolates. 
# the input should be transformed genotype data and the phenotype data 
# Written and maintained by Tung Trinh
# May 7th 2025
# Ver 1.0
# For more information, please contact to tungts@oucru.org 
####################################################################################################
plot_corr_pheno_geno <- function(pheno_input = "path to the metadata file (extracted from ACORN)",
                                 pheno_id = "variable contains the id within phenotype data",
                                 geno_id = "variable contains the id within genotype data",
                                 var_geno = "variable contains the genotype variable",
                                 geno_input = "data contains the genotype data",
                                 optional_text = NULL,
                                 output_dir = "path to the output folder",
                                 prefix = NULL){
  ##### required packages and built-in functions #####
  require(dplyr)
  require(magrittr)
  require(stringr)
  require(purrr)
  require(qdapTools)
  require(reshape2)
  require(ggplot2)
  require(alot)
  ##### AMR and WHONET database #####
  ## WHONET code 
  # list the available data and search for the latest one 
  list_whonet <- list.files(path = "Data/AMRdb/",pattern = "WHONET_antibiotic_codes",full.names = T)
  file_whonet <- file.info(list_whonet) %>% arrange(desc(mtime)) %>% slice(1) %>% rownames(.)
  whonet <- readRDS(file_whonet)
  ## AMR database 
  list_amrdb <- list.files(path = "Data/AMRdb/", pattern = "amr_genotype_db",full.names = T)
  file_amr <- file.info(list_amrdb) %>% arrange(desc(mtime)) %>% slice(1) %>% rownames(.)
  amrdb <- readRDS(file_amr)
  ##### check and process the data #####
  ## check whether the phenodata has required variables
  if(any(grepl("_ND|_NM|_NE|_EE|_ED|_EM",colnames(pheno_input))) == T){
    print("Found the target phenotype variables:")
    print(colnames(pheno_input)[grep("_ND|_NM|_NE|_EE|_ED|_EM",colnames(pheno_input))])
  }else{
    stop("please input the data with phenotype variables")
  }
  ## check whether the genotype data has requried variables
  if(any(grepl(paste0(var_geno),colnames(geno_input))) == T){
    print("Found the target Genotype variables:")
    print(colnames(geno_input)[grep(paste0(var_geno),colnames(geno_input))])
  }else{
    stop("please input the data with genotype variable")
  }
  # change the variable name of genotype 
  geno_input %<>% rename(target_geno = matches(var_geno))
  ## check whether the phenodata has required id variable 
  if(any(grepl(paste0(pheno_id),colnames(pheno_input))) == T){
    print("Found the target id varible:")
    print(colnames(pheno_input)[grep(paste0(pheno_id),colnames(pheno_input))])
  }else{
    stop("please input the data with id variable")
  }
  # change the variable name 
  pheno_input %<>% rename(target_id = matches(pheno_id))
  ## check whether the phenodata has required id variable 
  if(any(grepl(paste0(pheno_id),colnames(pheno_input))) == T){
    print("Found the target id varible in phenotype data:")
    print(colnames(pheno_input)[grep(paste0(pheno_id),colnames(pheno_input))])
  }else{
    stop("please input the pheno data with id variable")
  }
  # change the variable name 
  pheno_input %<>% rename(target_id = matches(pheno_id))
  ## check whether the genotype data has required id variable 
  if(any(grepl(paste0(geno_id),colnames(geno_input))) == T){
    print("Found the target id varible in genotype data:")
    print(colnames(geno_input)[grep(paste0(geno_id),colnames(geno_input))])
  }else{
    stop("please input the data with id variable")
  }
  # change the variable name 
  geno_input %<>% rename(target_id = matches(geno_id))
  # fix the id variable of the genotype data 
  geno_input <- geno_input %>% mutate(target_id = str_replace(target_id,".short|.fasta|.fa",""))
  ## check the matching id between two datasets 
  if(any(pheno_input$target_id %in% geno_input$target_id) == T){
    print("found matched sample between two datasets:")
    print(pheno_input$target_id[which(pheno_input$target_id %in% geno_input$target_id)])
  }else{
    stop("There is no matching samples between phenotype and genotype data")
  }
  ##### process the phenotype data #####
  # Firstly, filter the data based on th id of genotype 
  pheno <- pheno_input %>% filter(target_id %in% geno_input$target_id)
  # assign the id as the row name 
  rownames(pheno) <- pheno$target_id
  ## seperate the data into different datasets
  ## MIC 
  mic <- pheno %>% select(matches("_NM"),matches("_EM"))
  # discard the empty columns 
  to_be_excluded <- which((mic %>% map(~sum(is.na(.))) %>% unlist) == nrow(pheno))
  mic_w <- mic %>% select(-(to_be_excluded))
  ## Disk diffusion 
  nd <- pheno %>% select(matches("_ND"),matches("_ED"))
  # discard the empty columns 
  to_be_excluded <- which((nd %>% map(~sum(is.na(.))) %>% unlist) == nrow(pheno))
  nd_w <- nd %>% select(-(to_be_excluded))
  ## E test 
  etest <- pheno %>% select(matches("_NE$"),matches("_EE$"))
  # discard the empty columns 
  to_be_excluded <- which((etest %>% map(~sum(is.na(.))) %>% unlist) == nrow(pheno))
  etest_w <- etest %>% select(-(to_be_excluded))
  
}
