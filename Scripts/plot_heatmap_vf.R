# this script is used to plot the heatmap of the virulence factors of a particular pathogen. 
# the input should be the input matrix with row name is the sample ID/ sequence ID. 
# you can sort the order of the sample by another variable
# written and maintained by Tung Trinh
# May 14th 2025 
# version 1.0 
# For more information, please contact to tungts@oucru.org
####################################################################################################
plot_heatmap_vf <- function(input = "input data", var_id = NULL, metadata = NULL,
                            sort_by = NULL , outdir = "path to the output directory", optional_text = NULL,
                            title = NULL){
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
  source("Scripts/get_colour_code_func.R")
  ##### get the color code database #####
  if(!is.null(metadata) && !is.null(sort_by)){
  # find the latest colour code database 
  list_file <- list.files(path = "Data",pattern = "colour_code_db", full.names = T)
  file <- file.info(list_file) %>% arrange(desc(mtime)) %>% slice(1) %>% rownames()
  # read the file 
  clr_db <- readRDS(file = file)
  }
  # get the colour code 
  if(!is.null(metadata) && !is.null(sort_by)){
    if(sort_by %in% unique(clr_db$variable))
      # get color annotation for the column
      col_anno_clr <- get_colour_code(metadata,target_variable = sort_by)
  }
  
  ##### read the virulence factor dictionary and database  #####
  # read the virulence factor dictionary
  vf <- read_excel("Data/VFs.xls")
  colnames(vf) <- vf[1,] %>% as.vector()
  vf <- vf[-1,]
  # read data from vfdb 
  vfdb_detail <- read.table("Data/vfdb.tab", sep = "\t")
  vfdb_detail <- vfdb_detail %>% mutate(VFCID = str_extract(V14,"VF[0-9]+"))
  vfdb_detail <- vfdb_detail %>% select(Gene = V6,VFCID) %>% distinct(.keep_all = T)
  vfdb_detail <- left_join(vfdb_detail,vf, by = c("VFCID" = "VFID"))
  vfdb_detail <- vfdb_detail %>% mutate(Gene = str_replace(Gene,"/","."))
  ##### create color annotation list 
  row_anno_color <- list(Group = c("Adherence" = "#774762FF","Effector delivery system" = "#1B3A54FF",
                         "Nutritional/Metabolic factor" = "#D6BB3BFF","Immune modulation" = "#755028FF",
                         "Exotoxin" = "#F2DD78FF","Exoenzyme" = "#205F4BFF","Biofilm" = "#913914FF",
                         "Regulation" = "#585854FF" ,
                         "Antimicrobial activity/Competitive advantage" = "#F0A430FF","Invasion" = "#768048FF"))
  ##### process the data #####
  ### check the input
  if(is.data.frame(input) == F){
    stop("Please use the correct format of the dataset")
  }else{
    print(paste("Start processing the dataset"))
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
  print(paste("input matrix dimension:",dim(input)) )
  # rename the var_id 
  input %<>% rename(var_id = 1)
  # fix the var_id 
  # input %<>% mutate(var_id = str_replace_all(var_id,"assemblies|[/]|[.short.fasta]|[.long.fasta]",""))
  # rename the sort_by and var_id 
  if(!is.null(metadata) && !is.null(sort_by)){
    metadata %<>% dplyr::rename(sort_by = matches(sort_by))
    metadata %<>% dplyr::rename(var_id = one_of(var_id))
    # select the selected variable only from the metadata 
    to_sort_dat <- metadata %>% dplyr::select(var_id,sort_by)
  }
  ### join the virulence database with the input 
  if(!is.null(metadata)&& !is.null(sort_by)){
  to_work <- left_join(input,to_sort_dat)
  }else{
    to_work <- input 
  }
  # fix the matrix to plot
  out <- to_work %>% select(-var_id,-sort_by,-NUM_FOUND)
  out <- out %>% mutate_all(~as.character(.)) %>% as.matrix(.)
  out[out == "."] <- F
  out[out != F ] <- T
  # assign the name 
  row.names(out) <- to_work$var_id
  ##### plot the figure #####
  # column annotation 
  target_vf <- colnames(out)
  vf_annotate <- vfdb_detail %>% filter(Gene %in% target_vf) %>% select(Gene,VFcategory) %>%
    arrange(VFcategory) %>% distinct(.keep_all = T)
  col_order_hm <- vf_annotate %>% pull(Gene)
  # 
  colA <- data.frame(Group = vf_annotate$VFcategory)
  
  # row annotation
  rowA <- data.frame(Variable = to_work$sort_by) %>% arrange(Variable)
  rownames(rowA) <- to_work %>% arrange(sort_by) %>% pull(var_id)
  
  # create a list of color for column annotation 
  col_anno_cl <- list(Variable = (setNames(col_anno_clr,unique(metadata$sort_by) %>% sort())))
  annotated_row <- rowAnnotation(df = rowA, col = col_anno_cl)
  annotated_col <- HeatmapAnnotation(Group = colA %>% pull(Group), col = row_anno_color)
  
  row_order_hm <-  to_work %>% arrange(sort_by) %>% pull(var_id)
  out <- out[row_order_hm,col_order_hm]
  
  p<- Heatmap(out,column_title = title,col = c("gray80","brown"),
          # row_dend_reorder = F,
          row_order = row_order_hm,
          column_order = col_order_hm,
          clustering_method_rows = "ward.D",
          column_labels = rep("",(n = ncol(out))),
          top_annotation = annotated_col,
          # clustering_distance_columns =  "euclidean",
          heatmap_legend_param = list(
            title = "Genotype", labels = c("No","Yes")
          )
  )
  
##### output the figure #####
  jpeg(filename = paste0(outdir,"/heatmap_vf",optional_text,Sys.Date(),".jpg"), units = "in", res = 300, height = 6,width = 10)
  draw(p + annotated_row,annotation_legend_side = "top")
  dev.off()
}