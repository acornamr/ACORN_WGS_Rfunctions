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
  require(aplot)
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
  
  ## check whether the genotype data has required varible contain genotype
  if(any(grepl(paste0(var_geno),colnames(geno_input))) == T){
    print("Found the target genotype varible in genotype data:")
    print(colnames(geno_input)[grep(paste0(var_geno),colnames(geno_input))])
  }else{
    stop("please input the data with genotype")
  }
  # change the variable name 
  geno_input %<>% rename(target_geno = matches(var_geno))
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
  mic_w <- mic %>% select(-(to_be_excluded)) %>% mutate_all(as.numeric)
  # fix the column name 
  colnames(mic_w)<- str_replace(colnames(mic_w),"_NM|_EM","")
  ## Disk diffusion 
  nd <- pheno %>% select(matches("_ND"),matches("_ED")) %>% mutate_all(as.numeric)
  # discard the empty columns 
  to_be_excluded <- which((nd %>% map(~sum(is.na(.))) %>% unlist) == nrow(pheno))
  nd_w <- nd %>% select(-(to_be_excluded)) %>% mutate_all(as.numeric)
  ## E test 
  etest <- pheno %>% select(matches("_NE$"),matches("_EE$"))
  # discard the empty columns 
  to_be_excluded <- which((etest %>% map(~sum(is.na(.))) %>% unlist) == nrow(pheno))
  etest_w <- etest %>% select(-(to_be_excluded))
  
  ### transform the nd data 
  nd_dat <- colnames(nd_w) %>% str_replace_all(.,"[_]","") %>% str_split(.,"ND") %>% unlist %>% matrix(., ncol = 2, byrow = T) %>% as_tibble()
  colnames(nd_dat) <- c("Antibiotic","Threshold")
  nd_dat %<>% mutate(Threshold = as.numeric(Threshold))
  # create a proxy data frame 
  db <- data.frame(id = pheno$target_id)
  # create a for loop to calculate 
  for(ab in nd_d$Antibiotic){
    # search the target column 
    db$newid <- (nd_dat %>% filter(Antibiotic == ab) %>% pull(Threshold) %>% as.vector(.))/(nd_w[,grep(paste0(ab),colnames(nd_w))]) 
    colnames(db)[grep("newid", colnames(db))] <- ab
  }
  ### create the list of antibiotic tested 
  target_ab <- c(colnames(mic_w),nd_dat$Antibiotic) %>% unique()
  ### create the genotype table
  ## select AMR genotype 
  genotype <- geno_input %>% pull(target_geno) %>%  str_split(.,",") %>% 
    map(~str_replace(.," ","")) %>% map(~unique(.))   %>% mtabulate
  # get the list of the gene 
  target_gene <- colnames(genotype) %>% str_replace(.," ","")
  # filter the amr data based on the genotype 
  target_amrdb <- amrdb %>% mutate(gene_name = str_replace(Gene,"[_][1-9]$","")) %>% 
    filter(gene_name %in% target_gene) %>% select(Drug,gene_name,hydrolysis_group,amr_gene_group) %>% 
    distinct(., .keep_all = T)
  
  
  
  for(ab in target_ab){
    tryCatch({
      ab_name <- whonet %>% filter(Code == ab) %>% pull(Antibiotic)
      # get the ab result 
      db_ab <- db[,grep(ab,colnames(db))]
      db_ab$seqid <- pheno$target_id
      MIC_ab <- mic_w[,grep(ab,colnames(mic))]
      MIC_ab$seqid <- mic$seqid
      colnames(MIC_ab)[1] <- "antibiotic"
      # search the code from WHONET 
      name_ab <- whonet %>% filter(Code == ab) %>% pull(Antibiotic)
      # search the list of the genotype 
      if(ab_name == "Penicillin-G"){
        target_gene <- target_amrdb %>% filter(Drug == "benzylpenicillin") %>% pull(amr_gene_group)
      }else{
        target_gene <- target_amrdb %>% filter(Drug == tolower(gsub("/Sulbactam","",name_ab))) %>% pull(amr_gene_group)
        
      }
      # from target genotype to target gene group 
      to_plot <- genotype %>% select(seqid,contains(target_gene)) 
      # plot 
      p_mid  <- to_plot %>% melt() %>% mutate(value = as.factor(value)) %>% ggplot(data = ., aes(y = seqid,x = variable))+
        geom_tile(aes(fill = value))+
        theme_minimal()+ theme(legend.position = "none",
                               axis.text.y=element_blank(),
                               axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
                               panel.grid.major.x = element_blank()
                               # panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                               # panel.grid.minor=element_blank(),plot.background=element_blank()
        )+
        scale_fill_manual(values = c("grey80","brown"), labels = c("No","Yes"),name = "")+
        xlab("AMR Gene")+ylab("")+ggtitle(paste(ab_name,": Phenotype vs Genotype"))
      
      colnames(db_ab)[1] <- "antibiotic"
      if(ncol(db_ab) > 1){
        p_right <- ggplot(db_ab, aes(y = seqid , x = antibiotic)) + theme_minimal()+ 
          theme(legend.position = "top",
                axis.text.y=element_blank(),
                axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
                panel.grid.major.x = element_blank()
                # panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                # panel.grid.minor=element_blank(),plot.background=element_blank()
          )+
          xlab("Threshold/Discs Diffusion Ratio")+ylab("")+
          geom_col(fill= "steelblue")
      }
      if(ncol(MIC_ab) > 1){
        p_left <- ggplot(data = MIC_ab)+
          geom_col(aes(y = seqid, x= antibiotic, fill = antibiotic))+
          scale_x_reverse()+
          theme_minimal()+
          theme(legend.position = "none",
                axis.text.y=element_blank(),
                axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
                panel.grid.major.x = element_blank()
                # panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                # panel.grid.minor=element_blank(),plot.background=element_blank()
          )+
          scale_fill_continuous(high = "#132B43", low = "#56B1F7")+xlab("MIC value")
      }
      jpeg(file = paste0(output_dir,"/",prefix,ab,optional_text,Sys.Date(),".jpg"),units = "in", res = 300,height = 15,width = 15)
      if(ncol(MIC_ab) > 1 && ncol(db) >1){
        p_whole <- p_mid %>% insert_right(p_right) %>% insert_left(p_left)
        plot_list(p_whole) %>% print()
        dev.off()
      }else{
        if(ncol(MIC_ab) == 1){
          p_whole <- p_mid %>% insert_right(p_right)
          jpeg(file = paste0(output_dir,"/",prefix,ab,optional_text,Sys.Date(),".jpg"),units = "in", res = 300,height = 15,width = 15)
          plot_list(p_whole) %>% print()
          dev.off()
        }
        if(ncol(db_ab) == 1){
          
          p_whole <- p_mid  %>% insert_left(p_left)
          jpeg(file = paste0(output_dir,"/",prefix,ab,optional_text,Sys.Date(),".jpg"),units = "in", res = 300,height = 15,width = 15)
          plot_list(p_whole) %>% print()
          dev.off() 
        }
      }
    },error = function(e){print(paste("ERROR: ",ab,"is not included in the amrdb"))
  }
  
}
