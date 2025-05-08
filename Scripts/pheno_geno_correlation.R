# this script is used to analysis the correlation between AMR gene and MIC profiles of each isolation 
# the input is staramr profile and the metadata which consists of MIC value
# writen and maintained by Tung Trinh
# ver 1.0 
# Apr 21st 2025
# For more information, please contact to tungts@oucru.org
####################################################################################################
##### required packages and built-in packages #####
library(readxl)
library(stringr)
library(dplyr)
library(purrr)
library(qdapTools)
library(ggplot2)
library(aplot)
library(reshape2)
# library(ggtree)
##### get the data and preprocess to make them mergable #####
## amr profile - acine 
staramr_acine <- read_excel("Data/assemblies/Hanoi/staramr_acine/results.xlsx")
# create id 
staramr_acine <- staramr_acine %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+"))
## amr profile - kleb 
# staramr_kleb <- read_excel("Data/assemblies/Hanoi/staramr_kleb_plus/results.xlsx")
# # # create id 
# staramr_kleb <- staramr_kleb %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+"))
# ## amr profile - ecoli
# staramr_ecol <- read_excel("Data/assemblies/Hanoi/staramr_ecol/staramr_ecol_results.xlsx")
# # create id 
# staramr_ecol <- staramr_ecol %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+"))
# ## amr profile - staph
# staramr_staph <- read_excel("Data/assemblies/Hanoi/staramr_staph/results.xlsx")
# # create id 
# staramr_staph <- staramr_staph %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+")) 
## metadata 
metadata <- read_excel("Data/ACORN2_VN001_WGS_Cumulative_V2_AKP.xlsx")
## phenotype data 
whole_dat  <- read_excel("Data/acorn_data_2025-01-22_04H35.xlsx", sheet = "acorn_dta")
## WHONET code 
whonet <- readRDS("Data/AMRdb/WHONET_antibiotic_codes2025-04-21.rds")
## AMR database 
amrdb <- readRDS("Data/AMRdb/amr_genotype_db_updated_2025-04-10.rds")
##### processing the data #####

############################## Acinobacter subsp.  #################################################
# get the acorn id for acinobacter 
target_id_acine <- metadata %>% filter(`SEQ-ID`%in% staramr_acine$seqid) %>% pull(specid...15)
# target phenotype columns 
target_var <- colnames(whole_dat)[grep("[A-Z]+", colnames(whole_dat))]

# filter the phenotype data based on the id 
acorn_phenotype <- whole_dat %>% filter(specid %in% target_id_acine) %>% select(specid,contains(target_var))
# combine the result of duplicated ID 
acorn_phenotype <- acorn_phenotype %>% 
  group_by(specid) %>% 
  summarise(across(everything(), ~ ifelse(any(complete.cases(.x)),
                                          first(.x[!is.na(.x)]),
                                          NA)))
acorn_phenotype <- acorn_phenotype %>% left_join(metadata %>% filter(`SEQ-ID`%in% staramr_acine$seqid) %>% 
                                select(specid...15,`SEQ-ID`),., by = c("specid...15" = "specid")) %>% 
  mutate(seqid = `SEQ-ID`)
# remove variables contains no information
to_be_excluded <- which((colSums(is.na(acorn_phenotype)) %>% unlist) == nrow(acorn_phenotype))
acorn_phenotype <- acorn_phenotype %>% select(-to_be_excluded)
# 
## work on ND data 
# fix the column name of STX 
colnames(acorn_phenotype) <- gsub("SXT_ND1_2","SXT_ND12",colnames(acorn_phenotype))
ND <- colnames(acorn_phenotype)[grep("ND",colnames(acorn_phenotype))] %>% str_split(.,"_") %>%  
  unlist %>% matrix(., ncol = 2, byrow = T) %>% as_tibble()
colnames(ND) <- c("Antibiotic","Threshold")
# fix the Discs Diffusion 
ND <- ND %>% mutate(Threshold = gsub("ND","",Threshold) %>% as.numeric(.))
## calculate the ratio between Threshold and ND 
# select data for ND 
ND_acine <- acorn_phenotype %>% select(matches("_ND")) %>% mutate_all(as.numeric)
# create a proxy data frame 
db <- data.frame(id = acorn_phenotype$seqid)
# create a for loop to calculate 
for(ab in ND$Antibiotic){
  # search the target column 
  db$newid <- (ND %>% filter(Antibiotic == ab) %>% pull(Threshold) %>% as.vector(.))/(ND_acine[,grep(paste0(ab),colnames(ND_acine))]) 
  colnames(db)[grep("newid", colnames(db))] <- ab
}
# MIC 
mic <- acorn_phenotype %>% select(seqid,matches("NM")) %>% mutate_all(~gsub("[<>=]","",.)) %>% mutate(across(.cols = 2:24,.fns = as.numeric))
## select AMR genotype 
genotype <-
  staramr_acine %>% pull(Genotype) %>%  str_split(.,",") %>% 
  map(~str_replace(.," ","")) %>% map(~unique(.))   %>% mtabulate()
# get the list of the gene 
target_gene <- colnames(genotype) %>% str_replace(.," ","")
# filter the amr data based on the genotype 
target_amrdb <- amrdb %>% mutate(gene_name = str_replace(Gene,"[_][1-9]$","")) %>% 
  filter(gene_name %in% target_gene) %>% select(Drug,gene_name,hydrolysis_group,amr_gene_group) %>% 
  distinct(., .keep_all = T)

##### create a for loop to plot it all #### 
# create the list of target antibiotic 
target_ab <- colnames(acorn_phenotype)[grep("^[A-Z]",colnames(acorn_phenotype))] %>% str_sub(.,1,3) %>% unique(.)
target_ab <- target_ab[-1]
target_ab <- target_ab[-grep("FEP|CTX|CAZ|GEN|IPM|LVX|TZP|SXT|AMC|ERY|MFX|NIT|OXA",target_ab)]

for(ab in target_ab){
  # read transformed genotype data 
  t_genotype <- read_excel("Data/acorn_acine_transformed_amr_2025-04-24.xlsx", sheet = 2)
  
  ab_name <- whonet %>% filter(Code == ab) %>% pull(Antibiotic)
  # get the ab result 
  db_ab <- db[,grep(ab,colnames(db))]
  db_ab$seqid <- acorn_phenotype$seqid
  MIC_ab <- mic[,grep(ab,colnames(mic))]
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
  to_plot <- t_genotype %>% select(seqid,contains(target_gene)) 
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
  jpeg(file = paste0("Figures/acine_pheno_vs_geno_",ab,".jpg"),units = "in", res = 300,height = 15,width = 15)
  if(ncol(MIC_ab) > 1 && ncol(db) >1){
  p_whole <- p_mid %>% insert_right(p_right) %>% insert_left(p_left)
  plot_list(p_whole) %>% print()
  dev.off()
  }else{
    if(ncol(MIC_ab) == 1){
      p_whole <- p_mid %>% insert_right(p_right)
      jpeg(file = paste0("Figures/acine_pheno_vs_geno_",ab,".jpg"),units = "in", res = 300,height = 15,width = 15)
      plot_list(p_whole) %>% print()
      dev.off()
    }
    if(ncol(db_ab) == 1){
      
      p_whole <- p_mid  %>% insert_left(p_left)
      jpeg(file = paste0("Figures/acine_pheno_vs_geno_",ab,".jpg"),units = "in", res = 300,height = 15,width = 15)
      plot_list(p_whole) %>% print()
      dev.off() 
    }
  }
}


############################## Kleb  #################################################
# get the acorn id for acinobacter 
target_id_kleb <- metadata %>% filter(`SEQ-ID`%in% staramr_kleb$seqid) %>% pull(specid...15)
# target phenotype columns 
target_var <- colnames(whole_dat)[grep("[A-Z]+", colnames(whole_dat))]
# filter the phenotype data based on the id 
acorn_phenotype <- whole_dat %>% filter(specid %in% target_id_kleb) %>% select(specid,contains(target_var))
# combine the result of duplicated ID 
acorn_phenotype <- acorn_phenotype %>% 
  group_by(specid) %>% 
  summarise(across(everything(), ~ ifelse(any(complete.cases(.x)),
                                          first(.x[!is.na(.x)]),
                                          NA)))
acorn_phenotype <- acorn_phenotype %>% left_join(metadata %>% filter(`SEQ-ID`%in% staramr_kleb$seqid) %>% 
                                                   select(specid...15,`SEQ-ID`),., by = c("specid...15" = "specid")) %>% 
  mutate(seqid = `SEQ-ID`)
# remove variables contains no information
to_be_excluded <- which((colSums(is.na(acorn_phenotype)) %>% unlist) == nrow(acorn_phenotype))
acorn_phenotype <- acorn_phenotype %>% select(-to_be_excluded)
# 
## work on ND data 
# fix the column name of STX 
colnames(acorn_phenotype) <- gsub("SXT_ND1_2","SXT_ND12",colnames(acorn_phenotype))
ND <- colnames(acorn_phenotype)[grep("ND",colnames(acorn_phenotype))] %>% str_split(.,"_") %>%  
  unlist %>% matrix(., ncol = 2, byrow = T) %>% as_tibble()
colnames(ND) <- c("Antibiotic","Threshold")
# fix the Discs Diffusion 
ND <- ND %>% mutate(Threshold = gsub("ND","",Threshold) %>% as.numeric(.))
## calculate the ratio between Threshold and ND 
# select data for ND 
ND_kleb <- acorn_phenotype %>% select(matches("_ND")) %>% mutate_all(as.numeric)
# create a proxy data frame 
db <- data.frame(id = acorn_phenotype$seqid)
# create a for loop to calculate 
for(ab in ND$Antibiotic){
  # search the target column 
  db$newid <- (ND %>% filter(Antibiotic == ab) %>% pull(Threshold) %>% as.vector(.))/(ND_kleb[,grep(paste0(ab),colnames(ND_kleb))]) 
  colnames(db)[grep("newid", colnames(db))] <- ab
}
# MIC 
mic <- acorn_phenotype %>% select(seqid,matches("NM")) %>% mutate_all(~gsub("[<>=]","",.)) %>% mutate(across(.cols = 2:24,.fns = as.numeric))
## select AMR genotype 
genotype <- staramr_kleb %>% pull(Genotype) %>%  str_split(.,",") %>%
  map(~str_replace(.," ","")) %>% map(~unique(.))   %>% mtabulate()
# get the list of the gene 
target_gene <- colnames(genotype) %>% str_replace(.," ","")
# filter the amr data based on the genotype 
target_amrdb <- amrdb %>% mutate(gene_name = str_replace(Gene,"[_][1-9]$","")) %>% 
  filter(gene_name %in% target_gene) %>% select(Drug,gene_name,hydrolysis_group,amr_gene_group) %>% 
  distinct(., .keep_all = T)

##### create a for loop to plot it all #### 
# create the list of target antibiotic 
target_ab <- colnames(acorn_phenotype)[grep("^[A-Z]",colnames(acorn_phenotype))] %>% str_sub(.,1,3) %>% unique(.)
target_ab <- target_ab[-1]
target_ab <- target_ab[-grep("ATM|FEP|CTX|CAZ|IPM|LVX|TZP|SXT|COL|MFX|NIT|OXA|PEN|VAN|ETP|CLI",target_ab)]

for(ab in target_ab){
  # read transformed genotype data 
  t_genotype <- read_excel("Data/acorn_kleb_transformed_amr_2025-04-21.xlsx", sheet = 2)
  
  ab_name <- whonet %>% filter(Code == ab) %>% pull(Antibiotic)
  # get the ab result 
  db_ab <- db[,grep(ab,colnames(db))]
  db_ab$seqid <- acorn_phenotype$seqid
  MIC_ab <- mic[,grep(ab,colnames(mic))]
  MIC_ab$seqid <- mic$seqid
  colnames(MIC_ab)[1] <- "antibiotic"
  # search the code from WHONET 
  name_ab <- whonet %>% filter(Code == ab) %>% pull(Antibiotic) %>% gsub(" ","-",.)
  # search the list of the genotype 
  if(ab_name == "Penicillin-G"){
    target_gene <- target_amrdb %>% filter(Drug == "benzylpenicillin") %>% pull(amr_gene_group)
  }else{
  target_gene <- target_amrdb %>% filter(Drug == tolower(gsub("/Sulbactam","",name_ab))) %>% pull(amr_gene_group)
  
  }
  target_gene <- na.exclude(target_gene) %>% as.vector()
  # from target genotype to target gene group 
  to_plot <- t_genotype %>% select(seqid,contains(target_gene)) 
  # plot 
  p_mid  <- to_plot %>% melt() %>% mutate(value = as.factor(value)) %>% ggplot(data = ., aes(y = seqid,x = variable))+
    geom_tile(aes(fill = value))+
    theme_classic()+ theme(legend.position = "none",
                           axis.text.y=element_blank(),
                           axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
                           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                           panel.grid.minor=element_blank(),plot.background=element_blank())+
    scale_fill_manual(values = c("grey80","brown"), labels = c("No","Yes"),name = "")+
    xlab("AMR Gene")+ylab("")+ggtitle(paste(ab_name,": Phenotype vs Genotype"))
  
  colnames(db_ab)[1] <- "antibiotic"
  if(ncol(db_ab) > 1){
    p_right <- ggplot(db_ab, aes(y = seqid , x = antibiotic)) + theme_classic()+ 
      theme(legend.position = "top",
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank())+
      xlab("Threshold/Discs Diffusion Ratio")+ylab("")+
      geom_col(fill= "steelblue")
  }
  if(ncol(MIC_ab) > 1){p_left <- ggplot(data = MIC_ab)+
    geom_col(aes(y = seqid, x= antibiotic, fill = antibiotic))+
    scale_x_reverse()+
    theme_classic()+
    theme(legend.position = "none",
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())+
    scale_fill_continuous(high = "#132B43", low = "#56B1F7")+xlab("MIC value")
  }
  jpeg(file = paste0("Figures/kleb_pheno_vs_geno_",ab,".jpg"),units = "in", res = 300,height = 15,width = 15)
  if(ncol(MIC_ab) > 1 && ncol(db) >1){
    p_whole <- p_mid %>% insert_right(p_right) %>% insert_left(p_left)
    plot_list(p_whole) %>% print()
    dev.off()
  }else{
    if(ncol(MIC_ab) == 1){
      p_whole <- p_mid %>% insert_right(p_right)
      jpeg(file = paste0("Figures/kleb_pheno_vs_geno_",ab,".jpg"),units = "in", res = 300,height = 15,width = 15)
      plot_list(p_whole) %>% print()
      dev.off()
    }
    if(ncol(db_ab) == 1){
      
      p_whole <- p_mid  %>% insert_left(p_left)
      jpeg(file = paste0("Figures/kleb_pheno_vs_geno_",ab,".jpg"),units = "in", res = 300,height = 15,width = 15)
      plot_list(p_whole) %>% print()
      dev.off() 
    }
  }
}

############################## E coli  #################################################
# get the acorn id for acinobacter 
target_id_ecol <- metadata %>% filter(`SEQ-ID`%in% staramr_ecol$seqid) %>% pull(specid...15)
# target phenotype columns 
target_var <- colnames(whole_dat)[grep("[A-Z]+", colnames(whole_dat))]
# filter the phenotype data based on the id 
acorn_phenotype <- whole_dat %>% filter(specid %in% target_id_ecol) %>% select(specid,contains(target_var))
# combine the result of duplicated ID 
acorn_phenotype <- acorn_phenotype %>% 
  group_by(specid) %>% 
  summarise(across(everything(), ~ ifelse(any(complete.cases(.x)),
                                          first(.x[!is.na(.x)]),
                                          NA)))
acorn_phenotype <- acorn_phenotype %>% left_join(metadata %>% filter(`SEQ-ID`%in% staramr_ecol$seqid) %>% 
                                                   select(specid...15,`SEQ-ID`),., by = c("specid...15" = "specid")) %>% 
  mutate(seqid = `SEQ-ID`)
# remove variables contains no information
to_be_excluded <- which((colSums(is.na(acorn_phenotype)) %>% unlist) == nrow(acorn_phenotype))
acorn_phenotype <- acorn_phenotype %>% select(-to_be_excluded)
# 
## work on ND data 
# fix the column name of STX 
colnames(acorn_phenotype) <- gsub("SXT_ND1_2","SXT_ND12",colnames(acorn_phenotype))
ND <- colnames(acorn_phenotype)[grep("ND",colnames(acorn_phenotype))] %>% str_split(.,"_") %>%  
  unlist %>% matrix(., ncol = 2, byrow = T) %>% as_tibble()
colnames(ND) <- c("Antibiotic","Threshold")
# fix the Discs Diffusion 
ND <- ND %>% mutate(Threshold = gsub("ND","",Threshold) %>% as.numeric(.))
## calculate the ratio between Threshold and ND 
# select data for ND 
ND_ecol <- acorn_phenotype %>% select(matches("_ND")) %>% mutate_all(as.numeric)
# create a proxy data frame 
db <- data.frame(id = acorn_phenotype$seqid)
# create a for loop to calculate 
for(ab in ND$Antibiotic){
  # search the target column 
  db$newid <- (ND %>% filter(Antibiotic == ab) %>% pull(Threshold) %>% as.vector(.))/(ND_ecol[,grep(paste0(ab),colnames(ND_ecol))]) 
  colnames(db)[grep("newid", colnames(db))] <- ab
}
# MIC 
mic <- acorn_phenotype %>% select(seqid,matches("NM")) %>% mutate_all(~gsub("[<>=]","",.)) %>% mutate(across(.cols = 2:24,.fns = as.numeric))
## select AMR genotype 
genotype <- staramr_ecol %>% pull(Genotype) %>%  str_split(.,",") %>%
  map(~str_replace(.," ","")) %>% map(~unique(.))   %>% mtabulate()
# get the list of the gene 
target_gene <- colnames(genotype) %>% str_replace(.," ","")
# filter the amr data based on the genotype 
target_amrdb <- amrdb %>% mutate(gene_name = str_replace(Gene,"[_][1-9]$","")) %>% 
  filter(gene_name %in% target_gene) %>% select(Drug,gene_name,hydrolysis_group,amr_gene_group) %>% 
  distinct(., .keep_all = T)

##### create a for loop to plot it all #### 
# create the list of target antibiotic 
target_ab <- colnames(acorn_phenotype)[grep("^[A-Z]",colnames(acorn_phenotype))] %>% str_sub(.,1,3) %>% unique(.)
target_ab <- target_ab[-1]
target_ab <- target_ab[-grep("CTX|LVX|FEP|CAZ|IPM|MFX|NIT|OXA|PEN|TZP|SXT|VAN|ETP|CLI",target_ab)]
for(ab in target_ab){
  # read transformed genotype data 
  t_genotype <- read_excel("Data/acorn_ecoli_transformed_amr_2025-04-16.xlsx", sheet = 2)
  
  ab_name <- whonet %>% filter(Code == ab) %>% pull(Antibiotic)
  # get the ab result 
  db_ab <- db[,grep(ab,colnames(db))]
  db_ab$seqid <- acorn_phenotype$seqid
  MIC_ab <- mic[,grep(ab,colnames(mic))]
  MIC_ab$seqid <- mic$seqid
  colnames(MIC_ab)[1] <- "antibiotic"
  # search the code from WHONET 
  name_ab <- whonet %>% filter(Code == ab) %>% pull(Antibiotic) %>% gsub(" ","-",.)
  # search the list of the genotype 
  if(ab_name == "Penicillin-G"){
    target_gene <- target_amrdb %>% filter(Drug == "benzylpenicillin") %>% pull(amr_gene_group)
  }else{
    target_gene <- target_amrdb %>% filter(Drug == tolower(gsub("/Sulbactam","",name_ab))) %>% pull(amr_gene_group)
    
  }
  target_gene <- na.exclude(target_gene) %>% as.vector()
  # from target genotype to target gene group 
  to_plot <- t_genotype %>% select(seqid,contains(target_gene)) 
  # plot 
  p_mid  <- to_plot %>% melt() %>% mutate(value = as.factor(value)) %>% ggplot(data = ., aes(y = seqid,x = variable))+
    geom_tile(aes(fill = value))+
    theme_classic()+ theme(legend.position = "none",
                           axis.text.y=element_blank(),
                           axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
                           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                           panel.grid.minor=element_blank(),plot.background=element_blank())+
    scale_fill_manual(values = c("grey80","brown"), labels = c("No","Yes"),name = "")+
    xlab("AMR Gene")+ylab("")+ggtitle(paste(ab_name,": Phenotype vs Genotype"))
  
  colnames(db_ab)[1] <- "antibiotic"
  if(ncol(db_ab) > 1){
    p_right <- ggplot(db_ab, aes(y = seqid , x = antibiotic)) + theme_classic()+ 
      theme(legend.position = "top",
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank())+
      xlab("Threshold/Discs Diffusion Ratio")+ylab("")+
      geom_col(fill= "steelblue")
  }
  if(ncol(MIC_ab) > 1){p_left <- ggplot(data = MIC_ab)+
    geom_col(aes(y = seqid, x= antibiotic, fill = antibiotic))+
    scale_x_reverse()+
    theme_classic()+
    theme(legend.position = "none",
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())+
    scale_fill_continuous(high = "#132B43", low = "#56B1F7")+xlab("MIC value")
  }
  jpeg(file = paste0("Figures/ecol_pheno_vs_geno_",ab,".jpg"),units = "in", res = 300,height = 15,width = 15)
  if(ncol(MIC_ab) > 1 && ncol(db) >1){
    p_whole <- p_mid %>% insert_right(p_right) %>% insert_left(p_left)
    plot_list(p_whole) %>% print()
    dev.off()
  }else{
    if(ncol(MIC_ab) == 1){
      p_whole <- p_mid %>% insert_right(p_right)
      jpeg(file = paste0("Figures/ecol_pheno_vs_geno_",ab,".jpg"),units = "in", res = 300,height = 15,width = 15)
      plot_list(p_whole) %>% print()
      dev.off()
    }
    if(ncol(db_ab) == 1){
      
      p_whole <- p_mid  %>% insert_left(p_left)
      jpeg(file = paste0("Figures/ecol_pheno_vs_geno_",ab,".jpg"),units = "in", res = 300,height = 15,width = 15)
      plot_list(p_whole) %>% print()
      dev.off() 
    }
  }
}





