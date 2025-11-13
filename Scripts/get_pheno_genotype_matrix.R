# this script is be used to create the matrix for the genotype and phenotype correlation visualisation
# designed for ACORN project
# Nov 13 2025 
# Ver 1.0
# Written and maintained by Tung Trinh
# for more information, please contact to tungts@oucru.org
####################################################################################################
get_pheno_genotype_matrix <- function(input_genotype = "", input_phenotype = "",
                                      target_species = "", target_antibiotic = ""){
  ##### required packages and built-in functions ######
  require(dplyr)
  require(stringr)
  ##### get the database #####
  ab_class <- readRDS("Data/AMRdb/ab_class_db_updated_2025-11-13.rds")
  ## read the latest database
  # list the database 
  list_amr_genotype_db <- list.files("Data/AMRdb/", pattern = "amr_genotype",full.names = T)
  # find the latest file 
  file_amr_genotype_db <- file.info(list_amr_genotype_db) %>% arrange(desc(mtime)) %>% slice(1) %>% rownames()
  # read the database 
  amrdb <- readRDS(file_amr_genotype_db)
  ##### processing the data #####
  # first, filter the ab_class database 
  target_ab <- ab_class %>% filter(antibiotic_class == target_antibiotic & species == target_species) %>% pull(antibiotics)
  # get the target variables from the phenotype data
  ind_phe <- which(colnames(input_phenotype) %in% target_ab)
  out_pheno <- input_phenotype %>% select(ind_phe) 
  # get the related genes from genotype data
  amr_gene <- amrdb %>% filter(grepl(tolower(target_antibiotic),tolower(aware_subclass))) %>% pull(amr_gene_group) %>% unique()
  # find related gene from genotype data
  # get the indicator for column 
  ind_gen <- which(colnames(input_genotype) %in% amr_gene)
  # select data 
  out_geno <- input_genotype %>% select(seqid,ind_gen)
  ##### output ##### 
  # combine the two output 
  out <- list(out_geno,out_pheno)
  names(out) <- c("Genotype","Phenotype")
  # output 
  return(out)
}
