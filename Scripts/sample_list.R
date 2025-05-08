# this script is used to create sample list for snippy input 
# the input is the list of sample from remote server  
# written and maintained by Tung Trinh 
# ver 1.0
# March 3rd 2025
# for more information, please contact to tungts@oucru.org 
####################################################################################################
##### required packages and built-in functions #####
library(dplyr)
library(stringr)
library(data.table)
##### create sample list #####
### E.coli 
# read the input sample list 
ecol_sample_list <- read.table("Data/sample_list/sample_ecol.tab")
# create the id for each sample 
ecol_sample_list <- ecol_sample_list %>% mutate(id = str_extract(V1,"A[0-9]-[0-9]+"))
setDT(ecol_sample_list)[is.na(id), id:="A3-12"]
# create the path 
ecol_sample_list <- ecol_sample_list %>% mutate(path = paste0("assemblies/",V1))
# create sample list 
write.table(ecol_sample_list %>% select(id, path), "Data/sample_list/ecol_input_list.tab", 
            quote = F, col.names = F, sep = "\t", row.names = F)

### Acine.bau 
# read the input sample list 
acine_sample_list <- read.table("Data/sample_list/sample_acine.tab")
# create the id for each sample 
acine_sample_list <- acine_sample_list %>% mutate(id = str_extract(V1,"A[0-9]-[0-9]+"))
# create the path 
acine_sample_list <- acine_sample_list %>% mutate(path = paste0("assemblies/",V1))
# create sample list 
write.table(acine_sample_list %>% select(id, path), "Data/sample_list/acine_input_list.tab", 
            quote = F, col.names = F, sep = "\t", row.names = F)

### Kleb.pneu
# read the input sample list 
kleb_sample_list <- read.table("Data/sample_list/sample_kleb.tab")
# create the id for each sample 
kleb_sample_list <- kleb_sample_list %>% mutate(id = str_extract(V1,"A[0-9]-[0-9]+"))
# create the path 
kleb_sample_list <- kleb_sample_list %>% mutate(path = paste0("assemblies/",V1))
# create sample list 
write.table(kleb_sample_list %>% select(id, path), "Data/sample_list/kleb_input_list.tab", 
            quote = F, col.names = F, sep = "\t", row.names = F)

### Staph
# read the input sample list 
staph_sample_list <- read.table("Data/sample_list/sample_staph.tab")
# create the id for each sample 
staph_sample_list <- staph_sample_list %>% mutate(id = str_extract(V1,"A[0-9]-[0-9]+"))
# create the path 
staph_sample_list <- staph_sample_list %>% mutate(path = paste0("assemblies/",V1))
# create sample list 
write.table(staph_sample_list %>% select(id, path), "Data/sample_list/staph_input_list.tab", 
            quote = F, col.names = F, sep = "\t", row.names = F)

