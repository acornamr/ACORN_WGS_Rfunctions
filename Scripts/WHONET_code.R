# This script is used to convert WHONET antibiotic code into usable format 
# Written and maintained by Tung Trinh 
# April 21st 2025
# Version 1.0 
# For more information, please contact to tungts@oucru.org 
####################################################################################################
##### required packges and built-in functions #####
library(dplyr)
library(stringr)
library(pdftools)
library(purrr)
##### process the file #####
pdf <- pdf_text(pdf = "Data/AMRdb/WHONET_antibiotic_code.pdf")
pdf <- paste0(pdf, collapse = " ")
pdf <- gsub('\n', "  ", pdf)
extract <- strsplit(pdf, "\\s{2,}") %>% unlist() 
extract <- extract[-c(1:8)]
out <- extract %>% matrix(., ncol = 2, byrow = T)
out <- as_tibble(out)
colnames(out) <- c("Code","Antibiotic")
# output 
saveRDS(out,paste0("Data/AMRdb/WHONET_antibiotic_codes",Sys.Date(),".rds"))
