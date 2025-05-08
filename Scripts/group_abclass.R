# this script is used to check group up the antibiotic resistant genes by their hydrolysis class.
# Written and maintained by Tung Trinh
# March 27th 2025
# ver 1.1
# Updated on April 9th 2025
# for more information, please contact to tungts@oucru.org
####################################################################################################
###### required packages and built-in functions ######
library(dplyr)
library(stringr)
library(data.table)
###### read the data ##### 
amrdb <- readRDS("Data/AMRdb/amr_genotype_db_2025-02-14.rds")
####### modify the data #####
# fix the name of the gene
amrdb <- amrdb %>% mutate(gene_name = str_replace(Gene,"_1",""))
### create hydrolysis classes
amrdb$hydrolysis_group <- NA_character_
### create amr gene group 
amrdb$amr_gene_group <- NA_character_
### create group ESBLs
# first, create a vector of the listed gene
ESBLs <-  c("blaTEM-1A","blaTEM-1B","blaTEM-1C","blaTEM-1D","blaTEM-234","blaTEM-40","blaTEM-70",
            "blaSHV-100","blaSHV-108","blaSHV-148","blaSHV-182","blaSHV-187","blaSHV-190","blaSHV-26","blaSHV-42","blaSHV-59","blaSHV-61","blaSHV-67","blaSHV-80","blaSHV-81","blaSHV-89",
            "blaCTX-M-14","blaCTX-M-15","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-82",
            "blaOXA-1","blaOXA-10","blaLAP-2","blaSFO-1")
# then we will mark it as ESBLs
setDT(amrdb)[gene_name %in% ESBLs, hydrolysis_group := "ESBLs"]
## blaTEM group
blaTEM <- c("blaTEM-1A","blaTEM-1B","blaTEM-1C","blaTEM-1D","blaTEM-234","blaTEM-40","blaTEM-70")
setDT(amrdb)[gene_name %in% blaTEM, amr_gene_group := "blaTEM"]
## blaSHV
blaSHV <- c("blaSHV-100","blaSHV-108","blaSHV-148","blaSHV-182","blaSHV-187","blaSHV-190","blaSHV-26","blaSHV-42","blaSHV-59","blaSHV-61","blaSHV-67","blaSHV-80","blaSHV-81","blaSHV-89")
setDT(amrdb)[gene_name %in% blaSHV, amr_gene_group := "blaSHV"]
## blaCTX-M
blaCTX <- c("blaCTX-M-14","blaCTX-M-15","blaCTX-M-27","blaCTX-M-3","blaCTX-M-55","blaCTX-M-82")
setDT(amrdb)[gene_name %in% blaCTX, amr_gene_group := "blaCTX-M"]
## ESBL OXA 
esbl_oxa <- c("blaOXA-1","blaOXA-10")
setDT(amrdb)[gene_name %in% esbl_oxa, amr_gene_group := "OXA results ESBL"]
## Other ESBL 
esbl_other <- c("blaLAP-2","blaSFO-1")
setDT(amrdb)[gene_name %in% esbl_other, amr_gene_group := "Other genes result ESBL"]
### create group AmpC
# first, create a vector of listed gene 
AmpC <- c("blaCMY-145","blaCMY-2","blaCMY-42","blaDHA-1","blaADC-25")
# then we will mark it as AmpC
setDT(amrdb)[gene_name %in% AmpC, hydrolysis_group := "AmpC"]
## blaCMY
blaCMY <- c("blaCMY-145","blaCMY-2","blaCMY-42")
setDT(amrdb)[gene_name %in% blaCMY, amr_gene_group := "blaCMY"]
## blaDHA
blaDHA <- c("blaDHA-1")
setDT(amrdb)[gene_name %in% blaDHA, amr_gene_group := "blaDHA"]
## blaADC 
blaADC <- c("blaADC-25")
setDT(amrdb)[gene_name %in% blaADC, amr_gene_group := "blaADC"]
### create group Carbapenemase
# first, create a vector of listed gene 
carb <- c("blaOXA-48","blaOXA-181","blaOXA-484",
          "blaCARB-5","blaNDM-1","blaNDM-4","blaNDM-5","blaIMP-54","blaKPC-2")
# then we will mark it as Carbapenemase
setDT(amrdb)[gene_name %in% carb, hydrolysis_group := "Carbapenemase Resistance"]
## OXA-48
bla_oxa48 <- c("blaOXA-48","blaOXA-181","blaOXA-484")
setDT(amrdb)[gene_name %in% bla_oxa48, amr_gene_group := "blaOXA-48"]
## blaCARB
blaCARB <- c("blaCARB-5")
setDT(amrdb)[gene_name %in% blaCARB, amr_gene_group := "blaCARB"]
## blaNDM
blaNDM <- c("blaNDM-1","blaNDM-4","blaNDM-5")
setDT(amrdb)[gene_name %in% blaNDM, amr_gene_group := "blaNDM"]
## blaIMP
blaIMP <- c("blaIMP-54")
setDT(amrdb)[gene_name %in% blaIMP, amr_gene_group := "blaIMP"]
## blaKPC 
blaKPC <- c("blaKPC-2")
setDT(amrdb)[gene_name %in% blaKPC, amr_gene_group := "blaKPC"]
### create group CHDL
# first, create a vector of listed genes 
CHDL <- c("blaOXA-23","blaOXA-24","blaOXA-510","blaOXA-217","blaOXA-69","blaOXA-91","blaOXA-98","blaOXA-66","blaOXA-120",
          "blaOXA-58","blaOXA-143","blaOXA-235","blaOXA-417")
# then we mark them as CHDL
setDT(amrdb)[gene_name %in% CHDL, hydrolysis_group := "CHDLs"]
## oxa-23
oxa23 <- c("blaOXA-23")
setDT(amrdb)[gene_name %in% oxa23, amr_gene_group := "blaOXA-23"]
## oxa-24
oxa24 <- c("blaOXA-24")
setDT(amrdb)[gene_name %in% oxa24, amr_gene_group := "blaOXA-24"]
## oxa-51
oxa51 <- c("blaOXA-510","blaOXA-217","blaOXA-69","blaOXA-91","blaOXA-98","blaOXA-66","blaOXA-120")
setDT(amrdb)[gene_name %in% oxa51, amr_gene_group := "blaOXA-51"]
## oxa-58
oxa58 <- c("blaOXA-58")
setDT(amrdb)[gene_name %in% oxa58, amr_gene_group := "blaOXA-58"]
## oxa-143
oxa143 <- c("blaOXA-143")
setDT(amrdb)[gene_name %in% oxa143, amr_gene_group := "blaOXA-143"]
## oxa-235
oxa235 <- c("blaOXA-235")
setDT(amrdb)[gene_name %in% oxa235, amr_gene_group := "blaOXA-235"]
## oxa-417
oxa417 <- c("blaOXA-417")
setDT(amrdb)[gene_name %in% oxa417, amr_gene_group := "blaOXA-417"]

### create Aminoglycoside resistance group
# first, create a vector of listed genes
aminoglyco <- c("aac(3)-IIa","aac(3)-IId","aac(3)-IV","aac(6')-Ib","aac(6')-Ib3","aac(6')-Ib-cr","aac(6')-Ib-Hangzhou",
                "aadA1","aadA16","aadA2","aadA22","aadA5","aadA8b","aadD",
                "ant(2'')-Ia","ant(3'')-Ia","ant(6)-Ia","ant(9)-Ia",
                "aph(2'')-Ia","aph(3')-Ia","aph(3'')-Ib","aph(3')-III","aph(3')-VIa","aph(4)-Ia","aph(6)-Id",
                "armA","rmtB")
# then we mark them as CHDL
setDT(amrdb)[gene_name %in% aminoglyco, hydrolysis_group := "Aminoglycoside Resistance"]

## aac 
aac <- c("aac(3)-IIa","aac(3)-IId","aac(3)-IV","aac(6')-Ib","aac(6')-Ib3","aac(6')-Ib-cr","aac(6')-Ib-Hangzhou")
setDT(amrdb)[gene_name %in% aac, amr_gene_group := "Aminoglycoside Modifying Enzymes (AMEs) - aac"]

## aad 
aad <- c("aadA1","aadA16","aadA2","aadA22","aadA5","aadA8b","aadD")
setDT(amrdb)[gene_name %in% aad, amr_gene_group := "Aminoglycoside Modifying Enzymes (AMEs) - aad"]

## ant 
aat <- c("ant(2'')-Ia","ant(3'')-Ia","ant(6)-Ia","ant(9)-Ia")
setDT(amrdb)[gene_name %in% aat, amr_gene_group := "Aminoglycoside Modifying Enzymes (AMEs) - aat"]

## aph
aph <- c("aph(2'')-Ia","aph(3')-Ia","aph(3'')-Ib","aph(3')-III","aph(3')-VIa","aph(4)-Ia","aph(6)-Id")
setDT(amrdb)[gene_name %in% aph, amr_gene_group := "Aminoglycoside Modifying Enzymes (AMEs) - aph"]

## amrA 
amrA <- c("armA")
setDT(amrdb)[gene_name %in% amrA, amr_gene_group := "16S RMTases - amrA"]

## rmtB
rmtB <- c("rmtB")
setDT(amrdb)[gene_name %in% rmtB, amr_gene_group := "16S RMTases - rmtB"]




### create Fluoroquinolones group
# first, create a vector of listed genes 
fluoro <- c("qnrB1","qnrB4","qnrB6","qnrB91","qnrS1","qnrVC6")
# then we mark them as flouroquinolones
setDT(amrdb)[gene_name %in% fluoro, hydrolysis_group := "Fluoroquinolones Resistances"]

## qnrB
qnrB <- c("qnrB1","qnrB4","qnrB6","qnrB91")
setDT(amrdb)[gene_name %in% qnrB, amr_gene_group := "qnrB"]

## qnrS
qnrS <- c("qnrS1")
setDT(amrdb)[gene_name %in% qnrS, amr_gene_group := "qnrS"]

## qnrVC
qnrVC <- c("qnrVC6")
setDT(amrdb)[gene_name %in% qnrVC, amr_gene_group := "qnrVC"]


### create Macrolides group 
# first, create a vector of listed genes 
macrolides <- c("erm(A)","erm(B)","erm(C)",
                "mph(A)","mph(E)","msr(E)")
# then we mark tham as macrolides resistance
setDT(amrdb)[gene_name %in% macrolides, hydrolysis_group := "Macrolides Resistances"]

## erm
erm <- c("erm(A)","erm(B)","erm(C)")
setDT(amrdb)[gene_name %in% erm, amr_gene_group := "erm"]

## mph
mph <- c("mph(A)","mph(E)")
setDT(amrdb)[gene_name %in% mph, amr_gene_group := "mph"]

## msr
msr <- c("msr(E)")
setDT(amrdb)[gene_name %in% msr, amr_gene_group := "msr"]

### create Trimethoprim/sulfamethoxazole group 
trime <- c("dfrA1","dfrA12","dfrA14","dfrA16","dfrA17","dfrA27","dfrA7","dfrE","dfrG",
           "sul1","sul2","sul3")
# mark the listed gene as trime 
setDT(amrdb)[gene_name %in% trime,hydrolysis_group := "Trimethoprim/Sulfamethoxazole Resistances"]

## dfr 
dfr <- c("dfrA1","dfrA12","dfrA14","dfrA16","dfrA17","dfrA27","dfrA7","dfrE","dfrG")
setDT(amrdb)[gene_name %in% dfr, amr_gene_group := "dfr"]

## sul
sul <- c("sul1","sul2","sul3")
setDT(amrdb)[gene_name %in% sul, amr_gene_group := "sul"]

### create Chloramphenicol 
chloram <- c("cat(pC221)","cat(pC233)","catA1","catB3",
             "floR","fexA","cmlA1","OqxA","OqxB")
# 
setDT(amrdb)[gene_name %in% chloram,hydrolysis_group := "Chloramphenicol Resistance"]

## cat 
cat <- c("cat(pC221)","cat(pC233)","catA1","catB3")
setDT(amrdb)[gene_name %in% cat, amr_gene_group := "cat"]

## floR
floR <- c("floR")
setDT(amrdb)[gene_name %in% floR, amr_gene_group := "floR"]

## fexA
fexA <- c("fexA")
setDT(amrdb)[gene_name %in% fexA, amr_gene_group := "fexA"]

## cml
cml <- c("cmlA1")
setDT(amrdb)[gene_name %in% cml, amr_gene_group := "cml"]

## Oqx
Oqx <- c("OqxA","OqxB")
setDT(amrdb)[gene_name %in% Oqx, amr_gene_group := "Oqx"]


### create Fosfomycin
fosformycin <- c("fosA","fosA3","fosA6")
setDT(amrdb)[gene_name %in% fosformycin,hydrolysis_group := "Fosformycin Resistance" ]

## fos
fos <- c("fosA","fosA3","fosA6")
setDT(amrdb)[gene_name %in% fos, amr_gene_group := "fos"]

### create Colistin
colistin <- c("mcr-1.1")
setDT(amrdb)[gene_name %in% colistin,hydrolysis_group := "Colistin Resistance" ]

## mcr
mcr <- c("mcr-1.1")
setDT(amrdb)[gene_name %in% mcr, amr_gene_group := "mcr"]

### create Rifampicin
rifampicin <- c("ARR-2","ARR-3")
setDT(amrdb)[gene_name %in% rifampicin,hydrolysis_group := "Rifampicin Resistance" ]

## ARR 
ARR <- c("ARR-2","ARR-3")
setDT(amrdb)[gene_name %in% ARR, amr_gene_group := "ARR"]

### create Tetracycline 
tetra <- c("tet(39)","tet(A)","tet(B)","tet(D)","tet(K)","tet(L)","tet(M)")
setDT(amrdb)[gene_name %in% tetra, hydrolysis_group := "Tetracycline Resistance"]

## tet 
tet <- c("tet(39)","tet(A)","tet(B)","tet(D)","tet(K)","tet(L)","tet(M)")
setDT(amrdb)[gene_name %in% tet, amr_gene_group := "tet"]




### create amr gene resistance exclusively to S.aureus 
saureus <- c("mecA","blaZ")
setDT(amrdb)[gene_name %in% saureus, hydrolysis_group := "AMR gene exclusive to S.aureus"]

## mecA 
mecA <- c("mecA")
setDT(amrdb)[gene_name %in% mecA, amr_gene_group := "mecA"]

## blaZ
blaZ <- c("blaZ")
setDT(amrdb)[gene_name %in% blaZ, amr_gene_group := "blaZ"]


##### export 
saveRDS(amrdb, file = paste0("Data/AMRdb/amr_genotype_db_updated_",Sys.Date(),".rds"))
