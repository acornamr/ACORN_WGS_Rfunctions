# this script is used to plot the viruence factors of four ACORN targeted pathogen 
# written and maintained by Tung Trinh
# April 19th 2025
# version 1.2
# modified on April 28th 2025
# for more informaiton, please contact to tungts@oucru.org
####################################################################################################
##### built-in function and required packages #####
library(dplyr)
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(data.table)
library(ComplexHeatmap)
##### read the data #####
### metadata 
meta <- read_excel("Data/ACORN2_VN001_WGS_Cumulative_V2_AKP.xlsx", sheet = 1 )
meta <- meta %>% mutate(clinical_condition = case_when(
  clinical_severity_score== 0 ~ "None",
  clinical_severity_score== 1 ~ "Mild",
  clinical_severity_score== 2 ~ "Moderate",
  clinical_severity_score== 3 ~ "Severe"
))
## amr profile - acine 
staramr_acine <- read_excel("Data/assemblies/Hanoi/staramr_acine/results.xlsx")
# create id 
staramr_acine <- staramr_acine %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+"))
## amr profile - stap 
staramr_kleb <- read_excel("Data/assemblies/Hanoi/staramr_kleb_plus/results.xlsx")
# # # create id 
staramr_kleb <- staramr_kleb %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+"))
# ## amr profile - ecoli
staramr_ecol <- read_excel("Data/assemblies/Hanoi/staramr_ecol/staramr_ecol_results.xlsx")
# # create id 
staramr_ecol <- staramr_ecol %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+"))
# ## amr profile - staph
staramr_stap <- read_excel("Data/assemblies/Hanoi/staramr_staph/results.xlsx")
# # create id 
staramr_stap <- staramr_staph %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+"))
### read data from vfdb 
vfdb_detail <- read.table("Data/vfdb.tab", sep = "\t")
vfdb_detail <- vfdb_detail %>% mutate(VFCID = str_extract(V14,"VF[0-9]+"))
vfdb <- read.table("Data/summary_vfdb.tab", sep = "\t", header = T)
vfdb <- vfdb %>% mutate(seqid = str_extract(FILE,"A[0-9]-[0-9]+")) 
# read the vfdb data dictionary 
vf <- read_excel("Data/VFs.xls")
colnames(vf) <- vf[1,] %>% as.vector()
vf <- vf[-1,]
vfdb_detail <- vfdb_detail %>% select(Gene = V6,VFCID) %>% distinct(.keep_all = T)
vfdb_detail <- left_join(vfdb_detail,vf, by = c("VFCID" = "VFID"))
vfdb_detail <- vfdb_detail %>% mutate(Gene = str_replace(Gene,"/","."))

# join with metadata 
vfdb_work <- vfdb %>% left_join(.,meta %>% select(`SEQ-ID`,surveillance_category,clinical_condition), by = c("seqid" = "SEQ-ID"))


################### ACINE 
acine_id <- staramr_acine %>% pull(seqid)

vf_acine <- vfdb_work %>% filter(seqid %in% acine_id) %>% mutate(title = "Acinetobacter baumannii") %>%
  left_join(., staramr_acine %>% select(seqid,ST = `Sequence Type`), by = c("seqid" = "seqid"))
ind_acine  <- which((vf_acine %>% map(~sum(. ==".")) %>% unlist() ) == nrow(vf_acine)) %>% as.vector()
vf_acine <- vf_acine %>% select(-all_of(ind_acine))

# subset virulence factor only 
work_acine <- vf_acine[,c(3:138)]
work_acine[work_acine == "."] <- F
work_acine[work_acine != F ] <- T

row.names(work_acine) <- vf_acine$seqid
acine_m <- work_acine %>% as.matrix()

# column annotation 
target_vf_acine <- colnames(work_acine)

anno_acine  <- vfdb_detail %>% filter(Gene %in% target_vf_acine) %>% select(Gene,VFcategory) %>% 
  distinct(.keep_all = T) %>% arrange(VFcategory) %>% slice(-84,-106)
col_order_hm <- anno_acine %>% pull(Gene)

colA_acine <- data.frame(Group = anno_acine$VFcategory)
rownames(colA_acine) <- anno_acine %>% pull(Gene)

# create a data for row annotation
rowA_abau <- data.frame(source = vf_acine$ST) %>% arrange(source)
rownames(rowA_abau) <- vf_acine %>% arrange(ST) %>% pull(seqid)

annotated_row <- rowAnnotation(df = rowA_abau)
annotated_col <- HeatmapAnnotation(df = colA_acine)

row_order_hm <-  vf_acine %>% arrange(ST) %>% pull(seqid)
acine_m <- acine_m[row_order_hm,col_order_hm]

p <- Heatmap(acine_m,column_title = "Acinobacter baum",col = c("gray80","brown"),
        # row_dend_reorder = F,
        row_order = row_order_hm,
        column_order = col_order_hm,
        clustering_method_rows = "ward.D",
        column_labels = rep("",(n = ncol(acine_m))),
        top_annotation = annotated_col,
        # clustering_distance_columns =  "euclidean",
        heatmap_legend_param = list(
          title = "Genotype", labels = c("No","Yes")
        )
)

jpeg(filename = "Figures/acine_vf.jpg", units = "in", res = 300, height = 6,width = 10)
draw(p + annotated_row,annotation_legend_side = "top")
dev.off()

######################## ECOL
ecol_id <- staramr_ecol %>% pull(seqid)

vf_ecol <- vfdb_work %>% filter(seqid %in% ecol_id) %>% mutate(title = "Escherichia coli") %>%
  left_join(., staramr_ecol %>% select(seqid,ST = `Sequence Type`), by = c("seqid" = "seqid"))
ind_ecol  <- which((vf_ecol %>% map(~sum(. ==".")) %>% unlist() ) == nrow(vf_ecol)) %>% as.vector()
vf_ecol <- vf_ecol %>% select(-all_of(ind_ecol))

# subset virulence factor only 
work_ecol <- vf_ecol[,c(3:209)]
work_ecol[work_ecol == "."] <- F
work_ecol[work_ecol != F ] <- T

row.names(work_ecol) <- vf_ecol$seqid
ecol_m <- work_ecol %>% as.matrix()

# column annotation 
target_vf_ecol <- colnames(work_ecol)

anno_ecol  <- vfdb_detail %>% filter(Gene %in% target_vf_ecol) %>% select(Gene,VFcategory) %>% 
  distinct(.keep_all = T) %>% arrange(VFcategory) %>% slice(-141,-142)
col_order_hm <- anno_ecol %>% pull(Gene)

colA_ecol <- data.frame(Group = anno_ecol$VFcategory)
rownames(colA_ecol) <- anno_ecol %>% pull(Gene)

# create a data for row annotation
rowA_ecol <- data.frame(source = vf_ecol$clinical_condition) %>% arrange(source)
rownames(rowA_ecol) <- vf_ecol %>% arrange(clinical_condition) %>% pull(seqid)

annotated_row <- rowAnnotation(df = rowA_ecol)
annotated_col <- HeatmapAnnotation(df = colA_ecol)

row_order_hm <-  vf_ecol %>% arrange(clinical_condition) %>% pull(seqid)
ecol_m <- ecol_m[row_order_hm,col_order_hm]

p <- Heatmap(ecol_m,column_title = "Escherichia coli",col = c("gray80","brown"),
             # row_dend_reorder = F,
             row_order = row_order_hm,
             column_order = col_order_hm,
             clustering_method_rows = "ward.D",
             column_labels = rep("",(n = ncol(ecol_m))),
             top_annotation = annotated_col,
             # clustering_distance_columns =  "euclidean",
             heatmap_legend_param = list(
               title = "Genotype", labels = c("No","Yes")
             )
)

jpeg(filename = "Figures/ecol_vf.jpg", units = "in", res = 300, height = 6,width = 10)
draw(p + annotated_row,annotation_legend_side = "top")
dev.off()


######################## KLeb
kleb_id <- staramr_kleb %>% pull(seqid)

vf_kleb <- vfdb_work %>% filter(seqid %in% kleb_id) %>% mutate(title = "Klebsiella pneumoniae") %>%
  left_join(., staramr_kleb %>% select(seqid,ST = `Sequence Type`), by = c("seqid" = "seqid"))
ind_kleb  <- which((vf_kleb %>% map(~sum(. ==".")) %>% unlist() ) == nrow(vf_kleb)) %>% as.vector()
vf_kleb <- vf_kleb %>% select(-all_of(ind_kleb))

# subset virulence factor only 
work_kleb <- vf_kleb[,c(3:138)]
work_kleb[work_kleb == "."] <- F
work_kleb[work_kleb != F ] <- T

row.names(work_kleb) <- vf_kleb$seqid
kleb_m <- work_kleb %>% as.matrix()

# column annotation 
target_vf_kleb <- colnames(work_kleb)

anno_kleb  <- vfdb_detail %>% filter(Gene %in% target_vf_kleb) %>% select(Gene,VFcategory) %>% 
  distinct(.keep_all = T) %>% arrange(VFcategory) %>% slice(-141,-142)
col_order_hm <- anno_kleb %>% pull(Gene)

colA_kleb <- data.frame(Group = anno_kleb$VFcategory)
rownames(colA_kleb) <- anno_kleb %>% pull(Gene)

# create a data for row annotation
rowA_kleb <- data.frame(source = vf_kleb$ST) %>% arrange(source)
rownames(rowA_kleb) <- vf_kleb %>% arrange(ST) %>% pull(seqid)

annotated_row <- rowAnnotation(df = rowA_kleb)
annotated_col <- HeatmapAnnotation(df = colA_kleb)

row_order_hm <-  vf_kleb %>% arrange(clinical_condition) %>% pull(seqid)
kleb_m <- kleb_m[row_order_hm,col_order_hm]

p <- Heatmap(kleb_m,column_title = "Klebsiella pneumoniae",col = c("gray80","brown"),
             # row_dend_reorder = F,
             row_order = row_order_hm,
             column_order = col_order_hm,
             clustering_method_rows = "ward.D",
             column_labels = rep("",(n = ncol(kleb_m))),
             top_annotation = annotated_col,
             # clustering_distance_columns =  "euclidean",
             heatmap_legend_param = list(
               title = "Genotype", labels = c("No","Yes")
             )
)

jpeg(filename = "Figures/kleb_vf.jpg", units = "in", res = 300, height = 6,width = 10)
draw(p + annotated_row,annotation_legend_side = "top")
dev.off()

######################## stpah
stap_id <- staramr_stap %>% pull(seqid)

vf_stap <- vfdb_work %>% filter(seqid %in% stap_id) %>% mutate(title = "Staphylococcus aureus") %>%
  left_join(., staramr_stap %>% select(seqid,ST = `Sequence Type`), by = c("seqid" = "seqid"))
ind_stap  <- which((vf_stap %>% map(~sum(. ==".")) %>% unlist() ) == nrow(vf_stap)) %>% as.vector()
vf_stap <- vf_stap %>% select(-all_of(ind_stap))

# subset virulence factor only 
work_stap <- vf_stap[,c(3:100)]
work_stap[work_stap == "."] <- F
work_stap[work_stap != F ] <- T

row.names(work_stap) <- vf_stap$seqid
stap_m <- work_stap %>% as.matrix()

# column annotation 
target_vf_stap <- colnames(work_stap)

anno_stap  <- vfdb_detail %>% filter(Gene %in% target_vf_stap) %>% select(Gene,VFcategory) %>% 
  distinct(.keep_all = T) %>% arrange(VFcategory) %>% slice(-141,-142)
col_order_hm <- anno_stap %>% pull(Gene)

colA_stap <- data.frame(Group = anno_stap$VFcategory)
rownames(colA_stap) <- anno_stap %>% pull(Gene)

# create a data for row annotation
rowA_stap <- data.frame(source = vf_stap$ST) %>% arrange(source)
rownames(rowA_stap) <- vf_stap %>% arrange(ST) %>% pull(seqid)

annotated_row <- rowAnnotation(df = rowA_stap)
annotated_col <- HeatmapAnnotation(df = colA_stap)

row_order_hm <-  vf_stap %>% arrange(clinical_condition) %>% pull(seqid)
stap_m <- stap_m[row_order_hm,col_order_hm]

p <- Heatmap(stap_m,column_title = "Staphylococcus aureus",col = c("gray80","brown"),
             # row_dend_reorder = F,
             row_order = row_order_hm,
             column_order = col_order_hm,
             clustering_method_rows = "ward.D",
             column_labels = rep("",(n = ncol(stap_m))),
             top_annotation = annotated_col,
             # clustering_distance_columns =  "euclidean",
             heatmap_legend_param = list(
               title = "Genotype", labels = c("No","Yes")
             )
)

jpeg(filename = "Figures/stap_vf.jpg", units = "in", res = 300, height = 6,width = 10)
draw(p + annotated_row,annotation_legend_side = "top")
dev.off()

