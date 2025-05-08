# this script is used to analysis the amr profile and virulence factors of ACORN target pathogen.
# Feb 17th 2025 
# version 1.0 
# Written and maintained by Tung Trinh
# For more information, please contact to tungts@oucru.org
####################################################################################################
##### required packages and built-in functions #####
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(data.table)
##### read the data #####
### staramr 
staramr <- read_excel("Data/staramr_hanoi/results.xlsx", sheet = 1)
### metadata 
meta <- read_excel("Data/ACORN2_VN001_WGS_Cumulative_V2_AKP.xlsx", sheet = 1 )
##### transform the data #####
# create the id 
staramr <- staramr %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+")) 
# select variables to include into the analysis 
# genotype 
amr <- staramr %>% select(seqid, Genotype,Scheme) %>% separate_longer_delim(Genotype,delim = ",") %>% 
  mutate(Genotype = str_replace(Genotype," ",""))
# select data to merge 
amr <- meta %>% select(`SEQ-ID`,surveillance_category,specgroup...16) %>% left_join(amr,., by = c("seqid" = "SEQ-ID"))
##### plot the result ######
# calculate the distribution
amr_dis <- amr %>% group_by(Genotype,surveillance_category,Scheme) %>% summarise(count = n()) %>% ungroup() %>% 
  filter(Genotype != "None")
# fix the name 
setDT(amr_dis)[Genotype == "OqxA", Genotype := "oqxA"]
setDT(amr_dis)[Genotype == "OqxB", Genotype := "oqxB"]
# split the result based on target pathogen
# ecol 
ecol_amr <- amr_dis %>% filter(Scheme == "ecoli_achtman_4") %>% group_by(Genotype) %>% 
  mutate(total = sum(count),
           y = 1.02) 
ecol_amr$title <- "Escherichia coli" 
# kleb 
kleb_amr <- amr_dis %>% filter(Scheme == "klebsiella") %>% group_by(Genotype) %>% 
  mutate(total = sum(count),
         y = 1.02) 
kleb_amr$title <- "Klebsiella pneumoniae"
# abaum
abau_amr <- amr_dis %>% filter(Scheme == "abaumannii_2") %>% group_by(Genotype) %>% 
  mutate(total = sum(count),
         y = 1.02) 
abau_amr$title <- "Acinetobacter baumannii"
# staph
staph_amr <- amr_dis %>% filter(Scheme == "saureus") %>% group_by(Genotype) %>% 
  mutate(total = sum(count),
         y = 1.02) 
staph_amr$title <- "Staphylococcus aureus"
### AMR genotype distribution 
# e coli
p1 <-
  ggplot(ecol_amr, aes(x=Genotype,y = count, fill=surveillance_category)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  # geom_bar(stat="identity", position="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  scale_fill_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"))+
  xlab("AMR Genotype")+ylab("Percentage")+
  geom_text(aes(x=Genotype, y=y, label=as.factor(total)), vjust=-0.1)+
facet_wrap(~title)
  
# kleb  
p2 <- ggplot(kleb_amr, aes(x=Genotype,y = count, fill=surveillance_category)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  scale_fill_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"))+
  geom_text(aes(x=Genotype, y=y, label=as.factor(total)), vjust=-0.1)+
  xlab("AMR Genotype")+ylab("Percentage")+facet_wrap(~title)
# abaum
p3 <- ggplot(abau_amr, aes(x=Genotype,y = count, fill=surveillance_category)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  scale_fill_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"))+
  geom_text(aes(x=Genotype, y=y, label=as.factor(total)), vjust=-0.1)+
  xlab("AMR Genotype")+ylab("Percentage")+facet_wrap(~title)
# staph
p4 <- ggplot(staph_amr, aes(x=Genotype,y = count, fill=surveillance_category)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  scale_fill_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"))+
  geom_text(aes(x=Genotype, y=y, label=as.factor(total)), vjust=-0.1)+
  xlab("AMR Genotype")+ylab("Percentage")+facet_wrap(~title)
# get the legend
p_legend <- ggplot(ecol_amr, aes(x=Genotype,y = count, fill=surveillance_category)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  # geom_bar(stat="identity", position="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")+
  scale_fill_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"), label = c("Community-acquired Infection (CAI)",
                    "Hospital-acquired Infection (HAI)",
                    "Healthcare-associated Infection (HCAI)"), name ="")+
  xlab("AMR Genotype")+ylab("Percentage")+facet_wrap(~title)

p1_legend <- cowplot::get_plot_component(p_legend,"guide-box-bottom",return_all = T)
combined_plot <- grid.arrange(p1,p2,p3,p4, ncol = 2)

jpeg(filename = "Figures/amr_dis.jpeg", res = 300, units = "in",height = 10,width = 22)
grid.arrange(combined_plot,p1_legend, nrow = 2, heights = c(10,1))
dev.off()


### predicted phenotype 
# read the database 
amr_db <- readRDS("Data/AMRdb/amr_genotype_db_2025-02-14.rds") %>% mutate(Gene = str_replace(Gene,"_[0-9]",""))
amr_db <- amr_db %>% add_row(ab_subclass = "quinolone",
                   Gene = "qnrB91",
                   Drug = "ciprofloxacin",
                   aware_subclass = "Fluoroquinolones",
                   `ATC code` = "J01MA02",
                   Category = "Watch")

amr_db <- amr_db %>% add_row(ab_subclass = "trimethoprim",
                             Gene = "dfrE",
                             Drug = "trimethoprim",
                             aware_subclass = "Trimethoprim-derivatives",
                             `ATC code` = "J01EA01",
                             Category = "Access")
amr_db <- amr_db %>% filter(route == "oral") %>% select(-Accession,-route) %>% distinct()
# merge with the data 
amr_dis_drug <- left_join(amr_dis,amr_db, by = c("Genotype" = "Gene")) %>% distinct()
# resummarise the table 
amr_dis_drug <- amr_dis_drug %>% group_by(surveillance_category,Scheme,Drug) %>% summarise(total = sum(count)) %>% 
  ungroup() %>% filter(!is.na(Drug))

# split the result based on target pathogen
# ecol 
ecol_ab <- amr_dis_drug %>% filter(Scheme == "ecoli_achtman_4") %>% group_by(Drug) %>% 
  mutate(n = sum(total),
         y = 1.02) 
ecol_ab$title <- "Escherichia coli"
# kleb 
kleb_ab <- amr_dis_drug %>% filter(Scheme == "klebsiella") %>% group_by(Drug) %>% 
  mutate(n = sum(total),
         y = 1.02) 
kleb_ab$title <- "Klebsiella pneumoniae"
# abaum
abau_ab <- amr_dis_drug %>% filter(Scheme == "abaumannii_2") %>% group_by(Drug) %>% 
  mutate(n = sum(total),
         y = 1.02) 
abau_ab$title <- "Acinetobacter baumannii"
# staph
staph_ab <- amr_dis_drug %>% filter(Scheme == "saureus") %>% group_by(Drug) %>% 
  mutate(n = sum(total),
         y = 1.02) 
staph_ab$title <- "Staphylococcus aureus"

### AMR genotype distribution 
# e coli
p1 <-
  ggplot(ecol_ab, aes(x=Drug,y = total, fill=surveillance_category)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  # geom_bar(stat="identity", position="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        # axis.title.x = element_text(vjust=-0.5),
        legend.position = "none")+
  scale_fill_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"))+
  geom_text(aes(x=Drug, y=y, label=as.factor(n)), vjust=-0.1)+
  xlab("AMR Predicted Phenotype")+ylab("Percentage")+facet_wrap(~title)

p2 <-
  ggplot(kleb_ab, aes(x=Drug,y = total, fill=surveillance_category)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  # geom_bar(stat="identity", position="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        # axis.title.x = element_text(vjust=-0.5),
        legend.position = "none")+
  scale_fill_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"))+
  geom_text(aes(x=Drug, y=y, label=as.factor(n)), vjust=-0.1)+
  xlab("AMR Predicted Phenotype")+ylab("Percentage")+facet_wrap(~title)

p3 <-
  ggplot(abau_ab, aes(x=Drug,y = total, fill=surveillance_category)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  # geom_bar(stat="identity", position="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        # axis.title.x = element_text(vjust=-0.5),
        legend.position = "none")+
  scale_fill_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"))+
  geom_text(aes(x=Drug, y=y, label=as.factor(n)), vjust=-0.1)+
  xlab("AMR Predicted Phenotype")+ylab("Percentage")+facet_wrap(~title)
p4 <-
  ggplot(staph_ab, aes(x=Drug,y = total, fill=surveillance_category)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  # geom_bar(stat="identity", position="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        # axis.title.x = element_text(vjust=-0.5),
        legend.position = "none")+
  scale_fill_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"))+
  geom_text(aes(x=Drug, y=y, label=as.factor(n)), vjust=-0.1)+
  xlab("AMR Predicted Phenotype")+ylab("Percentage")+facet_wrap(~title)

combined_plot <- grid.arrange(p1,p2,p3,p4, ncol = 2)

jpeg(filename = "Figures/amr_drug.jpeg", res = 300, units = "in",height = 11,width = 14)
grid.arrange(combined_plot,p1_legend, nrow = 2, heights = c(10,1))
dev.off()

#
staramr <- staramr %>% mutate(total_amr_gene = str_count(Genotype,"[,]"),
                   total_predicted_drug_resisted = str_count(`Predicted Phenotype`,"[,]")) 

dat <- meta %>% left_join(staramr,., by = c("seqid" = "SEQ-ID"))

dat <- dat %>% mutate(title = if_else(Scheme == "abaumannii_2","Acinetobacter baumannii","Escherichia coli"),
               title = if_else(Scheme == "saureus","Staphylococcus aureus",title),
               title = if_else(Scheme == "klebsiella","Klebsiella pneumoniae",title))
library(ggpubr)

comp_list <- list(c("HCAI", "CAI"),c("CAI","HAI"),c("HAI","HCAI"))

# total amr gene 
jpeg(filename = "Figures/total_amr_gene.jpg", unit = "in", res = 300, height = 10,width = 12)
ggviolin(dat, x = "surveillance_category", y = "total_amr_gene", 
         add = "boxplot", # add boxplot
         add.params = list(fill = "white"), # configure the inside colour of box plot
         fill = "surveillance_category",
         # palette = c("#6B58EFFF","#1F8F99FF","#D60C00FF")
) + stat_compare_means(comparisons = comp_list) +
  xlab("")+ylab("Total AMR gene have")+theme_bw()+
  theme(legend.position = "top")+
  scale_fill_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"),
                    label = c("Community-acquired Infection (CAI)",
                              "Hospital-acquired Infection (HAI)",
                              "Healthcare-associated Infection (HCAI)"), name ="")+
  stat_compare_means(label.x = 0.75, label.y = 75)+
  facet_wrap(~title)
dev.off()

jpeg(filename = "Figures/total_phenotype_gene.jpg", unit = "in", res = 300, height = 10,width = 12)
ggviolin(dat, x = "surveillance_category", y = "total_predicted_drug_resisted", 
         add = "boxplot", # add boxplot
         add.params = list(fill = "white"), # configure the inside colour of box plot
         fill = "surveillance_category",
         # palette = c("#6B58EFFF","#1F8F99FF","#D60C00FF")
) + stat_compare_means(comparisons = comp_list) +
  xlab("")+ylab("Total predicted antibiotics resisted to")+theme_bw()+
  theme(legend.position = "top")+
  scale_fill_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"),
                    label = c("Community-acquired Infection (CAI)",
                              "Hospital-acquired Infection (HAI)",
                              "Healthcare-associated Infection (HCAI)"), name ="")+
  stat_compare_means(label.x = 0.75, label.y = 75)+
  facet_wrap(~title)
dev.off()

ggviolin(dat, x = "ho_discharge_status", y = "total_predicted_drug_resisted", 
         add = "boxplot", # add boxplot
         add.params = list(fill = "white"), # configure the inside colour of box plot
         fill = "ho_discharge_status",
         palette = c("#240E31FF", "#CB6BCEFF", "#468892FF", "#74F3D3FF")
) + 
  stat_compare_means() +
  theme_bw()+
  theme(legend.position = "top")+
  facet_wrap(~Scheme)# fix the position of p value 

# ggviolin(dat, x = "ho_discharge_status", y = "total_amr_gene", 
#          add = "boxplot", # add boxplot
#          add.params = list(fill = "white"), # configure the inside colour of box plot
#          fill = "ho_discharge_status",
#          palette = c("#240E31FF", "#CB6BCEFF", "#468892FF", "#74F3D3FF")
# ) + 
#   stat_compare_means() +
#   theme_bw()+
#   theme(legend.position = "top")+
#   facet_wrap(~Scheme)# fix the position of p value 


####
# calculate the distribution
amr_sample <- amr %>% group_by(Genotype,specgroup...16,Scheme) %>% summarise(count = n()) %>% ungroup() %>% 
  filter(Genotype != "None")
# fix the name 
setDT(amr_sample)[Genotype == "OqxA", Genotype := "oqxA"]
setDT(amr_sample)[Genotype == "OqxB", Genotype := "oqxB"]


amr_sample_drug <- left_join(amr_sample,amr_db, by = c("Genotype" = "Gene")) %>% distinct()
# resummarise the table 
amr_dis_drug <- amr_sample_drug %>% group_by(specgroup...16,Scheme,Drug) %>% summarise(total = sum(count)) %>% 
  ungroup() %>% filter(!is.na(Drug))

# split the result based on target pathogen
# ecol 
ecol_amr <- amr_sample %>% filter(Scheme == "ecoli_achtman_4") %>% group_by(Genotype) %>% 
  mutate(n = sum(count),
         y = 1.02) 
ecol_amr$title <- "Escherichia coli"
# kleb 
kleb_amr <- amr_sample %>% filter(Scheme == "klebsiella") %>% group_by(Genotype) %>% 
  mutate(n = sum(count),
         y = 1.02) 
kleb_amr$title <- "Klebsiella pneumoniae"
# abaum
abau_amr <- amr_sample %>% filter(Scheme == "abaumannii_2") %>% group_by(Genotype) %>% 
  mutate(n = sum(count),
         y = 1.02) 
abau_amr$title <- "Acinetobacter baumannii"
# staph
staph_amr <- amr_sample %>% filter(Scheme == "saureus") %>% group_by(Genotype) %>% 
  mutate(n = sum(count),
         y = 1.02) 
staph_amr$title <- "Staphylococcus aureus"

### AMR genotype distribution 
# e coli
p1 <-
  ggplot(ecol_amr, aes(x=Genotype,y = count, fill=specgroup...16)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  # geom_bar(stat="identity", position="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("#F4F4A9FF", "#6D9A58FF","#4E4E36FF","#353525FF"))+
  geom_text(aes(x=Genotype, y=y, label=as.factor(n)), vjust=-0.1)+
  theme(legend.position = "none")+
  xlab("AMR Genotype")+ylab("Percentage")+facet_wrap(~title)
# kleb  
p2 <- ggplot(kleb_amr, aes(x=Genotype,y = count, fill=specgroup...16)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  # geom_bar(stat="identity", position="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("#F4F4A9FF", "#6D9A58FF","#4E4E36FF", "#F6F6F6FF", "#C8CC9FFF"))+
  geom_text(aes(x=Genotype, y=y, label=as.factor(n)), vjust=-0.1)+
  theme(legend.position = "none")+
  xlab("AMR Genotype")+ylab("Percentage")+facet_wrap(~title)
# abaum
p3 <- ggplot(abau_amr, aes(x=Genotype,y = count, fill=specgroup...16)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  # geom_bar(stat="identity", position="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("#F4F4A9FF", "#6D9A58FF","#4E4E36FF", "#F6F6F6FF", "#C8CC9FFF","#353525FF"))+
  geom_text(aes(x=Genotype, y=y, label=as.factor(n)), vjust=-0.1)+
  theme(legend.position = "none")+
  xlab("AMR Genotype")+ylab("Percentage")+facet_wrap(~title)
# staph
p4 <- ggplot(staph_amr, aes(x=Genotype,y = count, fill=specgroup...16)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  # geom_bar(stat="identity", position="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("#F4F4A9FF", "#E5CA28FF", "#4E4E36FF", "#F6F6F6FF", "#C8CC9FFF","#353525FF"))+
  geom_text(aes(x=Genotype, y=y, label=as.factor(n)), vjust=-0.1)+
  theme(legend.position = "none")+
  xlab("AMR Genotype")+ylab("Percentage")+facet_wrap(~title)
# get the legend
p_legend <- ggplot(amr_sample, aes(x=Genotype,y = count, fill=specgroup...16)) + 
  geom_bar(stat="identity", position="fill")+theme_bw()+
  # geom_bar(stat="identity", position="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("#F4F4A9FF", "#6D9A58FF", "#E5CA28FF", "#4E4E36FF", "#F6F6F6FF", "#C8CC9FFF","#353525FF"), name = "Sample Type")+
  # geom_text(aes(x=Genotype, y=y, label=as.factor(n)), vjust=-0.1)+
  theme(legend.position = "bottom")+
  xlab("AMR Genotype")+ylab("Percentage")

p1_legend <- cowplot::get_plot_component(p_legend,"guide-box-bottom",return_all = T)
combined_plot <- grid.arrange(p1,p2,p3,p4, ncol = 2)

jpeg(filename = "Figures/amr_sampletype.jpeg", res = 300, units = "in",height = 10,width = 22)
grid.arrange(combined_plot,p1_legend, nrow = 2, heights = c(10,1))
dev.off()


####################################################################################################
### read data from vfdb 
vfdb_detail <- read.table("Data/vfdb.tab", sep = "\t")
vfdb_detail <- vfdb_detail %>% mutate(VFCID = str_extract(V14,"VF[0-9]+"))
vfdb <- read.table("Data/summary_vfdb.tab", sep = "\t", header = T)
# read the vfdb data dictionary 
vf <- read_excel("Data/VFs.xls")
colnames(vf) <- vf[1,] %>% as.vector()
vf <- vf[-1,]
vfdb_detail <- vfdb_detail %>% select(Gene = V6,VFCID) %>% distinct(.keep_all = T)
vfdb_detail <- left_join(vfdb_detail,vf, by = c("VFCID" = "VFID"))
vfdb_detail <- vfdb_detail %>% mutate(Gene = str_replace(Gene,"/","."))


vfdb <- vfdb %>% mutate(seqid = str_extract(FILE,"A[0-9]-[0-9]+")) 
# create the id 
staramr <- staramr %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+")) 
# select variables to include into the analysis 
# join with data 
vfdb_work <- staramr %>% select(seqid,Scheme) %>% left_join(vfdb,., by = c("seqid" = "seqid"))
# join with metadata 
vfdb_work <- vfdb_work %>% left_join(.,meta %>% select(`SEQ-ID`,surveillance_category), by = c("seqid" = "SEQ-ID"))
# subset data by species
ecol_vfdb <- vfdb_work %>% filter(Scheme ==  "ecoli_achtman_4") %>% mutate(title = "Escherichia coli")
kleb_vfdb <- vfdb_work %>% filter(Scheme ==  "klebsiella") %>% mutate(title = "Klebsiella pneumoniae")
stap_vfdb <- vfdb_work %>% filter(Scheme ==  "saureus") %>% mutate(title = "Staphylococcus aureus")
abau_vfdb <- vfdb_work %>% filter(Scheme ==  "abaumannii_2") %>% mutate(title = "Acinetobacter baumannii")
# 
library(ComplexHeatmap)
### working on ecol 
# find the column to discards
ind_ecol <- which((ecol_vfdb %>% map(~sum(. ==".")) %>% unlist() ) == nrow(ecol_vfdb)) %>% as.vector()
ecol_vfdb <- ecol_vfdb %>% select(-all_of(ind_ecol))
# subset virulence factor only 
ecol_vf <- ecol_vfdb[,c(3:209)]
ecol_vf[ecol_vf == "."] <- F
ecol_vf[ecol_vf != F ] <- T
ecol_work <- ecol_vfdb %>% select(seqid) %>% cbind(.,ecol_vf)
row.names(ecol_work) <- ecol_work$seqid
ecol_m <- ecol_work %>% select(-seqid) %>% as.matrix()

# column annotation 
target_vf_ecol <- colnames(ecol_work)[-1]
vf_ecol <- vfdb_detail %>% filter(Gene %in% target_vf_ecol) %>% select(Gene,VFcategory) %>% 
  distinct(.keep_all = T) %>% arrange(VFcategory) %>% slice(-141,-142)
col_order_hm <- vf_ecol %>% pull(Gene)

colA_ecol <- data.frame(Group = vf_ecol$VFcategory)
rownames(colA_ecol) <- vf_ecol %>% pull(Gene)



# create a data for row annotation
rowA_ecol <- data.frame(source = ecol_vfdb$surveillance_category) %>% arrange(source)
rownames(rowA_ecol) <- ecol_vfdb %>% arrange(surveillance_category) %>% pull(seqid)

annotated_row <- rowAnnotation(df = rowA_ecol)
annotated_col <- HeatmapAnnotation(df = colA_ecol)

row_order_hm <-  ecol_vfdb %>% arrange(surveillance_category) %>% pull(seqid)
ecol_m <- ecol_m[row_order_hm,col_order_hm]
  
p <-
  Heatmap(ecol_m,column_title = "Escherichia coli",col = c("gray80","brown"),
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

jpeg(filename = "Figures/vf_ecol.jpg", units = "in", res = 300, height = 6,width = 10)
draw(p + annotated_row,annotation_legend_side = "top")
dev.off()

### working on kleb 
# find the column to discards
ind_kleb <- which((kleb_vfdb %>% map(~sum(. ==".")) %>% unlist() ) == nrow(kleb_vfdb)) %>% as.vector()
kleb_vfdb <- kleb_vfdb %>% select(-all_of(ind_kleb))
# subset virulence factor only 
kleb_vf <- kleb_vfdb[,c(3:138)]
kleb_vf[kleb_vf == "."] <- F
kleb_vf[kleb_vf != F ] <- T
kleb_work <- kleb_vfdb %>% select(seqid) %>% cbind(.,kleb_vf)
row.names(kleb_work) <- kleb_work$seqid
kleb_m <- kleb_work %>% select(-seqid) %>% as.matrix()
# column annotation 
target_vf_kleb <- colnames(kleb_work)[-1]
vf_kleb <- vfdb_detail %>% filter(Gene %in% target_vf_kleb) %>% select(Gene,VFcategory) %>% 
  distinct(.keep_all = T) %>% arrange(VFcategory) %>% slice(-86,-87)
col_order_hm <- vf_kleb %>% pull(Gene)

colA_kleb <- data.frame(Group = vf_kleb$VFcategory)
rownames(colA_kleb) <- vf_kleb %>% pull(Gene)

# create a data for row annotation
rowA_kleb <- data.frame(source = kleb_vfdb$surveillance_category) %>% arrange(source)
rownames(rowA_kleb) <- kleb_vfdb %>% arrange(surveillance_category) %>% pull(seqid)

annotated_row <- rowAnnotation(df = rowA_kleb)
annotated_col <- HeatmapAnnotation(df = colA_kleb)

row_order_hm <-  kleb_vfdb %>% arrange(surveillance_category) %>% pull(seqid)
kleb_m <- kleb_m[row_order_hm,col_order_hm]

p <-
  Heatmap(kleb_m,column_title = "Klebsiella pneumoniae",col = c("gray80","brown"),
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

jpeg(filename = "Figures/vf_kleb.jpg", units = "in", res = 300, height = 6,width = 10)
draw(p + annotated_row,annotation_legend_side = "top")
dev.off()

### working on staph 
# find the column to discards
ind_stap <- which((stap_vfdb%>% map(~sum(. ==".")) %>% unlist() ) == nrow(stap_vfdb)) %>% as.vector()
stap_vfdb <- stap_vfdb %>% select(-all_of(ind_stap))
# subset virulence factor only 
stap_vf <- stap_vfdb[,c(3:100)]
stap_vf[stap_vf == "."] <- F
stap_vf[stap_vf != F ] <- T
stap_work <- stap_vfdb %>% select(seqid) %>% cbind(.,stap_vf)
row.names(stap_work) <- stap_work$seqid
stap_m <- stap_work %>% select(-seqid) %>% as.matrix()
# column annotation 
target_vf_stap <- colnames(stap_work)[-1]
vf_stap <- vfdb_detail %>% mutate(Gene = str_replace(Gene,"-",".")) %>%  filter(Gene %in% target_vf_stap) %>% select(Gene,VFcategory) %>% 
  distinct(.keep_all = T) %>% arrange(VFcategory) 
# %>% slice(-86,-87)
col_order_hm <- vf_stap %>% pull(Gene)

colA_stap <- data.frame(Group = vf_stap$VFcategory)
rownames(colA_stap) <- vf_stap %>% pull(Gene)

# create a data for row annotation
rowA_stap <- data.frame(source = stap_vfdb$surveillance_category) %>% arrange(source)
rownames(rowA_stap) <- stap_vfdb %>% arrange(surveillance_category) %>% pull(seqid)

annotated_row <- rowAnnotation(df = rowA_stap)
annotated_col <- HeatmapAnnotation(df = colA_stap)

row_order_hm <-  stap_vfdb %>% arrange(surveillance_category) %>% pull(seqid)
stap_m <- stap_m[row_order_hm,col_order_hm]

p <-
  Heatmap(stap_m,column_title = "Staphylococcus aureus",col = c("gray80","brown"),
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

jpeg(filename = "Figures/vf_stap.jpg", units = "in", res = 300, height = 6,width = 10)
draw(p + annotated_row,annotation_legend_side = "top")
dev.off()

### working on a.bau
# find the column to discards
ind_abau <- which((abau_vfdb %>% map(~sum(. ==".")) %>% unlist() ) == nrow(abau_vfdb)) %>% as.vector()
abau_vfdb <- abau_vfdb %>% select(-all_of(ind_abau))
# subset virulence factor only 
abau_vf <- abau_vfdb[,c(3:142)]
abau_vf[abau_vf == "."] <- F
abau_vf[abau_vf != F ] <- T
abau_work <- abau_vfdb %>% select(seqid) %>% cbind(.,abau_vf)
row.names(abau_work) <- abau_work$seqid
abau_m <- abau_work %>% select(-seqid,-Scheme,-surveillance_category,-title) %>% select(-7,-10) %>% as.matrix()
# column annotation 
target_vf_abau <- colnames(abau_work)[-1]
vf_abau <- vfdb_detail %>% filter(Gene %in% target_vf_abau) %>% select(Gene,VFcategory) %>% 
  distinct(.keep_all = T) %>% arrange(VFcategory)  %>% slice(-84,-106)
col_order_hm <- vf_abau %>% pull(Gene)

colA_abau <- data.frame(Group = vf_abau$VFcategory)
rownames(colA_abau) <- vf_abau %>% pull(Gene)

# create a data for row annotation
rowA_abau <- data.frame(source = abau_vfdb$surveillance_category) %>% arrange(source)
rownames(rowA_abau) <- abau_vfdb %>% arrange(surveillance_category) %>% pull(seqid)

annotated_row <- rowAnnotation(df = rowA_abau)
annotated_col <- HeatmapAnnotation(df = colA_abau)

row_order_hm <-  abau_vfdb %>% arrange(surveillance_category) %>% pull(seqid)
abau_m <- abau_m[row_order_hm,col_order_hm]

# p <-
  Heatmap(abau_m,column_title = "Acinobacter baum",col = c("gray80","brown"),
          # row_dend_reorder = F,
          row_order = row_order_hm,
          column_order = col_order_hm,
          clustering_method_rows = "ward.D",
          column_labels = rep("",(n = ncol(abau_m))),
          top_annotation = annotated_col,
          # clustering_distance_columns =  "euclidean",
          heatmap_legend_param = list(
            title = "Genotype", labels = c("No","Yes")
          )
  )

jpeg(filename = "Figures/vf_kleb.jpg", units = "in", res = 300, height = 6,width = 10)
draw(p + annotated_row,annotation_legend_side = "top")
dev.off()





