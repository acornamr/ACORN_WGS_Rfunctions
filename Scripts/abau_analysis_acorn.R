# this script is used to analyse Acineobacter bauma and sub species of ACORN project
# Written and maintained by Tung Trinh
# ver 1.0 
# For more information, please contact to tungts@oucru.org
####################################################################################################
##### required packages and built-in functions #####
library(dplyr)
library(stringr)
library(ggtree)
library(treeio)
library(ape)
library(readxl)
library(ggnewscale)
library(qdapTools)
library(ggplot2)
##### load the data #####
### ML tree 
tree <- read.tree("Data/iqtree/abau_full_ml_iqtree.contree")
# drop the reference
tree <- root(tree, outgroup = "Reference",edgelabel = T)
tree <- drop.tip(tree,"Reference")
### staramr 
staramr <- read_excel("Data/assemblies/Hanoi/out_acine/results.xlsx", sheet = 1)
### metadata 
meta <- read_excel("Data/ACORN2_VN001_WGS_Cumulative_V2_AKP.xlsx", sheet = 1 )
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
# create seqid
vfdb <- vfdb %>% mutate(seqid = str_extract(FILE,"A[0-9]-[0-9]+")) 
vfdb_acine <- vfdb %>% filter(seqid %in% id)
vfdb_acine_id <- vfdb_acine$seqid
vfdb_acine <- vfdb_acine[,-1]
vfdb_acine[vfdb_acine != "."] <- "Yes"
excluded_cols<- which(vfdb_acine %>% map(~unique(.) %>% length) == 1) %>% names
vfdb_acine <- vfdb_acine %>% select(!excluded_cols) %>% as.matrix()
row.names(vfdb_acine) <- vfdb_acine_id
## process and visualise the data 
# create the id 
staramr <- staramr %>% mutate(seqid = str_extract(`Isolate ID`,"A[0-9]-[0-9]+")) 
# get the id from the tree 
target_id <- tree$tip.label
# tree$tip.label <- gsub("-","_",tree$tip.label)
# filter the data based on the target_id 
acine_staramr <- staramr %>% filter(seqid %in% target_id) 
acine_meta <- meta %>% filter(`SEQ-ID` %in% target_id)
acine_vfdb <- vfdb %>% filter(seqid %in% target_id)
# select data from metadata to include into the tree
acine_meta <- acine_meta %>% select(label = `SEQ-ID`,group = surveillance_category,d28_status,surveillance_diag)
acine_meta <- acine_staramr %>% select(label = seqid,st = `Sequence Type`) %>% 
   left_join(acine_meta,., by = "label") 
# %>% mutate(label = gsub("-","_",label))
# join with the tree 
tree_added <- full_join(tree,acine_meta, by = "label")

### heatmap 
sub_hm1 <- meta %>% filter(`SEQ-ID` %in% target_id) %>% 
  mutate(clinical_condition = case_when(
    clinical_severity_score== 0 ~ "None",
    clinical_severity_score== 1 ~ "Mild",
    clinical_severity_score== 2 ~ "Moderate",
    clinical_severity_score== 3 ~ "Severe"
  )) %>% select(clinical_condition) %>% as.data.frame()

row.names(sub_hm1) <- acine_meta$label 

sub_hm2 <- meta %>% filter(`SEQ-ID` %in% target_id) %>% 
  select(d28_status) %>% as.data.frame()
row.names(sub_hm2) <- acine_meta$label 

sub_hm3 <- meta %>% filter(`SEQ-ID` %in% target_id) %>% 
  select(surveillance_diag) %>% as.data.frame()

row.names(sub_hm3) <- acine_meta$label 

p <- ggtree(tree_added,aes(x,y))+
  geom_tippoint(aes(colour = group),  size = 2)+
  geom_tiplab(aes(label = st), align = T, offset = 0.01,linetype = 3, size = 3, colour = "gray50")+
  scale_colour_manual(values = c("#1F8F99FF","#D60C00FF","#6B58EFFF"),
                      label = c("Community-acquired Infection (CAI)",
                                "Hospital-acquired Infection (HAI)",
                                "Healthcare-associated Infection (HCAI)"), name ="")+
  theme_tree()+theme(legend.position = "top")+guides(fill=guide_legend(ncol=2))

p1 <- gheatmap(p,sub_hm1, colnames = F,width = 0.1,color = NULL, offset = 0.006)+
  scale_fill_manual(values = c("#0571b0","#f4a582","#92c5de","#ca0020"), name = "Clinical Condition")+
  theme_tree()+theme(legend.position = "top")+guides(fill=guide_legend(ncol=2))

p1 <- p1 + new_scale_fill()


p2 <- gheatmap(p1, sub_hm2, width = 0.1, colnames = F, offset = 0.0425,color = NULL)+
  scale_fill_viridis_d(option = "A", name = "Outcome", direction = -1,na.value = NA,drop = T)+
  theme_tree()+theme(legend.position = "top")+guides(fill=guide_legend(ncol=2))

p2 <- p2 + new_scale_fill()
p3 <- gheatmap(p2,sub_hm3, width = 0.1, colnames = F, offset = 0.08,color = NULL)+
  scale_fill_manual(values = c("#7fc97f",
                               "#beaed4",
                               "#fdc086",
                               "#ffff99",
                               "#386cb0",
                               "#f0027f",
                               "#bf5b17",
                               "#666666"), 
                    name = "Diagnosis", na.value = NA,drop = T)+
  theme_tree()+theme(legend.position = "top")+guides(fill=guide_legend(ncol=2))
p3 <- p3+new_scale_fill()
###################
# transformed amr group 
# acine_amr_transformed <- read_excel("Data/acorn_acine_transformed_amr_2025-04-24.xlsx",sheet = 1)
acine_amr_transformed <- read_excel("Data/acorn_acine_transformed_amr_2025-04-24.xlsx",sheet = 2)
id<- acine_amr_transformed$seqid
acine_amr_transformed <- acine_amr_transformed[,-1]
acine_amr_transformed <- acine_amr_transformed %>% mutate_all(as.character)
acine_amr_transformed[acine_amr_transformed == "1"] <- "Yes"
acine_amr_transformed[acine_amr_transformed == "0"] <- "No"
acine_amr_transformed <- as.matrix(acine_amr_transformed)
row.names(acine_amr_transformed) <- id
# colnames(acine_amr_transformed) <- word(colnames(acine_amr_transformed),1)
colnames(acine_amr_transformed) <- word(colnames(acine_amr_transformed),-1)
jpeg("Figures/acine_tree_amr_group.jpeg",height = 10,width = 18,res = 300, unit = "in")
# jpeg("Figures/acine_tree_hydrolysis_group.jpeg",height = 10,width = 18,res = 300, unit = "in")
gheatmap(p3,acine_amr_transformed, colnames = T, offset = 0.13,color = NULL,colnames_angle=-90,font.size=2,hjust =0)+
  scale_fill_manual(values = c("gray80","brown"), name = "Genotype")+theme(legend.position = "top")+vexpand(.05,-1)
dev.off()

jpeg("Figures/acine_tree_vfdb.jpeg",height = 10,width = 18,res = 300, unit = "in")
gheatmap(p3,vfdb_acine, colnames = T, offset = 0.13,color = NULL,colnames_angle=-90,font.size=2,hjust =0)+
  scale_fill_manual(values = c("gray80","steelblue"), name = "Genotype")+theme(legend.position = "top")+vexpand(.05,-1)
dev.off()

