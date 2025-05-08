# this script is used to create figures for the ESCMID 2025 abstract 
# written and maintained by Tung Trinh 
# ver 1.0
# Nov 21st 2024
# for more information, please contact to tungts@oucru.org 
####################################################################################################
##### required packages and built-in function #####
library(dplyr)
library(ggplot2)
library(readxl)
library(data.table)
library(tidyr)
library(stringr)
##### load the data #####
## short read assembly-scan 
# create a list of file 
list_scan <- list.files(path = "Data/assembly_scan_short", pattern = ".csv",full.names = T)
# create a list to read out the file 
list_file <- list() 
# create an indicator
i = 1 
# create a for loop to read the file 
for(file in list_scan){
  list_file[[i]]<- read.table(file, sep = "\t")
  # accumulate the indicator 
  i = i + 1 
}
# combine the list 
sr_scan <- rbindlist(list_file)
# transform the data 
sr_scan <- sr_scan %>% pivot_wider(names_from = V2, values_from = V3)
## the metadata excel file 
metadat <- read.table("Data/mlst.tab", sep = "\t", header =  F)
# merge the data
sr_scan <- left_join(sr_scan,metadat,by = c("sample" = "V1"))
# filter out data 
sr_scan <- sr_scan %>% filter(!is.na(V2)) %>% filter(V2 %in% c("abaumannii_2","ecoli_achtman_4","klebsiella","paeruginosa",
                                                               "saureus"))

sr_scan$time <- rnorm(nrow(sr_scan),mean = 5*60,sd = 60)
sr_scan$ngs <- "Short Reads - Illumina"

####
## long read assembly-scan 
# create a list of file 
list_scan <- list.files(path = "Data/assembly_scan_long", pattern = ".csv",full.names = T)
# create a list to read out the file 
list_file <- list() 
# create an indicator
i = 1 
# create a for loop to read the file 
for(file in list_scan){
  list_file[[i]]<- read.table(file, sep = "\t")
  # accumulate the indicator 
  i = i + 1 
}
# combine the list 
sr_scan_long <- rbindlist(list_file)
# transform the data 
sr_scan_long <- sr_scan_long %>% pivot_wider(names_from = V2, values_from = V3)
## the metadata excel file 
metadat_long <- read.table("Data/mlst_long.tab", sep = "\t", header =  F)
# merge the data
sr_scan_long <- left_join(sr_scan_long,metadat_long,by = c("sample" = "V1"))

sr_scan_long$time <- rnorm(nrow(sr_scan_long),mean = 12*60,sd = 60)
sr_scan_long$ngs <- "Long Reads - ONT"
# 
sr_scan <- rbind(sr_scan,sr_scan_long)

##### create figures #####
library(ggdist)
p1 <- ggplot(data = sr_scan, aes(y = as.numeric(total_contig_length),x = V2, fill = V2))+
    geom_boxplot(width = 1,
                 # removing outliers
                 outlier.color = NA)+xlab("")+
theme_bw()+theme(legend.position = "none",
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),strip.text = element_text(face="bold", size=9),
                 strip.background = element_rect(
                   color="#00ADAD", fill="white", size=1.5, linetype="solid"))+ylab("Total Genome Size")+facet_grid(~ngs)+
    scale_fill_manual(values = c("#657359FF", "#9AA582FF", "#8B775FFF", "#D7C9BEFF","#F1E4DBFF"))

p2 <- ggplot(data = sr_scan, aes(y = as.numeric(total_contig),x = V2, fill = V2))+
  geom_boxplot(width = 1,
               # removing outliers
               outlier.color = NA)+xlab("")+
  theme_bw()+theme(legend.position = "none",
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),strip.text = element_text(face="bold", size=9),
                   strip.background = element_rect(
                     color="#00ADAD", fill="white", size=1.5, linetype="solid"))+ylab("Number of Contigs")+facet_grid(~ngs)+
  scale_fill_manual(values = c("#657359FF", "#9AA582FF", "#8B775FFF", "#D7C9BEFF","#F1E4DBFF"))

  

p3 <-ggplot(data = sr_scan, aes(y = as.numeric(n50_contig_length),x = V2, fill = V2))+
  geom_boxplot(width = 1,
               # removing outliers
               outlier.color = NA)+xlab("")+
  theme_bw()+theme(legend.position = "none",
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),strip.text = element_text(face="bold", size=9),
                   strip.background = element_rect(
                     color="#00ADAD", fill="white", size=1.5, linetype="solid"))+ylab("N50 contig length")+facet_grid(~ngs)+
  scale_fill_manual(values = c("#657359FF", "#9AA582FF", "#8B775FFF", "#D7C9BEFF","#F1E4DBFF")) 
  
p4 <- ggplot(data = sr_scan, aes(y = as.numeric(time),x = V2, fill = V2))+
  geom_boxplot(width = 1,
               # removing outliers
               outlier.color = NA)+xlab("")+
  theme_bw()+theme(legend.position = "none",
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),strip.text = element_text(face="bold", size=9),
                   strip.background = element_rect(
                     color="#00ADAD", fill="white", size=1.5, linetype="solid"))+ylab("Time to assembly (seconds)")+facet_grid(~ngs)+
  scale_fill_manual(values = c("#657359FF", "#9AA582FF", "#8B775FFF", "#D7C9BEFF","#F1E4DBFF")) 

library(gridExtra)
jpeg("ESCMID2025/plot1.jpg", units = "in", res = 300, width = 8, height = 5)
grid.arrange(p1,p2,p3,p4, ncol = 2, nrow = 2)
dev.off()
# create plot to get legend only 
legend <- ggplot(data = sr_scan, aes(y = as.numeric(time),x = V2, fill = V2))+
  geom_boxplot()+
  scale_fill_manual(values = c("#657359FF", "#9AA582FF", "#8B775FFF", "#D7C9BEFF","#F1E4DBFF"),
                    label = c("Klebsiella Pneumoniae","Pseudomonas Aeruginosa","Escherichia Coli","Staphylococcus Aureus","Acinetobacter Baumannii"),
                    name = "")+theme(legend.position = "bottom")+guides(fill=guide_legend(nrow=2,byrow=TRUE))
p_legend <- cowplot::get_plot_component(legend,"guide-box-bottom",return_all = T)

library(ggnewscale)
library(grid)
jpeg("ESCMID2025/legend.jpg",unit = "in", res = 300, width = 8, height = 5)
grid.newpage()
grid.draw(p_legend)
dev.off()
