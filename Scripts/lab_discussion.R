dorado_summary <- read.table("Data/summary.tsv", header = T)


library(ggplot2)
library(ggpubr)
p1 <- ggdensity(summary,x = "sequence_length_template", fill = "steelblue",
          add = "mean", rug = TRUE, xlab = "Sequence Length")

p2 <- ggdensity(summary,x = "mean_qscore_template", fill = "steelblue",
          add = "mean", rug = TRUE, xlab = "QScore")
jpeg("Data/Figure1.jpg",units = "in",res = 300, height = 6, width = 8)
ggarrange(p1,p2, ncol = 1)
dev.off()

guppy_basecalling <- read.table("Data/barcoding_summary.txt", header = T, sep = "\t")

guppy_basecall2 <- read.table("Data/sequencing_summary.txt",header = T, sep = "\t")
