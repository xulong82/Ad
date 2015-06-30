library(dplyr)
library(ggplot2)

rm(list = ls())
setwd("~/Dropbox/GitHub/Load/")

trim <- read.delim("shell/trim.log", header = F)
rsem <- read.delim("shell/rsem.log", header = F)
trim <- gsub("^.*\\((.*)%\\) Forward.*", "\\1", trim$V1) %>% as.numeric
rsem <- gsub("^.*\\((.*)%\\)", "\\1", rsem$V1) %>% as.numeric

n.2013 <- nrow(rsem2013)
n.2014 <- nrow(rsem2014)
n.14new <- nrow(rsem14new)

total <- sum(n.2013, n.2014, n.14new)

gdt <- data.frame(value = c(trim, rsem), process = rep(c("Trim", "Bowtie"), each = 36))
gdt$process <- relevel(gdt$process, ref = "Trim")

pdf("pdf/qc.pdf", width = 3, height = 4)
ggplot(gdt, aes(x = process, y = value)) +
  geom_boxplot(aes(fill = process)) + 
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() + xlab("") + ylab("Percentage (%)") +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'grey30'),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, vjust = 2),
        legend.position = "none")
dev.off()
