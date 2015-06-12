# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project - Quality control
# Rev: Dec 15, 2014

library(ggplot2)

trim2013 <- read.delim("~/Dropbox/AD/QC/trim2013.csv", header = F)
trim2014.brain <- read.delim("~/Dropbox/AD/QC/trim2014.brain.csv", header = F)
trim2014.retina <- read.delim("~/Dropbox/AD/QC/trim2014.retina.csv", header = F)
trim2014 <- rbind(trim2014.brain, trim2014.retina)

rsem2013 <- read.delim("~/Dropbox/AD/QC/rsem2013.csv", header = F)
rsem2014.brain <- read.delim("~/Dropbox/AD/QC/rsem2014.brain.csv", header = F)
rsem2014.retina <- read.delim("~/Dropbox/AD/QC/rsem2014.retina.csv", header = F)
rsem2014 <- rbind(rsem2014.brain, rsem2014.retina)

trim14new <- read.delim("~/Dropbox/AD/QC/trim14new.csv", header = F)
rsem14new <- read.delim("~/Dropbox/AD/QC/rsem14new.csv", header = F)

n.2013 <- nrow(rsem2013)
n.2014 <- nrow(rsem2014)
n.14new <- nrow(rsem14new)

total <- sum(n.2013, n.2014, n.14new)

boxplot.dat <- data.frame(
  value = c(trim2013$V1, trim2014$V1, trim14new$V1, rsem2013$V1, rsem2014$V1, rsem14new$V1),
  year = c(rep(2013, n.2013 * 2), rep(2014, n.2014 * 2), rep("2014-new", n.14new * 2),
           rep(2013, n.2013), rep(2014, n.2014), rep("2014-new", n.14new)),
  process = factor(c(rep("Trim", total * 2), rep("Bowtie", total)), levels = c("Trim", "Bowtie"))
)

pdf("~/Dropbox/AD/QC/trim-bowtie.pdf")
ggplot(boxplot.dat, aes(x = process, y = value)) +
  geom_boxplot(aes(fill = process)) + facet_grid(. ~ year) +
# scale_fill_manual(values = c("#1f78b4", "#fb9a99")) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() + 
  xlab("") + ylab("Percentage (%)") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold", vjust = 1)) +
  theme(legend.position = "none")
dev.off()
