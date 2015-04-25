# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project - general QC
# Rev: Dec 16, 2014

library(pheatmap)
library(ggplot2)

rm(list = ls())

load("~/Dropbox/AD/R/data1314.rdt")  # tpm per sample
load("~/Dropbox/AD/R/data14new.rdt")  # tpm per sample

table(brain2013.tpm$sym == brain2014.tpm$sym)
table(brain2014.tpm$sym == brain14new.tpm$sym)

# --- build data frame 
brain.tpm <- data.frame(row.names = brain2013.tpm$sym, brain2013.tpm[, -1], brain2014.tpm[, -1], brain14new.tpm[, -1])

retina.tpm <- data.frame(row.names = retina2014.tpm$sym, retina2014.tpm[, -1], retina14new.tpm[, -1])
retina.tpm <- retina.tpm[, -grep("1484L", colnames(retina.tpm))]  # trim-bowtie.R

brain.tpm  <- brain.tpm[apply(brain.tpm, 1, function (x) {max(x) > 5}), ]  # maximal
retina.tpm <- retina.tpm[apply(retina.tpm, 1, function (x) {max(x) > 10}), ]  # maximal

#--- median correlation
cor_retina <- cor(retina.tpm, method = "spearman")
heatmap(cor_retina)

cor_brain <- cor(brain.tpm, method = "spearman")
heatmap(cor_brain)

median_cor_brain <- apply(cor_brain, 1, function(x) {median(sort(x, decreasing = T)[-1])})
median_cor_retina <- apply(cor_retina, 1, function(x) {median(sort(x, decreasing = T)[-1])})

median_cor.df <- data.frame(value = c(median_cor_brain, median_cor_retina),
  tissue = c(rep("Brain", length(median_cor_brain)), rep("Retina", length(median_cor_retina))))
pdf("~/Dropbox/AD/QC/median_cor.pdf", width = 5)
ggplot(median_cor.df, aes(x = value, fill = tissue)) +
  geom_histogram(binwidth = 0.01) + facet_grid(. ~ tissue) +
# scale_fill_manual(values = c("#1f78b4", "#fb9a99")) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() + 
  xlab("Median Correlation Coefficient (Spearman)") + ylab("Count") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold", vjust = -0.5),
        axis.title.y = element_text(size = 12, face = "bold", vjust = 1)) +
  theme(legend.position = "none")
dev.off()

# tpm1.dat1 <- brain.tpm1[, grep("APP6m", colnames(brain.tpm1))]  # brain
# tpm1.dat1 <- retina.tpm1[, grep("APP2m", colnames(retina.tpm1))]  # retina
# pdf("~/Dropbox/AD/Figures/pair2.pdf")
# pairs(log2(tpm1.dat1 + 1))
# dev.off()

brain.tpm <- brain.tpm[, -grep("226941", colnames(brain.tpm))]  # delete the outlier sample
retina.tpm <- retina.tpm[, -grep("3310", colnames(retina.tpm))]  # delete the outlier sample

# --- App, Prnp
spikein <- c("App", "Prnp")
spikein.brain <- brain.tpm[spikein, ]
spikein.retina <- retina.tpm[spikein, ]

type <- rep("WT", time = ncol(spikein.brain))
type[grep("APP", colnames(spikein.brain))] <- "APP"
spikein.df <- data.frame(value = c(as.matrix(spikein.brain)), 
  sample = rep(colnames(spikein.brain), each = 2),
  gene = rep(c("App", "Prnp"), ncol(spikein.brain)),
  type = factor(rep(type, each = 2), levels = c("WT", "APP")))

type <- rep("WT", time = ncol(spikein.retina))
type[grep("APP", colnames(spikein.retina))] <- "APP"
spikein.df <- data.frame(value = c(as.matrix(spikein.retina)), 
  sample = rep(colnames(spikein.retina), each = 2),
  gene = rep(c("App", "Prnp"), ncol(spikein.retina)),
  type = factor(rep(type, each = 2), levels = c("WT", "APP")))

pdf(file = "~/Dropbox/AD/QC/spikein-brain.pdf", width = 10)
pdf(file = "~/Dropbox/AD/QC/spikein-retina.pdf", width = 3)
ggplot(spikein.df, aes(x = sample, y = value, fill = type)) + 
  geom_bar(width = 0.75, stat = "identity") + facet_grid(gene ~ .) + theme_bw() +
  scale_fill_manual(values = c("blue", "red")) +
  xlab("") + ylab("TPM") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 6, angle = -90, face = "bold"),
        axis.title = element_text(size = 10, face = "bold", vjust = 1),
        strip.text = element_text(size = 10, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

brain.tpm <- brain.tpm[, -grep("265182", colnames(brain.tpm))]

save(brain.tpm, retina.tpm, file = "~/Dropbox/AD/R/complete_tpm.rdt")
