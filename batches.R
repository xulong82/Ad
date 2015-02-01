# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project - batch effects
# Rev: Dec 30, 2014

library(pheatmap)
library(ggplot2)

rm(list = ls())
load(file = "~/Dropbox/AD/R/complete_tpm.rdt")

# --- optional: 2014 samples only
brain.tpm <- brain.tpm[, c(grep("2014", colnames(brain.tpm)), grep("mouse", colnames(brain.tpm)))]

# --- PCA
brain.tpm_ts <- log2(brain.tpm + 1)
brain.tpm_ts <- brain.tpm_ts - apply(brain.tpm_ts, 1, mean)
pca_brain <- prcomp(brain.tpm_ts)

barplot(pca_brain$rotation[, 1])

genotype <- rep("WT", ncol(brain.tpm))
genotype[grep("APP", colnames(brain.tpm))] <- "APP"
graph.dt <- data.frame(value = c(pca_brain$rotation[, 1:3]), 
                       pc = rep(paste("PC", 1:3, sep = ""), each = ncol(brain.tpm)), 
                       sample = rep(colnames(brain.tpm), 3),
                       type = rep(genotype, 3))
graph.dt$order <- reorder(colnames(brain.tpm), 1:ncol(brain.tpm))

# col.manual <- c("grey70", "firebrick1", "chartreuse3", "dodgerblue3", "gold1", "darkorchid2")
pdf("~/Dropbox/AD/Figures/pca_brain.pdf", height = 7, width = 10)
pdf("~/Dropbox/AD/Figures/pca_brain2014.pdf", height = 5, width = 10)
ggplot(graph.dt, aes(x = order, y = value, fill = type)) + 
  geom_bar(stat = "identity", width = .75) + facet_grid(pc ~ .) +
  scale_fill_manual(values = c("firebrick1", "chartreuse3")) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  xlab("") + ylab("") +
  theme(panel.border = element_rect(size = 1, color = "grey30")) +
  theme(axis.text.x = element_text(size = 6, angle = -90, face = "bold", color = "grey30"),
        axis.text.y = element_text(size = 6, face = "bold"),
        axis.title = element_text(size = 8, face = "bold"),
        strip.text = element_text(size = 8, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 8, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank())
dev.off()
