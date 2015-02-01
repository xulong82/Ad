# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project - B6xC3H
# Rev: June 4, 2014

library(VennDiagram)
library(ggplot2)
#-------------
rm(list = ls())
load("~/Dropbox/AD/R/brain.tpm1.rdt")  # brain TPM without batch correction
load("~/Dropbox/AD/R/brain.tpm2.rdt")  # brain TPM with batch correction
#-------------
tpm1 <- brain.tpm2
id1 <- which(apply(tpm1, 1, function(x) {min(x) < 0}))  # negative TPM from correction
tpm1[id1, ] <- tpm1[id1, ] + abs(apply(tpm1[id1, ], 1, min))
tpm2 <- log2(tpm1 + 1)  # log2 transformation

gids <- c("Arpp21", "Rrp9", "Manf", "Pfkfb4")
data1 <- tpm2[gids, ]
data2 <- data.frame(value = c(data1),
                    gene = rep(rownames(data1), time = ncol(data1)),
                    sample = rep(colnames(data1), each = nrow(data1)),
                    type = rep(gsub("[2456]m.*", "", colnames(data1)), each = nrow(data1)))
pdf("~/Dropbox/AD/B6xC3H/genes.chr9.pdf")
ggplot(data2, aes(x = sample, y = value, fill = type)) + 
  geom_bar(stat = "identity") + facet_grid(gene ~ .) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  xlab("") + ylab("") +
  theme(panel.border = element_rect(size = 2, color = "black")) +
  theme(axis.text.x = element_text(size = 7, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 7, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank())
dev.off()
