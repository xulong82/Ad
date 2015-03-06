# Xulong Wang (xulong.wang@jax.org)

library(ape)
library(amap)
library(ggdendro)

rm(list = ls())
setwd("~/Dropbox/GitHub/Ad")
load("data/brain_bc2014.rdt")
load("data/glm_brain.rdt")
source("function.R")

geneId.app = glm$less$app$symbol
geneId.age = glm$less$age$symbol
dt.hc <- dt.bc[unique(c(geneId.app, geneId.age)), ]
dim(dt.hc)

hc1 <- hcluster(t(dt.hc), method = "pearson", link = "average")
hc2 <- hcluster(dt.hc, method = "pearson", link = "average")

mycol <- rep("grey50", ncol(dt.hc))
mycol[grep("APP", colnames(dt.hc))] <- "firebrick1"
plot(as.phylo(hc1), edge.width = 2, font = 2, cex = 0.7, label.offset = 1e-4, tip.color = mycol, direction = "downward")

group1<- c("APP5m1558.2014", "APP5m2751.2014", "APP5m1633.2014","APP5m1636.2014","APP6m1684.2014","APP5m1648.2014", "APP5m1623.2014")
group2<- c("mouse_3346_6m_APP", "mouse_3440_6m_APP", "mouse_3351_6m_APP", "mouse_6585_6m_APP", "APP5m1647.2014", "APP5m1738.2014")

mycol[colnames(dt.hc) %in% c(group1, group2)] <- "blue"
plot(as.phylo(hc1), edge.width = 2, font = 2, cex = 0.7, label.offset = 1e-4, tip.color = mycol, direction = "downward")
myhc <- as.phylo(hc1)

dt.g1 <- dt.bc[, group1]
dt.g2 <- dt.bc[, group2]
colnames(dt.g1) <- paste("group1", colnames(dt.g1), sep = "_")
colnames(dt.g2) <- paste("group2", colnames(dt.g2), sep = "_")
dt.g12 <- cbind(dt.g1, dt.g2)
dim(dt.g12)

fc <- rowMeans(dt.g1) - rowMeans(dt.g2) 
treat <- gsub("^.*(group1|group2).*", "\\1", colnames(dt.g12))
tt.pval <- apply(dt.g12, 1, function(x) pairwise.t.test(x, treat, p.adj = "none")$p.value)
tt.qval <- p.adjust(tt.pval, method = "fdr")

idx <- which(tt.qval < 0.05 & abs(fc) > 0.2)
geneId <- rownames(dt.g12)[idx]
head(geneId)

gk.five <- myGK(geneId)
head(gk.five$BP[, c("Term", "Pvalue")])
head(gk.five$MF[, c("Term", "Pvalue")])
head(gk.five$CC[, c("Term", "Pvalue")])
head(gk.five$KEGG[, c("Term", "Pvalue")])

five <- list()
five$gdt <- myhc
five$gdt$mycol <- mycol
five$symbol <- geneId
five$gk <- gk.five
save(five, file = "data/five.rdt")
