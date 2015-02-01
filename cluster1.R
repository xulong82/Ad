# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project - clustering
# Rev: Dec 30, 2014

library(ape)
library(amap)
library(lattice)
library(ggplot2)
library(ggdendro)
library(pheatmap)

rm(list = ls())

setwd("~/Dropbox/AD")
load(file = "./R/complete_tpm.rdt")
col.manual <- c("grey70", "firebrick1", "chartreuse3", "dodgerblue3", "gold1", "darkorchid2")

# --- 2014 samples only
brain.tpm <- brain.tpm[, c(grep("2014", colnames(brain.tpm)), grep("mouse", colnames(brain.tpm)))]
brain.tpm <- brain.tpm[apply(brain.tpm, 1, function (x) max(x) > 5), ]
brain.tpm <- log2(brain.tpm + 1)
write.table(brain.tpm, file = "./Cluster3/brain2014complete.txt", sep = "\t", row.names = T, col.names = T, quote = F)

mycol <- rep("chartreuse3", ncol(brain.tpm))
mycol[grep("mouse", colnames(brain.tpm))] <- "firebrick"
hc1 <- hcluster(t(brain.tpm), method = "pearson", link = "average")  # clust samples
pdf("./Graphs/brain2014dendro.pdf", width = 12)
plot(as.phylo(hc1), type = "unrooted", tip.color = mycol, cex = .5, font = 2, lab4ut = "axial")
dev.off()

# --- mouse 1559
brain.tpm1 <- brain.tpm[apply(brain.tpm, 1, function(x) max(x) - min(x) > 1), ]
pdf(file = "~/Dropbox/AD/Graphs/mouse1559.pdf", height = 6, width = 9)
par(mfrow = c(2, 1), mar = c(2, 4, 1, 2))
plot(brain.tpm1$WT6m1484.2014 - brain.tpm1$APP5m1558.2014, ylim = c(-6, 6), xlab = "", ylab = "")
title("WT6m1484 - APP5m1558", line = -2, col.main = "firebrick1")
abline(0, 0, lwd = 1, col = "firebrick1")
plot(brain.tpm1$APP5m1559.2014 - brain.tpm1$APP5m2633.2014, ylim = c(-6, 6), xlab = "", ylab = "")
title("APP5m1559 - APP5m2633", line = -2, col.main = "firebrick1" )
abline(0, 0, lwd = 1, col = "firebrick1")
dev.off()

brain.tpm <- brain.tpm[, -grep("APP5m1559.2014", colnames(brain.tpm))]
save(brain.tpm, file = "~/Dropbox/AD/R/brain2014.rdt")

# --------------------------------------------------------------------
rm(list = ls())
load("./R/brain2014.rdt")
dt <- brain.tpm

# --- batch correction
batch <- rep("2014", ncol(dt))
batch[grep("mouse", colnames(dt))] <- "2014-new"

dt.mean <- apply(dt, 1, mean)
dt.center <- dt - dt.mean
dt.res <- t(apply(dt.center, 1, function (x) {lm(x ~ as.factor(batch) - 1)$res}))
dt <- dt.res + dt.mean

save(dt, file = "./R/batch2014.rdt")

mycol <- rep(col.manual[2], ncol(dt))
mycol[grep("4m", colnames(dt))] <- col.manual[3]
mycol[grep("5m", colnames(dt))] <- col.manual[4]
mycol[grep("6m", colnames(dt))] <- col.manual[5]

dt <- dt[apply(dt, 1, function(x) max(x) - min(x) > 1), ] 
hc1 <- hcluster(t(dt), method = "pearson", link = "average")  # clust samples
pdf("./Graphs/batch2014dendro.pdf")
plot(as.phylo(hc1), type = "unrooted", tip.color = mycol, cex = .5, font = 2, lab4ut = "axial")
dev.off()

load("./R/batch2014.rdt")
# --- ANOVA --- ANOVA's return confusing
treat <- gsub("^.*(WT|APP).*", "\\1", colnames(dt))
month <- gsub("^.*(2m|4m|5m|6m).*", "\\1", colnames(dt))
aov.pval <- NULL
for (i in 1:nrow(dt)) {
  if (i %% 1e3 == 0) cat(i, nrow(dt), "\n")
  aov1 = aov(as.matrix(dt)[i, ] ~ month * treat)
  aov.pval <- c(aov.pval, min(summary(aov1)[[1]][["Pr(>F)"]], na.rm = T))
}
dt.aov <- dt[aov.pval < 0.05, ]  # Select genes

dt.fc <- dt[apply(dt, 1, function(x) max(x) - min(x) > 1), ]  # FC works better

dt <- dt.aov
dt <- dt.fc

# ---
dt <- dt[, grep("6m", colnames(dt))]
hc1 <- hcluster(t(dt), method = "pearson", link = "average")  # clust samples
plot(as.phylo(hc1), type = "unrooted", cex = .5, font = 2, lab4ut = "axial")
# ---

complete2m <- dt[, grep("2m", colnames(dt))]
complete456m <- dt[, -grep("2m", colnames(dt))]

hc1 <- hcluster(t(complete456m), method = "pearson", link = "average")
plot(as.phylo(hc1), type = "unrooted", cex = .5, font = 2, lab4ut = "axial")

wt2m <- rowMeans(complete2m[, grep("WT", colnames(complete2m))])
app2m <- rowMeans(complete2m[, grep("APP", colnames(complete2m))])
wt456m <- complete456m[, grep("WT", colnames(complete456m))]
app456m <- complete456m[, grep("APP", colnames(complete456m))]
wt456m <- wt456m - wt2m
app456m <- app456m - app2m

dt <- cbind(wt456m, app456m)
col.manual <- c("firebrick1", "chartreuse3", "dodgerblue3", "gold1", "darkorchid2")
hc1 <- hcluster(t(dt), method = "pearson", link = "average")
clusts = cutree(hc1, 4)
pdf(file = "~/Dropbox/AD/Graphs/cluster1norm2dendro.pdf")
plot(as.phylo(hc1), type = "unrooted", tip.color = col.manual[clusts], cex = .5, font = 2, lab4ut = "axial")
dev.off()

save(dt, file = "~/Dropbox/AD/R/cluster1.rdt")
