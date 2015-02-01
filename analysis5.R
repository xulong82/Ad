# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: AD project - Clustering II
# Rev: June 14, 2014

library(lattice)
library(ggplot2)
library(ggdendro)
library(amap)
library(ape)

rm(list = ls())
load("~/Dropbox/AD/R/brain.tpm1.rdt")  # brain TPM without batch correction
load("~/Dropbox/AD/R/brain.tpm2.rdt")  # brain TPM with batch correction
#-------------
tpm1 <- brain.tpm2
id1 <- which(apply(tpm1, 1, function(x) {min(x) < 0}))  # negative TPM from correction
tpm1[id1, ] <- tpm1[id1, ] + abs(apply(tpm1[id1, ], 1, min))
tpm2 <- log2(tpm1 + 1)  # log2 transformation

#--- Gene selection: ANOVA
tpm3 <- tpm2[, -grep("5m", colnames(tpm2))]  # 2, 4, 6 month
group <- gsub("m.*", "m", colnames(tpm3))
month <- gsub("WT", "", gsub("APP", "", group))
treat <- gsub("[2456]m", "", group)

aov.pval <- rep(NA, nrow(tpm3))
for (i in 1:nrow(tpm3)) {
  if (i %% 1e3 == 0) cat(i, "in", nrow(tpm3), "\n")
  dat1 <- data.frame(tpm = as.matrix(tpm3)[i, ], month, treat)
  aov1 = aov(tpm ~ month, data = dat1)
  aov.pval[i] <- min(summary(aov1)[[1]][["Pr(>F)"]], na.rm = T)
}

tpm4 <- tpm3[aov.pval < 0.05, ]  # Select genes

#--- Nearest shrunken centroid
cs.id <- c("2m", "4m", "6m")
cs.n1 <- length(cs.id)  # cluster number
cs.n2 <- rep(NA, cs.n1)  # sample number per cluster
cs.mean1 <- rowMeans(tpm4)  # total mean
cs.mean2 <- matrix(NA, nrow = nrow(tpm4), ncol = cs.n1, dimnames = list(rownames(tpm4), cs.id))
cs.var <- matrix(NA, nrow = nrow(tpm4), ncol = cs.n1, dimnames = list(rownames(tpm4), cs.id))

for (i in 1:cs.n1) {
  cat(i, "in", cs.n1, "\n")
  month = cs.id[i]
  dat1 <- tpm4[, grep(month, colnames(tpm4))]
  
  cs.n2[i] <- ncol(dat1)
  cs.mean2[, i] <- rowMeans(dat1)
  cs.var[, i] <- (apply(dat1, 1, sd))^2
}

cs.mk <- sqrt(1/cs.n2 + 1/sum(cs.n2))
cs.si <- sqrt(rowSums(cs.var %*% (cs.n2 - 1)) / (sum(cs.n2) - cs.n1))  # pooled sd
cs.so <- median(cs.si)

dik1 <- t(apply(cs.mean2 - cs.mean1, 1, function(x) {x / cs.mk})) / cs.si

shrink <- 0
shrink <- 1
shrink <- 2
shrink <- 3
shrink <- 4
shrink <- 5

dik2 <- apply(dik1, 2, function(x) {ifelse((abs(x) - shrink) > 0, sign(x) * (abs(x) - shrink), 0)})
centroid1 <- t(apply(cs.si * dik2, 1, function (x) {x * cs.mk}))  # shrink = 0
centroid2 <- t(apply(cs.si * dik2, 1, function (x) {x * cs.mk}))  # shrink = 5
centroid3 <- cs.mean2 + centroid2

cs.id1 <- rownames(centroid3)[which(apply(centroid2, 1, function (x) {max(abs(x)) > 0}))]
centroid3 <- centroid3[cs.id1, ]
sample.al <- tpm4[cs.id1, ]  # 2,4,6 month
sample.5m <- tpm2[cs.id1, grep("5m", colnames(tpm2))]  # 5 month
dat1 <- cbind(centroid3, sample.al)
dat2 <- cbind(centroid3, sample.5m)

dat1.cor <- cor(dat1, method = "pearson")[1:3, -c(1:3)]
dat1.dis <- as.matrix(dist(t(dat1), method = "euclidean"))[1:3, -c(1:3)]
dat2.cor <- cor(dat2, method = "pearson")[1:3, -c(1:3)]
dat2.dis <- as.matrix(dist(t(dat2), method = "euclidean"))[1:3, -c(1:3)]
pdf("~/Dropbox/AD/Figures/cluss2.pdf", height = 8, width = 12)
levelplot(t(dat2.dis), main = "", xlab   = "", ylab = "",
          scales = list(x = list(rot = 90), y = list()))
dev.off()

pdf("~/Dropbox/AD/Figures/nsc2.pdf")
par(mfrow = c(1, cs.n1), mar = c(5, 4, 4, 2))
for (i in 1:cs.n1) {
  name <- paste(paste("Cluster", i, sep = " "), cs.id[i], sep = ": ")
  plot(centroid2[, i], 1:nrow(centroid2), 
       type = "l", col = "red", 
       main = name, xlab = "Averaged expression", ylab = "Gene") 
}
dev.off()

tpm5 <- tpm2[, ]
tpm5 <- tpm4[cs.id1, ]
tpm5 <- tpm5 - apply(tpm5, 1, mean)
hc1 <- hcluster(t(tpm5), method = "pearson", link = "average")  # clust samples
hc2 <- hcluster(tpm5, method = "pearson", link = "average")  # clust genes

colors = rainbow(7, s = 1, v = .7)[c(1, 3, 5)]
clusts = cutree(hc1, 3)
pdf("~/Dropbox/AD/Figures/dendro1.pdf", height = 4)
# plot(as.phylo(hc1), type = "unrooted", tip.color = colors[clusts], cex = .5, font = 2, lab4ut = "axial")
plot(as.phylo(hc1), tip.color = colors[clusts], direction = "downwards", cex = .5, font = 2)
dev.off()

colors = rainbow(5, s = 1, v = .7)
clusts = cutree(hc2, 5)
pdf("~/Dropbox/AD/Figures/dendro4.pdf", width = 15)
plot(as.phylo(hc2), direction = "upwards", tip.color = colors[clusts], cex = .5, font = 2)
dev.off()

names(clusts[clusts == 1])
names(clusts[clusts == 2])
names(clusts[clusts == 3])
names(clusts[clusts == 4])
names(clusts[clusts == 5])
for (i in 1:5) {
  write.table(names(clusts[clusts == i]), row.names = FALSE, col.names = FALSE, quote = FALSE, 
              file = paste(paste("~/Dropbox/AD/List/gene.cs", i, sep = ""), "txt", sep = "."))
}

tile2.dat1 <- tpm5[hc2$order, hc1$order]
tile2.dat2 <- t(apply(tile2.dat1, 1, rank))
tile2.dat3 <- data.frame(value = c(tile2.dat2), 
                         gene = rep(rownames(tile2.dat2), time = ncol(tile2.dat2)),
                         sample = rep(colnames(tile2.dat2), each = nrow(tile2.dat2)))
tile2.dat3$sample <- factor(tile2.dat3$sample, levels = colnames(tile2.dat1))
tile2.dat3$gene <- factor(tile2.dat3$gene, levels = rownames(tile2.dat1))
pdf("~/Dropbox/AD/Figures/tile2.pdf", width = 15, height = 5)
ggplot(tile2.dat3, aes(x = gene, y = sample, fill = value)) + 
  geom_tile() + scale_fill_gradient(low="green", high="red") +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5),
        axis.ticks = element_blank()) +
  theme(legend.position = "none")
dev.off()

#-------
