library(dplyr)
library(quantro)
library(preprocessCore)
library(genefilter)
library(VennDiagram)
library(png)
library(scales)
library(tidyr)

rm(list = ls())
setwd("~/Dropbox/GitHub/Load/")
source("../../X/kegg.R")

load("data/myList.rdt")
for(obj in names(myList)) assign(obj, myList[[obj]])

grp1 <- colnames(ge)
grp2 <- c("ApoE4", "Apoe", "Bin1", "Cd2ap", "Clu") 
grp3 <- c("B6", "ApoE4", "Apoe", "Bin1", "Cd2ap", "Clu") 
grp4 <- factor(gsub("-.*", "", grp1), levels = grp3)

ge <- log2(ge + 1)

pdf("pdf/boxplot.pdf", width = 8, height = 5)
par(mar = c(8, 4, 4, 2))
matboxplot(ge, groupFactor = grp4)
dev.off()

# Average median of distributions
objectMedians <- apply(ge, 2, median)
anovaFit <- anova(lm(objectMedians ~ grp4))
anovaPval <- (anovaFit$`Pr(>F)`[1] < 0.05)

objectNorm <- sweep(ge, 2, objectMedians, FUN = "-")
Fnik = apply(objectNorm, 2, sort) 

Fndotk = sapply( grp3, function(x){ 
  rowMeans(Fnik[, which(grp4 %in% x)]) 
} )

Fndotdot = rowMeans(Fnik)

(betweenDiff = colMeans( (Fndotk - replicate(6, Fndotdot))^2 ))
boxplot(betweenDiff)
MSb = sum(betweenDiff * table(grp4)) / (6 - 1)

withinDiff <- do.call(cbind, lapply(grp3, function(x) Fnik[, which(grp4 %in% x) ] - Fndotk[, x]))
withinDiff <- colMeans((withinDiff)^2)
boxplot(withinDiff ~ gsub("-.*", "", names(withinDiff)))
MSe <- sum(withinDiff) / (length(grp4) - length(grp3))	

(quantroStat <- MSb / MSe)
        
(qtest <- quantro(object = ge, groupFactor = grp4, B = 1e3))  # quantro package

pdf("pdf/quantro.pdf", width = 8, height = 5)
quantroPlot(qtest)
dev.off()

ge[rowMax(ge) > 1e4, ]  # suspicious

vsB6 <- ge[, grep("B6|Apoe", grp1)]
vsB6 <- ge[, grep("B6|ApoE", grp1)]
vsB6 <- ge[, grep("B6|Bin1", grp1)]
vsB6 <- ge[, grep("B6|Cd2ap", grp1)]
vsB6 <- ge[, grep("B6|Clu", grp1)]

vsB6.norm <- normalize.quantiles(vsB6, copy = TRUE)
dimnames(vsB6.norm) <- dimnames(vsB6)
group <- factor(gsub("-.*", "", colnames(vsB6)), levels = c("B6", "Apoe"))

pdf("pdf/qqplot.pdf", width = 8, height = 5)
par(mfrow = c(1, 2))
qqplot(vsB6[, 1], vsB6[, 7], xlab = "Apoe", ylab = "B6")
abline(0, 1)
qqplot(vsB6.norm[, 1], vsB6.norm[, 7], xlab = "Apoe", ylab = "B6")
abline(0, 1)
dev.off()

ttest.raw = rowttests(vsB6, group)
ttest.norm = rowttests(vsB6.norm, group)

ttest.raw$sig <- with(ttest.raw, abs(dm) > 0.1 & p.value < 0.05)
ttest.norm$sig <- with(ttest.norm, abs(dm) > 0.1 & p.value < 0.05)

inter <- kegg(rownames(vsB6)[intersect(which(ttest.raw$sig), which(ttest.norm$sig))])
diff1 <- kegg(rownames(vsB6)[setdiff(which(ttest.raw$sig), which(ttest.norm$sig))])
diff2 <- kegg(rownames(vsB6)[setdiff(which(ttest.norm$sig), which(ttest.raw$sig))])

vennList <- list(raw = which(ttest.raw$sig), norm = which(ttest.norm$sig))
venn.diagram(vennList, cex = 2, fill = c("grey70", "dodgerblue3"), imagetype = "png", file = "temp/temp.png")

pdf("pdf/vocano2.pdf", width = 8, height = 5)
# ggplot(ttest.raw, aes(x = dm, y = -log10(p.value))) + 
ggplot(ttest.norm, aes(x = dm, y = -log10(p.value))) + 
  geom_point(aes(color = as.factor(sig))) +
  scale_color_manual(values = c("grey30", "firebrick1")) +
  theme_bw() + xlab("Effect") + ylab("-log10(pvalue)") +
  theme(axis.line = element_line(color = "grey30"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = "none")
dev.off()

# PCA: mathematical interpration of biological system
fit <- apply(ge, 1, function(x) lm(x ~ grp4))
anova <- lapply(fit, anova)
anova.pval <- sapply(anova, function(x) x["grp4", "Pr(>F)"])

pca.ge <- ge[anova.pval < 0.05, ]

var.ge <- rowVars(ge)
pca.ge <- ge[var.ge > median(var.ge), ]

# pca.ge <- sapply(grp3, function(x) rowMeans(pca.ge[, grp4 == x]))
# pca.ge <- pca.ge - rowMeans(pca.ge)

pca.ge <- pca.ge - rowMeans(pca.ge[, grep("B6", colnames(pca.ge))])
pca.ge <- cbind(pca.ge[, grepl("B6", grp4)], pca.ge[, ! grepl("B6", grp4)])

pca <- prcomp(pca.ge)

pdf("pdf/var.pdf", width = 8, height = 5)
barplot(pca$sdev / sum(pca$sdev), col = c(rep("red", 6), rep("grey", 29)), 
        xlab = "PC", ylab = "% variants explained")
dev.off()

rotation <- -pca$rotation
rownames(rotation) = gsub("-.*", "", rownames(rotation))
group <- factor(rownames(rotation), levels = grp3)

pdf("pdf/pcs.pdf", width = 10, height = 7)
par(mfrow = c(1, 2), mar = c(10, 4, 4, 2))
lapply(1:12, function(x) {
  barplot(rotation[, x], col = as.numeric(group), xaxt = "n")
  if (x %in% 1:4) legend("topleft", levels(grp4), col = seq_along(group), pch = 19, cex = 0.7)
}); dev.off()

pc4 <- pca$x[, 4]
hist(pca$x[, 4])
pc4 <- pc4[abs(pc4) > quantile(abs(pc4), 0.95)]

dat <- sapply(grp3, function(x) rowMeans(ge[names(pc4), grp4 == x])) 
dat <- apply(dat, 1, scale) %>% t %>% as.data.frame
colnames(dat) <- grp3

dat1 <- dat[with(dat, ApoE4 > Apoe), ]
dat2 <- dat[with(dat, ApoE4 < Apoe), ]

kegg1 <- kegg(rownames(dat1))
kegg2 <- kegg(rownames(dat2))

pdf("pdf/line2.pdf", width = 8, height = 5)
plot(colMeans(dat1), type = "b", lwd = 5, col = "red", ylim = c(-3, 3), xaxt = "n", xlab = "", ylab = "")
axis(1, at = 1:6, labels = levels(grp4))
apply(dat1, 1, function(x) lines(x, col = alpha("grey", 0.5)))
dev.off()

dat <- mutate(dat, gene = rownames(dat))
dat <- gather(dat, group, value, B6:Clu)

barplot(svd$v[, 1], col = as.numeric(grp4))
barplot(svd$v[, 2], col = as.numeric(grp4))

svd <- svd(pca.ge)
z <- svd$d * t(svd$v)
plot(z[1, ], z[2, ], col = as.numeric(grp4))
legend("topleft", levels(grp4), col = seq_along(grp4), pch=1)
