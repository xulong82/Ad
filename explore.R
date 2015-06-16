library(dplyr)
library(quantro)
library(preprocessCore)
library(VennDiagram)
library(png)

rm(list = ls())
setwd("~/Dropbox/GitHub/Load/")

load("data/myList.rdt")
for(obj in names(myList)) assign(obj, myList[[obj]])

grp1 <- colnames(ge)
grp2 <- c("ApoE4", "Apoe", "Bin1", "Cd2ap", "Clu") 
grp3 <- c("B6", "ApoE4", "Apoe", "Bin1", "Cd2ap", "Clu") 
grp4 <- factor(gsub("-.*", "", grp1), levels = grp3)

# Greg in PCA: mathematical interpration of biological system
var.ge <- apply(ge, 1, var)
pca.ge <- ge[var.ge > median(var.ge), ]
pca.ge <-  pca.ge - rowMeans(pca.ge[, grep("B6", colnames(pca.ge))])
pca <- prcomp(pca.ge)
barplot(pca$sdev)
par(mfrow = c(1, 2))
barplot(pca$rotation[, 1], las = 2)
barplot(pca$rotation[, 2], las = 2)
plot(pca$rotation[, 1], pca$rotation[, 2], xlim = c(-0.6, 1.0))

# Rafa in normalization
(qtest <- quantro(object = ge, groupFactor = grp4))
qtest <- quantro(object = ge, groupFactor = grp4, B = 1e3)
matboxplot(ge, groupFactor = grp4)
matboxplot(log2(ge + 1), groupFactor = grp4)

ge[rowMax(ge) > 1e4, ]  # suspicious

vsB6 <- ge[, grep("B6|Apoe", grp1)]
vsB6 <- ge[, grep("B6|ApoE", grp1)]
vsB6 <- ge[, grep("B6|Bin1", grp1)]
vsB6 <- ge[, grep("B6|Cd2ap", grp1)]
vsB6 <- ge[, grep("B6|Clu", grp1)]

vsB6.norm <- normalize.quantiles(vsB6, copy = TRUE)
dimnames(vsB6.norm) <- dimnames(vsB6)

matboxplot(log2(vsB6 + 1), factor(gsub("-.*", "", colnames(vsB6))))
matboxplot(log2(vsB6.norm + 1), factor(gsub("-.*", "", colnames(vsB6))))

(qtest <- quantro(vsB6, factor(gsub("-.*", "", colnames(vsB6))), B = 1e3))
quantroPlot(qtest)

(qtest <- quantro(vsB6.norm, factor(gsub("-.*", "", colnames(vsB6.norm))), B = 1e3))
quantroPlot(qtest)

vsB6_ttest <- function(vsB6) {  # GLM fit: vs B6
  grp <- relevel(factor(gsub("-.*", "", colnames(vsB6))), "B6")
  fit <- apply(log2(vsB6 + 1), 1, function (x) lm(x ~ grp))
  fit.pv <- sapply(fit, function(x) summary(x)$coefficients[2, "Pr(>|t|)"])
  fit.et <- sapply(fit, function(x) summary(x)$coefficients[2, "Estimate"])
  which(fit.pv < 0.05 & abs(fit.et) > 0.1) %>% names %>% return
}

raw = vsB6_ttest(vsB6)
norm = vsB6_ttest(vsB6.norm)

venn.diagram(list(raw = raw, norm = norm), imagetype = "png", file = "temp/temp.png")
example(readPNG)
