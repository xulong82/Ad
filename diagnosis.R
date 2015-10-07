# Diagnosis

library(xlsx)
library(reshape)
library(gplots)
library(ggplot2)
library(genefilter)
library(dplyr)

options(stringsAsFactors = F)

rm(list = ls())
setwd("~/Dropbox/GitHub/Load/")

load("data/tpmAll.rdt")

Bin1 = tpmAll[grep("^B", names(tpmAll))]
Bin1 = Bin1[rowMax(as.matrix(Bin1)) > 10, ]
Bin1 = log2(Bin1 + 1)
geno = factor(gsub("-.*", "", names(Bin1)), levels = c("B6", "Bin1"))

gene = "Bin1"

mutate(melt(Bin1[gene, ]), group = gsub("-.*", "", variable)) %>%
  ggplot(aes(group, value, label = variable)) + geom_boxplot(width = 0.5) +
  theme_bw() + xlab("") + ylab("") + geom_text()

Bin1_remove = within(Bin1, "Bin1-6962" <- NULL)
geno_remove = factor(gsub("-.*", "", names(Bin1_remove)), levels = c("B6", "Bin1"))

ttest = rowttests(as.matrix(Bin1), geno)
ttest_remove = rowttests(as.matrix(Bin1_remove), geno_remove)

sum(ttest$p.value < 0.05)
sum(ttest_remove$p.value < 0.05)

listA = rownames(ttest)[ttest$p.value < 0.05]
listB = rownames(ttest_remove)[ttest_remove$p.value < 0.05]

report_1 = read.xlsx("data/Bin1.xlsx", sheetIndex = 1)
report_2 = read.xlsx("data/vsB6.xlsx", sheetName = "Bin1")

listC = report_1$symbol
listD = report_2$symbol

venn(list(A = listA, C = listC))
venn(list(B = listB, D = listD))

venn(list(listA, listB, listC))

Bin1[setdiff(listA, listC), ]
ttest[setdiff(listA, listC), ]
ttest_remove[setdiff(listD, listB), ]

tpmAll_remove = within(tpmAll, "Bin1-6962" <- NULL)
grp1 <- c("B6", "ApoE4", "Apoe", "Bin1", "Cd2ap", "Clu") 
grp2 <- factor(gsub("-.*", "", names(tpmAll_remove)), levels = grp1)

tpmAll_remove = tpmAll_remove[rowMax(as.matrix(tpmAll_remove)) > 10, ]

y = as.matrix(log2(tpmAll_remove + 1))
hist(y[, 1])

fit <- apply(y, 1, function (x) lm(x ~ grp2))

fit.pv <- sapply(fit, function(x) summary(x)$coefficients[-1, "Pr(>|t|)"]) %>% t
fit.et <- sapply(fit, function(x) summary(x)$coefficients[-1, "Estimate"]) %>% t

geneList = rownames(fit.pv)[fit.pv[, "grp2Bin1"] < 0.05]
gk = mmGK(geneList)

summary(lm(y[gene, ] ~ grp2))

gene = "Ubr3"

join = intersect(rownames(fit.pv), rownames(ttest_remove))
join = cbind(fit.pv[join, ], ttest_remove[join, ])

summary(lm(y[gene, ] ~ grp2))
summary(lm(as.matrix(Bin1_remove)[gene, ] ~ geno_remove))

join[gene, ]

Bin1 = tpmAll[grep("^B", names(tpmAll))]
yy = Bin1[geneList, ]
yy$"Bin1-6962" = NULL
yy = cbind(yy, B6_avg = rowMeans(yy[, grep("B6", names(yy))]))
yy = cbind(yy, Bin1_avg = rowMeans(yy[, grep("Bin1", names(yy))]))
yy$log2fc = log2(yy$Bin1_avg + 1) - log2(yy$B6_avg + 1)
yy$fit.pv = fit.pv[geneList, "grp2Bin1"]

write.xlsx(yy, file = "temp.xlsx", sheetName = "new_2", append = T)
write.xlsx(gk$GO$CC, file = "temp.xlsx", sheetName = "CC", append = T)
write.xlsx(gk$KEGG, file = "temp.xlsx", sheetName = "KEGG", append = T)
write.xlsx(gk$GO$BP, file = "temp.xlsx", sheetName = "BP", append = T)
write.xlsx(gk$GO$MF, file = "temp.xlsx", sheetName = "MF", append = T)
write.xlsx(gk$GO$CC, file = "temp.xlsx", sheetName = "CC", append = T)
