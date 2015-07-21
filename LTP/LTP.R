library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggvis)
library(KEGGREST)
library(xlsx)

rm(list = ls())
setwd("~/Dropbox/GitHub/Ad")
load("data/brain2014.rdt") 

glm.dt <- brain.tpm
cutoff <- quantile(c(as.matrix(glm.dt)), 0.25)  # TPM level
glm.dt <- glm.dt[apply(glm.dt, 1, function(x) max(x) > cutoff & sum(x > 0) > round(ncol(glm.dt) / 10)), ]
sample <- colnames(glm.dt)
age <- factor(gsub("^.*(2m|4m|5m|6m).*", "\\1", sample), levels = c("2m", "4m", "5m", "6m"))
geno <- factor(gsub("^.*(WT|APP).*", "\\1", sample), levels = c("WT", "APP"))
batch <- factor(gsub("^.*(2014|mouse).*", "\\1", sample), levels = c("2014", "mouse"))
uid <- paste(age, geno, sep = "_")
grp <- intersect(c("2m_WT", "2m_APP", "4m_WT", "4m_APP", "5m_WT", "5m_APP", "6m_WT", "6m_APP"), uid)
nid <- factor(uid, levels = grp)
glm.dt_mean <- sapply(grp, function(x) rowMeans(glm.dt[, uid == x]))

# Model the data with GLM 
glm.fit <- apply(glm.dt, 1, function (x) summary(lm(x ~ age + geno + batch + age * geno)))
fit.batch <- sapply(glm.fit, function (x) x$coefficients["batchmouse", "Estimate"])
glm.dt_bc <- glm.dt - as.matrix(fit.batch) %*% (as.numeric(batch) - 1)
glm.fit <- apply(glm.dt_bc, 1, function (x) lm(x ~ age + geno + age * geno))
fit.r2 <- sapply(glm.fit, function (x) summary(x)$r.squared)
fit.fs <- sapply(glm.fit, function (x) summary(x)$fstatistic)
fit.pval <- apply(fit.fs, 2, function (x) pf(x[1], x[2], x[3], lower.tail = F)) 

glm.fit <- glm.fit[fit.pval < 0.05 & fit.r2 > 0.3]
fit.eff <- do.call(rbind, lapply(glm.fit, function (x) summary(x)$coefficients[, "Estimate"]))
fit.pval <- do.call(rbind, lapply(glm.fit, function (x) summary(x)$coefficients[, "Pr(>|t|)"]))

# Identify signal: ANOVA
aov.fit <- lapply(glm.fit, anova)
aovGene <- lapply(c("age", "geno", "age:geno"), function(x1) {
  y1 = sapply(aov.fit, function(x2) x2[x1, "Pr(>F)"]); y1[y1 < 0.05]
}); names(aovGene) <- c("age", "app", "age:app")

# Identify signal: GLM estiamtes
binary <- apply(fit.pval, 2, function (x) x < 0.05) & apply(fit.eff, 2, function (x) abs(x) > 0.1)
binary <- binary[apply(binary[, -1], 1, any), -1]
signal <- apply(binary, 1, function (x) paste(x, collapse = "-"))
signalTable <- sort(table(signal), decreasing = T)
(signalTable <- signalTable[signalTable > 30])
signalGene <- sapply(names(signalTable), function(x) names(signal)[signal == x])

op <- par(mar = c(5, 20, 4, 10))
bar <- barplot(signalTable, xlim = c(0, 750), axes = F, border = NA, horiz = T, las = 1, space = 0.75, cex.name = 0.5)
abline(v = 0, lwd = 1, col = "black")
text(y = bar, x = signalTable + 30, labels = signalTable)

binaryApp <- binary[, grep("APP", colnames(binary))]
binaryApp <- binaryApp[apply(binaryApp, 1, any), ]
signalApp <- apply(binaryApp, 1, function (x) paste(x, collapse = "-"))
signalTableApp <- sort(table(signalApp), decreasing = T)
(signalTableApp <- signalTableApp[signalTableApp > 30])
signalGeneApp <- sapply(names(signalTableApp), function(x) names(signalApp)[signalApp == x])

keggFind("pathway", "Long-term potentiation")  # KEGG: LTR
mmu04720 <- keggGet("mmu04720")[[1]]$GENE
mmu04720 <- gsub(";.*", "", matrix(mmu04720, byrow = T, ncol = 2)[, 2])

lapply(mmu04720, function(x) summary(glm.fit[[x]]))
lapply(signalGene, function(x) intersect(x, mmu04720))
lapply(signalGeneApp, function(x) intersect(x, mmu04720))

write(keggGet("mmu04720", "kgml"), file = "mmu04720.xml")  # Cytoscape
writePNG(keggGet("mmu04720", "image"), "mmu04720.png")

ggvis1 <- function(x) # Visualize expression
  data.frame(sample, uid, age, geno, batch, nid) %>% mutate(value = c(as.matrix(glm.dt[x, ]))) %>%
  ggvis(~as.numeric(nid), ~value) %>% layer_boxplots(fill=~nid, width = 0.5)

bgMs <- gsub("[,|;].*", "", keggList("mmu"))
myhyper <- function(g1, g2) {  # Hypergeometric
  if(length(intersect(g1, g2)) == 0) return(1)
  1 - phyper(length(intersect(g1, g2)) - 1, length(g2), length(setdiff(bgMs, g2)), length(g1))
}  # Pr(count >= length(intersect(g1, g2)))

intersect(rownames(binaryApp), mmu04720)
myhyper(rownames(binaryApp), mmu04720)
lapply(signalGene, function(x) myhyper(x, mmu04720))
lapply(signalGeneApp, function(x) myhyper(x, mmu04720))

xx <- glm.dt_mean[signalGene[[3]], grep("WT", colnames(glm.dt_mean))]
xx <- glm.dt_mean[signalGene[[3]], grep("APP", colnames(glm.dt_mean))]
xx <- log2(xx + 1) - log2(xx[, 1] + 1)

down = rownames(xx)[fit.eff[rownames(xx), "age5m"] < 0]
xx <- xx[down, ]

plot(1:4, rep(0, 4), type = "l", col = "blue", xlab = "", xaxt = "n", ylab = "", ylim = c(-0.1, 0.02))
lapply(1:nrow(xx), function(i) lines(xx[i, ], type = "l", col = "grey30"))
lines(colMeans(xx[, ]), type = "b", lwd = 3, col = "red")
axis(1, at=1:4, labels = paste0(c(2, 4, 5, 6), "m"))

source("../../X/kegg.R")
source("../../X/function.R")
mykegg <- kegg(rownames(xx))
mygk <- myGK(down)

ggvis1("Gnaq")
plot(glm.dt_mean["Gnaq", grep("WT", colnames(glm.dt_mean))], type = "b")
lines(glm.dt_mean["Gnaq", grep("APP", colnames(glm.dt_mean))], type = "b")

write.xlsx(glm.dt_mean[down, ], file = "LTP/ltp.xlsx", append = T)
write.xlsx(mygk[[2]], file = "LTR/ltp.xlsx", sheetName = "KEGG", append = T)
write.xlsx(mygk[[1]][1], file = "LTP/ltp.xlsx", sheetName = "BP", append = T)
write.xlsx(mygk[[1]][2], file = "LTP/ltp.xlsx", sheetName = "MF", append = T)
write.xlsx(mygk[[1]][3], file = "LTP/ltp.xlsx", sheetName = "CC", append = T)
write.xlsx(mmu04720, file = "LTP/ltp.xlsx", sheetName = "LTP", append = T)
